# --- Working directory & config -------------------------------------------
# Run from the GAMBLR.data checkout (holds inst/extdata and data/); override
# with the GAMBLR_DATA_ROOT env var. No hard-coded paths.
PKG_ROOT <- Sys.getenv("GAMBLR_DATA_ROOT", unset = getwd())
stopifnot(dir.exists(PKG_ROOT))
setwd(PKG_ROOT)

library(readxl)
library(GAMBLR)
library(parallel)
library(tidyverse)

# GAMBLR.results locates GSC data via config.yml (repo_base / project_base).
# check_config_and_value() looks first in the working directory, then falls back
# to the copy shipped in GAMBLR.results/inst/extdata. As a last resort (and to
# cover any direct config::get() calls), point R_CONFIG_FILE at that shipped
# config when nothing else provides one.
if (!file.exists("config.yml") && !nzchar(Sys.getenv("R_CONFIG_FILE"))) {
    Sys.setenv(R_CONFIG_FILE = system.file("extdata", "config.yml",
                                            package = "GAMBLR.results"))
}

# get_gambl_metadata()'s min_corrected_cov QC filter (default 15x) can drop
# samples that were previously part of a released bundle out of a fresh one
# (e.g. 01-16433_tumorA/B: FFPE genomes with coverage under the bar, silently
# excluded, taking all their SNVs with them). Disable it for bundle assembly
# so a rebuild doesn't lose samples that are already in the released dataset.
# Shadows the package function for the rest of this script, so every call
# site (there are ~13) picks this up in one place.
get_gambl_metadata <- function(...) {
    GAMBLR.results::get_gambl_metadata(..., min_corrected_cov = 0)
}

# Global variables definition
colnames_for_bundled_meta <- c(
    "patient_id",
    "sample_id",
    "Tumor_Sample_Barcode",
    "seq_type",
    "sex",
    "COO_consensus",
    "lymphgen",
    "genetic_subgroup",
    "EBV_status_inf",
    "cohort",
    "pathology",
    "reference_PMID"
)

pmids <- list(
    "Dreval_FL" = 37084389,
    "Grande_BL" = 30617194,
    "Thomas_BL" = 36201743,
    "Reddy_DLBCL" = 28985567,
    "Schmitz_DLBCL" = 29641966,
    "Chapuy_DLBCL" = 29713087,
    "Chapuy_other" = 22343534,
    "Arthur_DLBCL" = 30275490,
    "Hilton_DLBCL" = 37319384
)

maf_columns_to_keep <- c(
    "RefSeq",
    "Protein_position"
)

all_cols <- c(
    names(GAMBLR.helpers:::maf_header[1:45]),
    maf_columns_to_keep
)

# restrict to the most inclusive DLBCL gene list. Expanded to include known
# aliases (e.g. old/new HGNC histone names) so the Hugo_Symbol %in%
# all_lymphoma_genes filters below don't drop rows annotated under a gene's
# other name.
all_lymphoma_genes <- GAMBLR.utils::expand_gene_aliases(lymphoma_genes_comprehensive$Gene)

# Built once, reused by every get_ssm_by_regions() call below that wants only
# lymphoma-gene mutations. Restricting via tabix -R (region) instead of
# pulling every mutation for a sample set is the performance win; every call
# site still adds `filter(Hugo_Symbol %in% all_lymphoma_genes)` afterward on
# the now-small result, because gene_to_region()'s coordinates and VEP's
# Hugo_Symbol assignment don't always agree at gene boundaries -- without the
# post-filter, mutations from overlapping/neighbouring genes and non-coding
# loci (e.g. AC/AL/AF-prefixed lncRNA transcripts) leak in as false "gains".
# GENE_PAD_BP widens the tabix window itself so real target-gene mutations
# just outside gene_to_region()'s exact span (promoter/UTR/annotation-source
# discrepancies) aren't lost before the Hugo_Symbol filter even sees them;
# the filter makes over-padding cheap (extra I/O, not incorrect inclusion).
# 10kb is a 2x margin over VEP's default 5kb upstream/downstream annotation
# window (the likely source of Hugo_Symbol on these rows) -- e.g. an ID3
# variant 4.2kb upstream of its TSS (well inside VEP's 5kb default, but
# outside a 2kb pad) was confirmed lost under the old 2000bp value.
GENE_PAD_BP <- 10000
lymphoma_genes_bed_grch37 <- create_bed_data(
    gene_to_region(gene_symbol = all_lymphoma_genes, projection = "grch37",
                   return_as = "bed", pad_length = GENE_PAD_BP),
    genome_build = "grch37"
)
lymphoma_genes_bed_hg38 <- create_bed_data(
    gene_to_region(gene_symbol = all_lymphoma_genes, projection = "hg38",
                   return_as = "bed", pad_length = GENE_PAD_BP),
    genome_build = "hg38"
)

# Wraps an expression, printing its wall-clock time; returns the expression's
# value unchanged. Used below to benchmark get_ssm_by_samples()
# (subset_from_merge TRUE vs FALSE, per-sample loop) against
# get_ssm_by_regions() (tabix -R) at each call site, to decide which
# approach to standardize on for each use case.
time_it <- function(label, expr) {
    t0 <- Sys.time()
    result <- expr
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    message(sprintf("[TIMING] %-55s %8.1fs", label, elapsed))
    result
}

##### Importing SLMS-3 variants
# Reused by the single consolidated SLMS-3 pull in Phase 4 below (and by
# nothing else -- every cohort-specific SLMS-3 pull that used to call this
# once per cohort block has been removed; see Phase 4).
pull_data <- function(
        pull_meta,
        pull_projection = "grch37"
    ){
    lymphoma_genes_bed <- if(pull_projection == "grch37") lymphoma_genes_bed_grch37 else lymphoma_genes_bed_hg38
    slms3 <- get_ssm_by_regions(
        regions_bed = lymphoma_genes_bed,
        these_samples_metadata = pull_meta,
        basic_columns = FALSE,
        projection = pull_projection
    ) %>%
    filter(Hugo_Symbol %in% all_lymphoma_genes) %>%
    select(
        all_of(all_cols)
    )
    return(slms3)
}


# ============================================================================
# Phase 1: metadata (and, where genuinely distinct per-paper data exists,
# Publication-pipeline SSM calls) for every cohort. No SLMS-3 pulling here --
# that used to happen per-cohort-block, which silently duplicated mutation
# rows whenever a sample belonged to more than one cohort (confirmed: e.g.
# the Dreval FL x Hilton trios overlap). Every cohort below also builds a
# small (sample_id, study, study_id, reference_PMID) frame, assembled into
# the new sample_study table in Phase 2 -- this replaces the old inline
# `Study =` tagging on mutation rows entirely, for every pipeline, so cohort
# membership is tracked once per sample instead of once per mutation row.
# ============================================================================

# --- Thomas BL (grch37, genome) --------------------------------------------
# It has more patients and also contains sample ids, not just patient ids
bl_data <- list()

bl_data$meta <- read_xlsx("inst/extdata/studies/BL_Thomas.xlsx", sheet = 1)

bl_data$meta_to_bundle <- data.frame(
    bl_data$meta$`Patient barcode`,
    bl_data$meta$`Genome sample id`,
    bl_data$meta$`Genome sample id`,
    rep("genome", nrow(bl_data$meta)),
    bl_data$meta$Sex,
    rep(NA, nrow(bl_data$meta)),
    rep(NA, nrow(bl_data$meta)),
    rep(NA, nrow(bl_data$meta)),
    bl_data$meta$`EBV status`,
    rep("BL_Thomas", nrow(bl_data$meta)),
    rep("BL", nrow(bl_data$meta)),
    rep(pmids$Thomas_BL, nrow(bl_data$meta))
)

colnames(bl_data$meta_to_bundle) <- colnames_for_bundled_meta

bl_data$meta_to_bundle <- bl_data$meta_to_bundle %>%
    filter(! sample_id == "NA") %>%
    arrange(sample_id)

bl_data$meta_to_bundle <- read_xlsx(
        "inst/extdata/studies/BL_Thomas.xlsx",
        sheet = 12
    ) %>%
    rename("sample_id" = "Patient barcode") %>%
    filter(sample_id %in% bl_data$meta_to_bundle$sample_id) %>%
    select(sample_id, Subgroup) %>%
    right_join(
        bl_data$meta_to_bundle,
        .
    ) %>%
    mutate(genetic_subgroup = Subgroup) %>%
    select(- Subgroup)

bl_data$ssm_to_bundle <- read_xlsx(
        "inst/extdata/studies/BL_Thomas.xlsx",
        sheet = 6
    ) %>%
    rename(
        "Tumor_Sample_Barcode" = "tumor_biospecimen_id",
        "Matched_Norm_Sample_Barcode" = "normal_biospecimen_id"
    ) %>%
    select(names(GAMBLR.helpers:::maf_header[1:45])) %>%
    filter(Tumor_Sample_Barcode %in% bl_data$meta_to_bundle$sample_id)


bl_data$cnv_to_bundle <- read_xlsx(
        "inst/extdata/studies/BL_Thomas.xlsx",
        sheet = 4
    ) %>%
    select(-normal_biospecimen_id) %>%
    rename(
        "ID" = "tumor_biospecimen_id",
        "log.ratio" = "depth.ratio"
    ) %>%
    mutate(LOH_flag = NA, .after = end) %>%
    filter(ID %in% bl_data$meta_to_bundle$sample_id) %>%
	mutate(CN = round(2 * 2^log.ratio))

# study_id mirrors sample_id here (not patient_id): the xlsx's own
# "Genome sample id" column was already adopted directly as sample_id, and
# patient_id is patient-level, not sample-level -- ambiguous for any cohort
# where a patient could have more than one sample (see hilton_study below
# for a case where that ambiguity is real).
thomas_bl_study <- data.frame(
    sample_id = bl_data$meta_to_bundle$sample_id,
    study = "Thomas",
    study_id = as.character(bl_data$meta_to_bundle$sample_id),
    reference_PMID = pmids$Thomas_BL
)


# --- Dreval FL (grch37, genome) ---------------------------------------------
fl_data <- list()

fl_data$meta <- read_xlsx("inst/extdata/studies/FL_Dreval.xlsx", sheet = 1)

fl_data$meta_to_bundle <- data.frame(
    fl_data$meta$`Patient barcode`,
    fl_data$meta$`Genome sample id`,
    fl_data$meta$`Genome sample id`,
    rep("genome", nrow(fl_data$meta)),
    fl_data$meta$Sex,
    rep(NA, nrow(fl_data$meta)),
    rep(NA, nrow(fl_data$meta)),
    fl_data$meta$`cFL/dFL label`,
    rep(NA, nrow(fl_data$meta)),
    rep("FL_Dreval", nrow(fl_data$meta)),
    fl_data$meta$Pathology,
    rep(pmids$Dreval_FL, nrow(fl_data$meta))
)

colnames(fl_data$meta_to_bundle) <- colnames_for_bundled_meta

fl_data$meta_to_bundle <- fl_data$meta_to_bundle %>%
    filter(! sample_id == "NA") %>%
    arrange(sample_id)


fl_data$ssm_to_bundle <- read_xlsx(
        "inst/extdata/studies/FL_Dreval.xlsx",
        sheet = 2
    )

difference <- setdiff(
    names(GAMBLR.helpers:::maf_header[1:45]),
    colnames(fl_data$ssm_to_bundle)
)

new_cols <- setNames(rep(NA, length(difference)), difference)

fl_data$ssm_to_bundle  <- fl_data$ssm_to_bundle %>%
    mutate(!!! new_cols) %>%
    select(names(GAMBLR.helpers:::maf_header[1:45]))


fl_data$cnv_to_bundle <- read_xlsx(
        "inst/extdata/studies/FL_Dreval.xlsx",
        sheet = 3
    ) %>%
    rename(
        "ID" = "Tumor_Sample_Barcode"
    ) %>%
    mutate(CN = round(2 * 2^log.ratio))

# study_id mirrors sample_id, not patient_id -- see thomas_bl_study above.
dreval_study <- data.frame(
    sample_id = fl_data$meta_to_bundle$sample_id,
    study = "Dreval",
    study_id = as.character(fl_data$meta_to_bundle$sample_id),
    reference_PMID = pmids$Dreval_FL
)


# --- Thomas DLBCL (hg38, genome) --------------------------------------------
dlbcl_data <- list()

dlbcl_data$meta <- read_xlsx("inst/extdata/studies/BL_Thomas.xlsx", sheet = 2)

dlbcl_data$meta_to_bundle <- data.frame(
    dlbcl_data$meta$`Patient barcode`,
    dlbcl_data$meta$`Genome sample id`,
    dlbcl_data$meta$`Genome sample id`,
    rep("genome", nrow(dlbcl_data$meta)),
    dlbcl_data$meta$Sex,
    rep(NA, nrow(dlbcl_data$meta)),
    rep(NA, nrow(dlbcl_data$meta)),
    rep(NA, nrow(dlbcl_data$meta)),
    dlbcl_data$meta$`EBV status`,
    rep("DLBCL_Thomas", nrow(dlbcl_data$meta)),
    rep("DLBCL", nrow(dlbcl_data$meta)),
    rep(pmids$Thomas_BL, nrow(dlbcl_data$meta))
)

colnames(dlbcl_data$meta_to_bundle) <- colnames_for_bundled_meta

dlbcl_data$meta_to_bundle <- dlbcl_data$meta_to_bundle %>%
    filter(! sample_id == "NA") %>%
    arrange(sample_id)

dlbcl_data$ssm_to_bundle <- read_xlsx(
        "inst/extdata/studies/BL_Thomas.xlsx",
        sheet = 6
    ) %>%
    rename(
        "Tumor_Sample_Barcode" = "tumor_biospecimen_id",
        "Matched_Norm_Sample_Barcode" = "normal_biospecimen_id"
    ) %>%
    select(names(GAMBLR.helpers:::maf_header[1:45])) %>%
    filter(Tumor_Sample_Barcode %in% dlbcl_data$meta_to_bundle$sample_id)

dlbcl_data$cnv_to_bundle <- read_xlsx(
        "inst/extdata/studies/BL_Thomas.xlsx",
        sheet = 4
    ) %>%
    select(-normal_biospecimen_id) %>%
    rename(
        "ID" = "tumor_biospecimen_id",
        "log.ratio" = "depth.ratio"
    ) %>%
    mutate(LOH_flag = NA, .after = end) %>%
    filter(ID %in% dlbcl_data$meta_to_bundle$sample_id) %>%
	mutate(CN = round(2 * 2^log.ratio))

# study_id mirrors sample_id, not patient_id -- see thomas_bl_study above.
thomas_dlbcl_study <- data.frame(
    sample_id = dlbcl_data$meta_to_bundle$sample_id,
    study = "Thomas",
    study_id = as.character(dlbcl_data$meta_to_bundle$sample_id),
    reference_PMID = pmids$Thomas_BL
)


# --- DLBCL capture cohorts: Reddy, Schmitz, Chapuy, Golub -------------------
reddy_data <- list()
schmitz_data <- list()
chapuy_data <- list()
golub_data <- list()

# Importing metadata from Reddy et al and updating IDs to be consistent with GAMBL metadata
reddy_meta_full <- read_excel(
        "inst/extdata/studies/DLBCL_Reddy.xlsx",
        sheet = 1
    ) %>%
    mutate(
        study_id = `Sample  ID`,
        patient_id = paste0(
            "Reddy_",
            `Sample  ID`
        ),
        sample_id = paste0(
            "Reddy_",
            `Sample  ID`,
            "T"
        ),
        Tumor_Sample_Barcode = sample_id
    ) %>%
    dplyr::rename("COO_consensus" = "ABC GCB (RNAseq)") %>%
    mutate(
        COO_consensus = ifelse(
            COO_consensus == "Unclassified",
            "UNCLASS",
            COO_consensus
        )
    )

# Reddy's own paper-specific sample identifier ("Sample  ID"), captured
# before it's dropped by the select() below -- see reddy_study.
reddy_study <- data.frame(
    sample_id = reddy_meta_full$sample_id,
    study = "Reddy",
    study_id = as.character(reddy_meta_full$study_id),
    reference_PMID = pmids$Reddy_DLBCL
)

reddy_meta <- reddy_meta_full %>%
    dplyr::select(
        sample_id,
        patient_id,
        Tumor_Sample_Barcode,
        sex = Gender,
        COO_consensus
    )

# (config is resolved via config.yml / the GAMBLR.results fallback; no chdir to
# a gambl repo is needed to locate data — repo_base is absolute)

# patient_id is deliberately retained here (unlike the original version of
# this select(), which dropped it) so all-capture samples get a real
# patient_id in sample_data$meta instead of NA.
reddy_meta_gambl <- get_gambl_metadata() %>%
    dplyr::filter(cohort == "dlbcl_reddy") %>%
    dplyr::select(
        sample_id, patient_id, lymphgen, EBV_status_inf, cohort, pathology,
        seq_type, unix_group, genome_build, pairing_status, normal_sample_id
    ) %>%
    mutate(reference_PMID = pmids$Reddy_DLBCL)

reddy_data$meta_to_bundle <- left_join(
    reddy_meta,
    reddy_meta_gambl
)

schmitz_data$meta <- get_gambl_metadata() %>%
    dplyr::filter(cohort == "dlbcl_schmitz") %>%
    mutate(reference_PMID = pmids$Schmitz_DLBCL) %>%
    mutate(
        genetic_subgroup = lymphgen_wright,
        lymphgen = lymphgen_wright
    )

chapuy_data$meta <- get_gambl_metadata() %>%
    dplyr::filter(cohort == "dlbcl_chapuy") %>%
    mutate(reference_PMID = pmids$Chapuy_DLBCL) %>%
    mutate(genetic_subgroup = lymphgen) %>%
    filter(!sample_id == "DLBCL-RICOVER_148-Tumor")

golub_data$meta <- get_gambl_metadata() %>%
    dplyr::filter(cohort == "NCI_DLBCL_Golub") %>%
    mutate(reference_PMID = pmids$Chapuy_other) %>%
    mutate(genetic_subgroup = lymphgen)

all_capture_meta <- bind_rows(
        schmitz_data$meta,
        chapuy_data$meta,
        golub_data$meta
    ) %>%
    select(
        all_of(
            colnames(reddy_meta_gambl)
        )
    ) %>%
    bind_rows(
        .,
        reddy_meta_gambl
    )

# No paper-specific ID has been parsed anywhere in this script for these
# three cohorts (unlike Reddy/Thomas/Dreval/Arthur, they're pulled straight
# from GAMBL's own metadata, never cross-referenced against each paper's own
# supplementary sample table) -- study_id is NA until that sourcing is done.
schmitz_study <- data.frame(
    sample_id = schmitz_data$meta$sample_id,
    study = "Schmitz",
    study_id = NA_character_,
    reference_PMID = pmids$Schmitz_DLBCL
)
chapuy_study <- data.frame(
    sample_id = chapuy_data$meta$sample_id,
    study = "Chapuy",
    study_id = NA_character_,
    reference_PMID = pmids$Chapuy_DLBCL
)
golub_study <- data.frame(
    sample_id = golub_data$meta$sample_id,
    study = "NCI_Golub",
    study_id = NA_character_,
    reference_PMID = pmids$Chapuy_other
)

# Add data from Reddy paper's own original variant calls (distinct from the
# SLMS-3 recall above -- this is the as-published set, kept as its own
# Publication-pipeline pull, not deduplicated against SLMS-3).
reddy_original_maf <- read_tsv(
    "inst/extdata/studies/reddy_original_variants_with_VAF.maf.gz"
) %>%
    select(any_of(all_cols))


# --- DLBCL cell lines --------------------------------------------------------
cell_lines_data <- list()

cell_lines_data$meta <- get_gambl_metadata(seq_type_filter = "genome") %>%
    filter(sample_id %in% c(
        "DOHH-2", "SU-DHL-10", "OCI-Ly10", "OCI-Ly3", "SU-DHL-4"
    )) %>%
    arrange(sample_id)

cell_lines_data$meta_to_bundle <- data.frame(
    cell_lines_data$meta$patient_id,
    cell_lines_data$meta$sample_id,
    cell_lines_data$meta$sample_id,
    cell_lines_data$meta$seq_type,
    cell_lines_data$meta$sex,
    rep(NA, nrow(cell_lines_data$meta)),
    rep(NA, nrow(cell_lines_data$meta)),
    rep(NA, nrow(cell_lines_data$meta)),
    cell_lines_data$meta$EBV_status_inf,
    rep("DLBCL_cell_lines", nrow(cell_lines_data$meta)),
    rep("DLBCL", nrow(cell_lines_data$meta)),
    rep(NA, nrow(cell_lines_data$meta))
)

colnames(cell_lines_data$meta_to_bundle) <- colnames_for_bundled_meta

cell_lines_data$grch37$cnv_to_bundle <- get_cn_segments(
    these_samples_metadata = cell_lines_data$meta,
    projection="grch37"
) %>%
    dplyr::select(all_of(c("ID","chrom","start","end","LOH_flag","log.ratio","CN","seg_seq_type")))

cell_lines_data$hg38$cnv_to_bundle <- get_cn_segments(
    these_samples_metadata = cell_lines_data$meta,
    projection="hg38"
) %>%
    dplyr::select(all_of(c("ID","chrom","start","end","LOH_flag","log.ratio","CN","seg_seq_type")))

# No sample_study row for cell lines -- no published cohort applies to them
# (matches the previous Study=NA behaviour for this group).


# --- Arthur (grch37, genome) -------------------------------------------------
# Previously bundled from a raw, completely unfiltered flat file
# (inst/extdata/studies/DLBCL_Arthur.maf.gz, tagged Pipeline="strelka") with
# no sample-list gate and no lymphoma-gene-panel restriction at all --
# 2.83M rows / 65,121 distinct genes, ~18,480 "mutations" per sample vs
# ~137 for a properly gene-panel-restricted SLMS-3 sample. That block is
# removed entirely; Arthur's samples now flow through the same consolidated
# SLMS-3 pull as every other cohort (Phase 4).
# Case ID -> patient_id is done as an explicit join (not a
# `patient_id %in% arthur_case_ids$\`Case ID\`` filter) so a type mismatch
# between the xlsx's Case ID and GAMBL's patient_id (e.g. numeric vs
# character -- the same class of bug already fixed for study_id/Reddy's
# Sample ID elsewhere in this script; read_xlsx() can silently type a
# leading-zero ID like "08-15460" as numeric) can't silently drop matching
# patients from a %in% comparison. transmute() keeps only what's needed for
# the join, since nothing else from this sheet is used downstream.
arthur_case_ids <- read_xlsx(
    "inst/extdata/studies/DLBCL_Arthur.xlsx",
    sheet = 1
) %>% filter(`WGS data` == 1) %>%
    transmute(patient_id = as.character(`Case ID`))

# sample_id/Tumor_Sample_Barcode are GAMBL's own real values throughout --
# never overwritten to the paper's own patient-level ID (as the removed
# code used to do). Arthur's own "Case ID" is captured separately as
# study_id in arthur_study below, joined to metadata by sample_id like any
# other study-specific identifier.
#
# No `!grepl("tumor", sample_id)` filter here (the previous version of this
# block had one): it excluded every sample for any patient with more than
# one tumor biopsy (e.g. ..._tumorA/..._tumorB -- the same multi-sample-
# per-patient pattern Hilton has), dropping those patients' SLMS-3 coverage
# entirely. No other cohort in this script excludes samples this way.
arthur_meta <- get_gambl_metadata() %>%
    mutate(patient_id = as.character(patient_id)) %>%
    inner_join(arthur_case_ids, by = "patient_id") %>%
    filter(seq_type == "genome") %>%
    mutate(
        cohort = "DLBCL_Arthur",
        reference_PMID = pmids$Arthur_DLBCL
    )

arthur_study <- data.frame(
    sample_id = arthur_meta$sample_id,
    study = "Arthur",
    study_id = as.character(arthur_meta$patient_id),
    reference_PMID = pmids$Arthur_DLBCL
)


# --- Hilton trios (genome + capture, both projections) ----------------------
trios_samples <- read_xlsx(
    "inst/extdata/studies/DLBCL_Hilton.xlsx"
) %>%
drop_na(DNAseq_sample_id)

# (config is resolved via config.yml / the GAMBLR.results fallback; no chdir to
# a gambl repo is needed to locate data — repo_base is absolute)

trios_meta <- get_gambl_metadata() %>%
    filter(
        seq_type %in% c("genome", "capture"),
        sample_id %in% trios_samples$DNAseq_sample_id
    ) %>%
    mutate(
        cohort = "DLBCL_Hilton",
        reference_PMID = pmids$Hilton_DLBCL
    )

# study_id uses sample_id, not patient_id: Hilton is a trios study, so a
# single patient can have multiple samples (e.g. LY_RELY_116_tumorA and
# LY_RELY_116_tumorB) -- patient_id would collapse them to the same
# study_id, making the two rows ambiguous. sample_id already carries the
# distinguishing suffix and matches what the paper itself would call each
# sample.
hilton_study <- data.frame(
    sample_id = trios_meta$sample_id,
    study = "Hilton",
    study_id = as.character(trios_meta$sample_id),
    reference_PMID = pmids$Hilton_DLBCL
)


# ============================================================================
# Phase 2: assemble sample_data$meta and sample_data$sample_study once, now
# that every cohort's metadata (including Arthur and Hilton, previously
# added much later -- after the main aSHM pull had already run without them)
# is available from the start.
# ============================================================================
sample_data <- list()

sample_data$meta <- bind_rows(
    bl_data$meta_to_bundle,
    fl_data$meta_to_bundle,
    dlbcl_data$meta_to_bundle,
    cell_lines_data$meta_to_bundle
)

sample_data$meta <- sample_data$meta %>%
    select(-COO_consensus, -lymphgen, -EBV_status_inf) %>%
    left_join(
        .,
        get_gambl_metadata() %>%
            select(
                sample_id,
                COO_consensus,
                lymphgen,
                EBV_status_inf
            )
    ) %>%
    select(all_of(colnames_for_bundled_meta))

sample_data$meta <- bind_rows(
    sample_data$meta,
    all_capture_meta,
    arthur_meta,
    trios_meta
)

### begin metadata fixing
# This preserves the original cohort column and ensures there are no duplicates
# in the metadata
fix <- sample_data$meta
fix <- fix %>% rename(study = cohort)

fix <- left_join(
    fix,
    get_gambl_metadata() %>%
        select(sample_id, seq_type, cohort)
)

fix <- fix %>% filter(!is.na(study))

fix <- distinct(fix)

sample_data$meta <- fix
### end metadata fixing

# sample_study: a many-to-many bridge table (one row per (sample_id, study)
# pair, not one row per sample), replacing the old inline `Study =` tagging
# on mutation rows. A sample belonging to N studies is simply N rows here --
# no schema change needed as multi-study overlap becomes more common.
sample_data$sample_study <- bind_rows(
    thomas_bl_study,
    thomas_dlbcl_study,
    dreval_study,
    reddy_study,
    schmitz_study,
    chapuy_study,
    golub_study,
    arthur_study,
    hilton_study
) %>% distinct()


#add SVs
# Adding the manta SVs for published studies. Scoped to sample_data$meta (the
# metadata just finalized above, for this run) rather than
# GAMBLR.data::sample_data$meta (whatever happens to be installed) so the SV
# sample scope is self-consistent with the rest of this bundle instead of
# drifting with installed-package state across runs.
full_genome_meta <- get_gambl_metadata(seq_type_filter = "genome")

bundled_meta <- full_genome_meta %>%
    filter(
        sample_id %in% sample_data$meta$sample_id
    )

full_sv_to_bundle <- get_manta_sv(
        these_samples_metadata = bundled_meta,
        projection = "hg38"
    )

annotated_sv_to_bundle <- annotate_sv(
    full_sv_to_bundle,
    genome_build = "hg38"
)
annotated_sv_to_bundle <- annotated_sv_to_bundle %>%
    filter(!is.na(partner)) %>%
    mutate(
        chrom1 = paste0("chr", chrom1),
        chrom2 = paste0("chr", chrom2)
    )

# drop all annotation columns to restore original data subset just to the putative driver SVs
annotated_sv_keep <- left_join(
    full_sv_to_bundle,
    annotated_sv_to_bundle,
    by = c(
        "CHROM_A" = "chrom1",
        "CHROM_B" = "chrom2",
        "START_A" = "start1",
        "tumour_sample_id")
    ) %>%
    dplyr::filter(!is.na(partner)) %>%
    select(c(1:16))

# Now same for the grch37 projection
full_sv_to_bundle_grch37 <- get_manta_sv(
    these_samples_metadata = bundled_meta
)

annotated_sv_to_bundle_grch37 <- annotate_sv(
    full_sv_to_bundle_grch37
)

annotated_sv_to_bundle_grch37 <- annotated_sv_to_bundle_grch37 %>%
    filter(!is.na(partner))

# drop all annotation columns to restore original data subset just to the putative driver SVs
annotated_sv_keep_grch37 <- left_join(
    full_sv_to_bundle_grch37,
    annotated_sv_to_bundle_grch37,
    by = c(
        "CHROM_A" = "chrom1",
        "CHROM_B" = "chrom2",
        "START_A" = "start1",
        "tumour_sample_id")
    ) %>%
    filter(!is.na(partner)) %>%
    select(c(1:16))

sample_data$grch37$bedpe <- annotated_sv_keep_grch37
sample_data$hg38$bedpe <- annotated_sv_keep

sample_data$grch37$seg <- bind_rows(
    fl_data$cnv_to_bundle,
    cell_lines_data$grch37$cnv_to_bundle
)

sample_data$hg38$seg <- bind_rows(
    bl_data$cnv_to_bundle,
    dlbcl_data$cnv_to_bundle,
    cell_lines_data$hg38$cnv_to_bundle
)


# ============================================================================
# Phase 3: Publication-pipeline SSM data. Genuinely distinct per-paper data
# (not a redundant SLMS-3 recall), kept as separate, un-deduplicated pulls --
# tagged Pipeline="Publication" only; cohort membership lives in
# sample_data$sample_study, not an inline Study column.
# ============================================================================

# Publication-pipeline data is read directly from each paper's own
# supplementary file, which isn't guaranteed to use GAMBL's own
# Tumor_Sample_Barcode convention -- Arthur's now-removed raw dump was the
# clearest example, using bare patient-style IDs instead of GAMBL's real
# per-sample naming. Defensive relabel: for any row whose
# Tumor_Sample_Barcode matches a study's own identifier
# (sample_study$study_id) rather than GAMBL's real sample_id, replace it
# with the real sample_id. A no-op wherever Tumor_Sample_Barcode already IS
# the real sample_id -- true for Thomas/Dreval/Hilton, where study_id
# mirrors sample_id (see thomas_bl_study above).
relabel_to_sample_id <- function(df, study_name) {
    lookup <- sample_data$sample_study %>%
        filter(study == study_name, !is.na(study_id)) %>%
        select(study_id, .gambl_sample_id = sample_id)
    df %>%
        left_join(lookup, by = c("Tumor_Sample_Barcode" = "study_id")) %>%
        mutate(Tumor_Sample_Barcode = coalesce(.gambl_sample_id, Tumor_Sample_Barcode)) %>%
        select(-.gambl_sample_id)
}

# This is needed for the proteinpainter compatibility. Reads from the
# *installed* GAMBLR.data::sample_data (not the locally-built sample_data
# above) -- a pre-existing, out-of-scope inconsistency, left as-is.
GAMBLR.data::sample_data$meta$cohort %>% table

selected_columns <- c(
        "Tumor_Sample_Barcode", "Hugo_Symbol",
        "NCBI_Build", "Chromosome", "Start_Position", "End_Position",
        "Tumor_Seq_Allele2", maf_columns_to_keep
)

these_samples <- GAMBLR.data::sample_data$meta %>%
    filter(cohort %in% c("BL_Thomas")) %>%
    pull(sample_id)

these_samples_dlbcl <- GAMBLR.data::sample_data$meta %>%
    filter(cohort %in% c("DLBCL_Thomas", "DLBCL_cell_lines")) %>%
    pull(sample_id)

coding_maf <- read_tsv("/projects/adult_blgsp/results_manuscript/BL.hg38.CDS.maf") %>% # get from flat maf file to show SSM in hg38 coordinates similar to the original manuscript
    filter(Tumor_Sample_Barcode %in% these_samples & # drop BL58 cell line
           ! str_detect(Tumor_Sample_Barcode, "^SP|^06")) %>% # drop ICGC and 1 LLMPP case
    select(
        all_of(selected_columns)
    )

coding_maf_dlbcl <- get_coding_ssm(
        these_samples_metadata =  get_gambl_metadata() %>%
            filter(sample_id %in% these_samples_dlbcl),
        this_seq_type = "genome",
        projection = "hg38",
        basic_columns = FALSE
    ) %>%
    select(
        all_of(selected_columns)
    )

coding_maf <- bind_rows(
    coding_maf,
    coding_maf_dlbcl
)

# hg38 Publication rows (Thomas BL + Thomas DLBCL), enriched via the
# proteinpainter-compatibility join above. This used to be applied directly
# to sample_data$hg38$maf once it existed early in the script; now it's
# applied to the Publication-only rows before Phase 6 combines them with the
# consolidated SLMS-3 pull.
hg38_publication_rows <- bind_rows(
    bl_data$ssm_to_bundle %>% mutate(Pipeline = "Publication"),
    dlbcl_data$ssm_to_bundle %>% mutate(Pipeline = "Publication")
) %>%
    relabel_to_sample_id("Thomas") %>%
    left_join(coding_maf)

this_study_samples <- GAMBLR.data::sample_data$meta %>%
    filter(cohort %in% c("FL_Dreval", "DLBCL_cell_lines")) %>%
    pull(sample_id)

# FLs in grch37
coding_maf <- time_it("coding_maf get_ssm_by_samples", {
    get_ssm_by_samples(
        these_samples_metadata =  get_gambl_metadata() %>%
                filter(sample_id %in% this_study_samples),
        basic_columns = FALSE) %>%
        select(
            all_of(selected_columns)
        )
})

fl_data$ssm_to_bundle <- fl_data$ssm_to_bundle %>%
    relabel_to_sample_id("Dreval") %>%
    dplyr::left_join(
        coding_maf
    ) %>%
    distinct()

# grch37 Publication rows (Dreval FL + Reddy's original variants).
grch37_publication_rows <- bind_rows(
    fl_data$ssm_to_bundle %>% mutate(Pipeline = "Publication"),
    reddy_original_maf %>% relabel_to_sample_id("Reddy") %>% mutate(Pipeline = "Publication")
)

setwd(PKG_ROOT)


# ============================================================================
# Phase 4: one consolidated SLMS-3 pull for every sample in the bundle,
# instead of the previous per-cohort blocks (all-capture, cell lines, Hilton
# trios, Arthur) each independently deciding which samples to pull for --
# that pattern silently duplicated mutation rows whenever a sample belonged
# to more than one cohort block. Pulling once per (sample, genome_build),
# keyed off the now-complete deduplicated sample_data$meta, makes
# duplication structurally impossible regardless of how many cohorts a
# sample belongs to.
# ============================================================================
# Cell lines are excluded here -- they get their own separate, genome-wide
# pull below (not restricted to the lymphoma gene panel), so they must not
# also go through this panel-restricted pull or they'd be pulled twice.
#
# sample_data$meta is used ONLY to decide which sample_ids belong in this
# pull -- never passed directly to a GAMBLR.results call. Cohorts built by
# hand (BL_Thomas, DLBCL_Thomas, FL_Dreval, DLBCL_cell_lines) don't have
# every column a fresh get_gambl_metadata() pull would (e.g. unix_group is
# NA for all of them), which silently caused get_ssm_by_regions() to miss
# real coding-classified calls for those samples (confirmed: a raw pull for
# one such sample had real Missense_Mutation/Nonsense_Mutation/Silent rows
# that never made it into the assembled maf table). Re-fetching complete,
# live metadata for exactly this sample_id set avoids that entirely.
all_slms3_meta <- sample_data$meta %>%
    filter(seq_type %in% c("genome", "capture"),
           ! sample_id %in% cell_lines_data$meta$sample_id)
all_slms3_meta <- get_gambl_metadata() %>%
    filter(sample_id %in% all_slms3_meta$sample_id)

slms3_grch37 <- bind_rows(
    time_it("SLMS-3 genome grch37", pull_data(all_slms3_meta %>% filter(seq_type == "genome"))),
    time_it("SLMS-3 capture grch37", pull_data(all_slms3_meta %>% filter(seq_type == "capture")))
) %>% mutate(Pipeline = "SLMS-3")
print("Done collecting grch37 SLMS-3")

slms3_hg38 <- bind_rows(
    time_it("SLMS-3 genome hg38", pull_data(all_slms3_meta %>% filter(seq_type == "genome"), "hg38")),
    time_it("SLMS-3 capture hg38", pull_data(all_slms3_meta %>% filter(seq_type == "capture"), "hg38"))
) %>% mutate(Pipeline = "SLMS-3")
print("Done collecting hg38 SLMS-3")

# Cell lines get their own, separate, genome-WIDE SNV pull -- not restricted
# to the lymphoma gene panel like every other cohort above. This matches
# their original pre-refactor behaviour, which was lost when they were first
# folded into the panel-restricted consolidated pull above (confirmed via
# compare_bundle_changes.R: cell lines showed ~50-65k "lost" rows per sample
# against the old bundle, correctly diagnosed as a real scope reduction, not
# a bug worth working around).
cell_lines_ssm_grch37 <- time_it("cell lines get_ssm_by_samples grch37", {
    get_ssm_by_samples(
        these_samples_metadata = cell_lines_data$meta,
        basic_columns = FALSE
    ) %>% select(all_of(all_cols)) %>% mutate(Pipeline = "SLMS-3")
})

cell_lines_ssm_hg38 <- time_it("cell lines get_ssm_by_samples hg38", {
    get_ssm_by_samples(
        these_samples_metadata = cell_lines_data$meta,
        projection = "hg38",
        basic_columns = FALSE
    ) %>% select(all_of(all_cols)) %>% mutate(Pipeline = "SLMS-3")
})

slms3_grch37 <- bind_rows(slms3_grch37, cell_lines_ssm_grch37)
slms3_hg38 <- bind_rows(slms3_hg38, cell_lines_ssm_hg38)
print("Done collecting cell line SLMS-3 (genome-wide)")


# ============================================================================
# Phase 5: one consolidated aSHM pull, using the now-complete sample_data$meta
# (includes Arthur + Hilton from the start). This replaces both the old main
# aSHM block (which ran before Arthur/Hilton were added to sample_data$meta,
# so neither got any aSHM coverage from it) and Hilton's separate dedicated
# aSHM block (which existed only to compensate for that gap). Arthur gets
# aSHM coverage for the first time as a result.
# ============================================================================
# sample_data$meta is used only to fix the sample_id set -- see the same
# rationale next to all_slms3_meta above for why a fresh get_gambl_metadata()
# pull is used for the actual GAMBLR.results call instead of sample_data$meta
# directly.
ashm_pull_meta <- get_gambl_metadata() %>%
    filter(sample_id %in% sample_data$meta$sample_id)

regions_bed_grch37 <- create_bed_data(
    grch37_ashm_regions,
    fix_names = "concat",
    concat_cols = c("gene", "region"), sep = "-"
)

grch37_ashm <- time_it("grch37_ashm get_ssm_by_regions", {
    get_ssm_by_regions(
        these_samples_metadata = ashm_pull_meta,
        regions_bed = regions_bed_grch37,
        streamlined = FALSE,
        basic_columns = FALSE
    ) %>%
        select(any_of(all_cols))
})

grch37_ashm <- grch37_ashm %>%
    filter(Tumor_Sample_Barcode %in% sample_data$meta$Tumor_Sample_Barcode) %>%
    mutate(Pipeline = "SLMS-3")

regions_bed_hg38 <- create_bed_data(
    hg38_ashm_regions,
    fix_names = "concat",
    concat_cols = c("gene", "region"), sep = "-"
)

hg38_ashm <- time_it("hg38_ashm get_ssm_by_regions", {
    get_ssm_by_regions(
        these_samples_metadata = ashm_pull_meta,
        regions_bed = regions_bed_hg38,
        projection = "hg38",
        streamlined = FALSE,
        basic_columns = FALSE
    ) %>%
        select(any_of(all_cols))
})
hg38_ashm <- hg38_ashm %>%
    filter(Tumor_Sample_Barcode %in% sample_data$meta$Tumor_Sample_Barcode) %>%
    mutate(Pipeline = "SLMS-3")

sample_data$grch37$ashm <- grch37_ashm %>% distinct()
sample_data$hg38$ashm <- hg38_ashm %>% distinct()

print("done extracting aSHM mutations from GAMBLR.results")


# ============================================================================
# Phase 6: final maf assembly -- Publication rows + the consolidated SLMS-3
# pull, per genome build. No `Study =` in any mutate() anywhere in this
# script; cohort membership lives entirely in sample_data$sample_study.
# ============================================================================
sample_data$grch37$maf <- bind_rows(
    grch37_publication_rows,
    slms3_grch37
)

sample_data$hg38$maf <- bind_rows(
    hg38_publication_rows,
    slms3_hg38
)

print("done extracting all mutations in lymphoma genes with GAMBLR.results")


# --- Persist the assembled data -------------------------------------------
# Previously this bundled the multi-GB `sample_data.rda`. Instead we now write:
#   1. a lightweight, bundled `sample_metadata` object (metadata only), and
#   2. the large per-sample frames to gambl_mutations.db (a release asset,
#      NOT shipped in the package tarball).
# This is built straight from the in-memory `sample_data` above, so there is no
# sample_data.rda round-trip.

# 1. lightweight bundled metadata (replaces sample_data$meta lookups)
sample_metadata <- sample_data$meta
usethis::use_data(sample_metadata, overwrite = TRUE, compress = "xz")

print("Starting sqlite build")
# 2. large sample-level frames -> SQLite
source("data-raw/write_mutations_db.R")
write_mutations_db(
    sample_data,
    out_db = "gambl_mutations.db",
    source_desc = paste0("assemble_bundled_data.R @ ", format(Sys.Date()))
)

# During the transition you may still want the legacy monolithic object
# (e.g. until other GAMBLR packages are audited off sample_data). Flip to TRUE
# to also write data/sample_data.rda.
WRITE_LEGACY_SAMPLE_DATA <- FALSE
if (WRITE_LEGACY_SAMPLE_DATA) {
    usethis::use_data(sample_data, overwrite = TRUE, compress = "xz")
}

library(data.tree)

tree <- FromListSimple(sample_data)
tree
