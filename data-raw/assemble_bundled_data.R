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

# restrict to the most inclusive DLBCL gene list
all_lymphoma_genes <- lymphoma_genes_comprehensive$Gene

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


# Importing BL data from Thomas et al
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


# Importing FL data from Dreval et al
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



# Importing DLBCL data in hg38 from Thomas et al
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

##### Importing SLMS-3 variants
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


# Importing DLBCL capture data
reddy_data <- list()
schmitz_data <- list()
chapuy_data <- list()
golub_data <- list()

# Importing metadata from Reddy et al and updating IDs to be consistent with GAMBL metadata
reddy_meta <- read_excel(
        "inst/extdata/studies/DLBCL_Reddy.xlsx",
        sheet = 1
    ) %>%
    mutate(
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
    ) %>%
    dplyr::select(
        sample_id,
        patient_id,
        Tumor_Sample_Barcode,
        sex = Gender,
        COO_consensus
    )

# (config is resolved via config.yml / the GAMBLR.results fallback; no chdir to
# a gambl repo is needed to locate data — repo_base is absolute)

reddy_meta_gambl <- get_gambl_metadata() %>%
    dplyr::filter(cohort == "dlbcl_reddy") %>%
    dplyr::select(
        sample_id, lymphgen, EBV_status_inf, cohort, pathology, seq_type,
        unix_group, genome_build, pairing_status, normal_sample_id
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

all_capture_grch37_ssm_to_bundle <- time_it("pull_data(all_capture_meta, grch37)", pull_data(all_capture_meta))
print("Done collecting grch37")
all_capture_hg38_ssm_to_bundle <- time_it("pull_data(all_capture_meta, hg38)", pull_data(all_capture_meta, "hg38"))
print("Done collecting hg38")


# Importing DLBCL cell lines
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

cell_lines_data$grch37$ssm_to_bundle <- time_it("cell_lines get_ssm_by_samples grch37", {
    get_ssm_by_samples(
        these_samples_metadata = cell_lines_data$meta,
        basic_columns = FALSE
    ) %>% select(all_of(all_cols))
})

cell_lines_data$hg38$ssm_to_bundle <- time_it("cell_lines get_ssm_by_samples hg38", {
    get_ssm_by_samples(
        these_samples_metadata = cell_lines_data$meta,
        projection = "hg38",
        basic_columns = FALSE
    ) %>% select(all_of(all_cols))
})

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


# Manta SVs for published studies: moved to after sample_data$meta is
# finalized below (see "Adding the manta SVs for published studies"), so the
# sample scope is the freshly-assembled local metadata for this run rather
# than whatever GAMBLR.data happens to be installed.

# Combine everything together
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
    all_capture_meta
)

sample_data$hg38$maf <- bind_rows(
    bl_data$ssm_to_bundle %>% mutate(
        Pipeline = "Publication",
        Study = "Thomas"
    ),
    dlbcl_data$ssm_to_bundle %>% mutate(
        Pipeline = "Publication",
        Study = "Thomas"
    ),
    cell_lines_data$hg38$ssm %>% mutate(
        Pipeline = "SLMS-3",
        Study = NA
    ),
    all_capture_hg38_ssm_to_bundle %>% mutate(
        Pipeline = "SLMS-3",
        Study = case_when(
            Tumor_Sample_Barcode %in% reddy_data$meta$sample_id ~ "Reddy",
            Tumor_Sample_Barcode %in% schmitz_data$meta$sample_id ~ "Schmitz",
            Tumor_Sample_Barcode %in% chapuy_data$meta$sample_id ~ "Chapuy",
            Tumor_Sample_Barcode %in% golub_data$meta$sample_id ~ "NCI_Golub"
        )
    )
)

sample_data$grch37$seg <- bind_rows(
    fl_data$cnv_to_bundle,
    cell_lines_data$grch37$cnv_to_bundle
)

sample_data$hg38$seg <- bind_rows(
    bl_data$cnv_to_bundle,
    dlbcl_data$cnv_to_bundle,
    cell_lines_data$hg38$cnv_to_bundle
)

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

# This is needed for the proteinpainter compatibility
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

dim(GAMBLR.data::sample_data$hg38$maf)

sample_data$hg38$maf <- sample_data$hg38$maf %>%
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
    dplyr::left_join(
        coding_maf
    ) %>%
    distinct()
    
sample_data$grch37$maf <- bind_rows(
    fl_data$ssm_to_bundle %>% mutate(
        Pipeline = "Publication",
        Study = "Dreval"
    ),
    cell_lines_data$grch37$ssm %>% mutate(
        Pipeline = "SLMS-3",
        Study = NA
    ),
    all_capture_grch37_ssm_to_bundle %>% mutate(
        Pipeline = "SLMS-3",
        Study = case_when(
            Tumor_Sample_Barcode %in% reddy_data$meta$sample_id ~ "Reddy",
            Tumor_Sample_Barcode %in% schmitz_data$meta$sample_id ~ "Schmitz",
            Tumor_Sample_Barcode %in% chapuy_data$meta$sample_id ~ "Chapuy",
            Tumor_Sample_Barcode %in% golub_data$meta$sample_id ~ "NCI_Golub"
        )
    )
)

regions_bed_grch37 = GAMBLR.utils::create_bed_data(
    GAMBLR.data::grch37_ashm_regions ,
    fix_names = "concat",
    concat_cols = c("gene","region"),sep="-"
)

# Add aSHM mutations for the already released samples
grch37_ashm <- time_it("grch37_ashm get_ssm_by_regions", {
    get_ssm_by_regions(
        these_samples_metadata = sample_data$meta,
        regions_bed = regions_bed_grch37,
        streamlined = FALSE,
        basic_columns = FALSE
    ) %>%
        select(
            any_of(c(colnames(sample_data$grch37$maf), maf_columns_to_keep))
        )
})

grch37_ashm <- grch37_ashm %>%
    filter(Tumor_Sample_Barcode %in% sample_data$meta$Tumor_Sample_Barcode)

grch37_ashm <- grch37_ashm %>% mutate(Pipeline = "SLMS-3")

studies <- bind_rows(
    sample_data$grch37$maf %>%
        distinct(Tumor_Sample_Barcode, Study),
    sample_data$hg38$maf %>%
        distinct(Tumor_Sample_Barcode, Study)
) %>%
distinct()

grch37_ashm <- left_join(
    grch37_ashm,
    studies
)

hg38_ashm <- time_it("hg38_ashm get_ssm_by_regions", {
    get_ssm_by_regions(
        these_samples_metadata = sample_data$meta,
        regions_bed = GAMBLR.utils::create_bed_data(
            GAMBLR.data::hg38_ashm_regions,
            fix_names = "concat",
            concat_cols = c("gene","region"),sep="-"
        ),
        projection = "hg38",
        streamlined = FALSE,
        basic_columns = FALSE
    ) %>%
        select(
            any_of(c(colnames(sample_data$hg38$maf), maf_columns_to_keep))
        )
})
hg38_ashm <- hg38_ashm %>%
    filter(Tumor_Sample_Barcode %in% sample_data$meta$Tumor_Sample_Barcode)

hg38_ashm <- hg38_ashm %>% mutate(Pipeline = "SLMS-3")

hg38_ashm <- left_join(
    hg38_ashm,
    studies
)

sample_data$grch37$ashm <- grch37_ashm
sample_data$hg38$ashm <- hg38_ashm

print("done extracting aSHM mutations from GAMBLR.results")

# Now add the SLMS-3 calls in both projections for those samples that
# are bundled as publication data
publication_samples_grch37 <- sample_data$grch37$maf %>%
    filter(Pipeline == "Publication") %>%
    pull(Tumor_Sample_Barcode) %>%
    unique %>% sort

publication_samples_hg38 <- sample_data$hg38$maf %>%
    filter(Pipeline == "Publication") %>%
    pull(Tumor_Sample_Barcode) %>%
    unique %>% sort

publication_samples <- c(
    publication_samples_grch37,
    publication_samples_hg38
)
print("extracting grch37 mutations in lymphoma genes with GAMBLR.results")

sample_data$grch37$maf <- time_it("publication-samples get_ssm_by_regions grch37", {
    get_ssm_by_regions(
        regions_bed = lymphoma_genes_bed_grch37,
        these_samples_metadata = get_gambl_metadata() %>%
            filter(sample_id %in% publication_samples),
        basic_columns = FALSE) %>%
        filter(Hugo_Symbol %in% all_lymphoma_genes) %>%
        mutate(Pipeline = "SLMS-3") %>%
        left_join(
            .,
            studies
        ) %>%
        select(colnames(sample_data$grch37$maf)) %>%
        bind_rows(
            .,
            sample_data$grch37$maf
        )
})
print("extracting hg38 mutations in lymphoma genes with GAMBLR.results")

sample_data$hg38$maf <- time_it("publication-samples get_ssm_by_regions hg38", {
    get_ssm_by_regions(
        regions_bed = lymphoma_genes_bed_hg38,
        these_samples_metadata = get_gambl_metadata() %>%
            filter(sample_id %in% publication_samples),
        projection = "hg38",
        basic_columns = FALSE) %>%
        filter(Hugo_Symbol %in% all_lymphoma_genes) %>%
        mutate(Pipeline = "SLMS-3") %>%
        left_join(
            .,
            studies
        ) %>%
        select(colnames(sample_data$hg38$maf)) %>%
        bind_rows(
            .,
            sample_data$hg38$maf
        )
})
print("done extracting all mutations in lymphoma genes with GAMBLR.results")

setwd(PKG_ROOT)

# Add data from Arthur paper
arthur_maf <- read_tsv(
    "inst/extdata/studies/DLBCL_Arthur.maf.gz"
)
sample_data$grch37$maf <- bind_rows(
    sample_data$grch37$maf,
    arthur_maf %>%
        mutate(
            Pipeline = "strelka",
            Study = "Arthur"
        )
)

arthur_meta <- read_xlsx(
    "inst/extdata/studies/DLBCL_Arthur.xlsx",
    sheet = 1
) %>% filter(`WGS data` == 1)


arthur_meta <- gambl_metadata %>%
    filter(
        patient_id %in% arthur_meta$`Case ID`,
        seq_type == "genome",
        ! grepl("tumor", sample_id)
    ) %>%
    mutate(
        sample_id = patient_id,
        Tumor_Sample_Barcode = patient_id,
        cohort = "DLBCL_Arthur",
        reference_PMID = pmids$Arthur_DLBCL
    )

sample_data$meta <- bind_rows(
    sample_data$meta,
    arthur_meta
)

# Add data from Hilton trios paper
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
    select(any_of(colnames(sample_data$meta))) %>%
    mutate(
        cohort = "DLBCL_Hilton",
        reference_PMID = pmids$Hilton_DLBCL
    )

sample_data$meta <- bind_rows(
    sample_data$meta,
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

# trios grch37 ssm
genome_trios_ssm_grch37 <- time_it("trios get_ssm_by_regions genome grch37", {
    get_ssm_by_regions(
        regions_bed = lymphoma_genes_bed_grch37,
        these_samples_metadata = trios_meta %>%
            filter(seq_type == "genome"),
        basic_columns = FALSE
    ) %>%
        filter(Hugo_Symbol %in% all_lymphoma_genes) %>%
        mutate(
            Pipeline = "SLMS-3",
            Study = "Hilton"
        ) %>%
        select(all_of(colnames(sample_data$grch37$maf)))
})
print("extracting mutations for Trios cohort")
capture_trios_ssm_grch37 <- time_it("trios get_ssm_by_regions capture grch37", {
    get_ssm_by_regions(
        regions_bed = lymphoma_genes_bed_grch37,
        these_samples_metadata = trios_meta %>%
            filter(seq_type == "capture"),
        basic_columns = FALSE
    ) %>%
        filter(Hugo_Symbol %in% all_lymphoma_genes) %>%
        mutate(
            Pipeline = "SLMS-3",
            Study = "Hilton"
        ) %>%
        select(all_of(colnames(sample_data$grch37$maf)))
})

trios_ssm_grch37 <- bind_rows(
    genome_trios_ssm_grch37,
    capture_trios_ssm_grch37
)

# trios hg38 ssm
genome_trios_ssm_hg38 <- time_it("trios get_ssm_by_regions genome hg38", {
    get_ssm_by_regions(
        regions_bed = lymphoma_genes_bed_hg38,
        these_samples_metadata = trios_meta %>%
            filter(seq_type == "genome"),
        basic_columns = FALSE,
        projection = "hg38"
    ) %>%
        filter(Hugo_Symbol %in% all_lymphoma_genes) %>%
        mutate(
            Pipeline = "SLMS-3",
            Study = "Hilton"
        ) %>%
        select(all_of(colnames(sample_data$hg38$maf)))
})

capture_trios_ssm_hg38 <- time_it("trios get_ssm_by_regions capture hg38", {
    get_ssm_by_regions(
        regions_bed = lymphoma_genes_bed_hg38,
        these_samples_metadata = trios_meta %>%
            filter(seq_type == "capture"),
        basic_columns = FALSE,
        projection = "hg38"
    ) %>%
        filter(Hugo_Symbol %in% all_lymphoma_genes) %>%
        mutate(
            Pipeline = "SLMS-3",
            Study = "Hilton"
        ) %>%
        select(all_of(colnames(sample_data$hg38$maf)))
})

trios_ssm_hg38 <- bind_rows(
    genome_trios_ssm_hg38,
    capture_trios_ssm_hg38
)

sample_data$grch37$maf <- bind_rows(
    sample_data$grch37$maf,
    trios_ssm_grch37
)

sample_data$hg38$maf <- bind_rows(
    sample_data$hg38$maf,
    trios_ssm_hg38
)

regions_bed <- create_bed_data(
    grch37_ashm_regions,
    fix_names = "concat",
    concat_cols = c("gene", "region"),
    sep = "-"
)

trios_ashm_grch37 <- time_it("trios_ashm get_ssm_by_regions grch37", {
    get_ssm_by_regions(
        these_samples_metadata = trios_meta,
        regions_bed = regions_bed,
        streamlined = FALSE,
        basic_columns = FALSE
    )
})

trios_ashm_grch37 <- trios_ashm_grch37 %>%
    filter(Tumor_Sample_Barcode %in% trios_meta$Tumor_Sample_Barcode)

trios_ashm_grch37 <- trios_ashm_grch37 %>%
    mutate(
        Pipeline = "SLMS-3",
        Study = "Hilton"
    ) %>%
    select(all_of(colnames(sample_data$grch37$maf)))


sample_data$grch37$ashm <- bind_rows(
    sample_data$grch37$ashm,
    trios_ashm_grch37
) %>% distinct



regions_bed <- create_bed_data(
    hg38_ashm_regions,
    fix_names = "concat",
    concat_cols = c("gene", "region"),
    sep = "-"
)

trios_ashm_hg38 <- time_it("trios_ashm get_ssm_by_regions hg38", {
    get_ssm_by_regions(
        these_samples_metadata = trios_meta,
        regions_bed = regions_bed,
        projection = "hg38",
        streamlined = FALSE,
        basic_columns = FALSE
    )
})
trios_ashm_hg38 <- trios_ashm_hg38 %>%
    filter(Tumor_Sample_Barcode %in% trios_meta$Tumor_Sample_Barcode)

trios_ashm_hg38 <- trios_ashm_hg38 %>%
    mutate(
        Pipeline = "SLMS-3",
        Study = "Hilton"
    ) %>%
    select(all_of(colnames(sample_data$hg38$maf)))

sample_data$hg38$ashm <- bind_rows(
    sample_data$hg38$ashm,
    trios_ashm_hg38
) %>% distinct

print("Done getting aSHM from Trios")

setwd(PKG_ROOT)

# Add data from Reddy paper
reddy_original_maf <- read_tsv(
    "inst/extdata/studies/reddy_original_variants_with_VAF.maf.gz"
) %>%
select(any_of(colnames(sample_data$grch37$maf)))

sample_data$grch37$maf <- bind_rows(
    sample_data$grch37$maf,
    reddy_original_maf %>%
        mutate(
            Pipeline = "Publication",
            Study = "Reddy"
        )
)

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
