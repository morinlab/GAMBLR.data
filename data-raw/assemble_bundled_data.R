setwd("~/my_dir/repos/GAMBLR.data/")

library(readxl)
library(GAMBLR)
library(parallel)
library(tidyverse)

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
    slms3 <- get_ssm_by_samples(
        these_samples_metadata = pull_meta,
        basic_columns = FALSE,
        projection = pull_projection
    ) %>%
    select(
        all_of(all_cols)
    ) %>%
    filter(
        Hugo_Symbol %in% all_lymphoma_genes
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

setwd("/projects/rmorin/projects/gambl-repos/gambl-kdreval/")

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

# warning: this is very slow!
all_capture_grch37_ssm_to_bundle <- pull_data(all_capture_meta)
all_capture_hg38_ssm_to_bundle <- pull_data(all_capture_meta, "hg38")

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

cell_lines_data$grch37$ssm_to_bundle <- get_ssm_by_samples(
    these_samples_metadata = cell_lines_data$meta,
    basic_columns = FALSE
) %>% select(all_of(all_cols))

cell_lines_data$hg38$ssm_to_bundle <- get_ssm_by_samples(
    these_samples_metadata = cell_lines_data$meta,
    projection = "hg38",
    basic_columns = FALSE
) %>% select(all_of(all_cols))

cell_lines_data$grch37$cnv_to_bundle <- get_sample_cn_segments(
    these_sample_ids = cell_lines_data$meta$sample_id
)

cell_lines_data$hg38$cnv_to_bundle <- get_sample_cn_segments(
    these_sample_ids = cell_lines_data$meta$sample_id,
    projection = "hg38",
    with_chr_prefix = TRUE
)

cell_lines_data$grch37$sv_to_bundle <- get_manta_sv(
    these_samples_metadata = cell_lines_data$meta,
)

cell_lines_data$hg38$sv_to_bundle <- get_manta_sv(
    these_samples_metadata = cell_lines_data$meta,
    projection = "hg38"
)

# Adding the manta SVs for published studies
full_genome_meta <- get_gambl_metadata(seq_type_filter = "genome")

bundled_meta <- full_genome_meta %>%
    filter(
        sample_id %in% GAMBLR.data::sample_data$meta$sample_id
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
    cell_lines_data$grch37$cnv
)

sample_data$hg38$seg <- bind_rows(
    bl_data$cnv_to_bundle,
    dlbcl_data$cnv_to_bundle,
    cell_lines_data$hg38$cnv
)

#add SVs
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
coding_maf <- get_ssm_by_samples(
    these_samples_metadata =  get_gambl_metadata() %>%
            filter(sample_id %in% this_study_samples),
    basic_columns = FALSE) %>%
    select(
        all_of(selected_columns)
    )

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


# Add aSHM mutations for the already released samples
grch37_ashm <- get_ssm_by_regions(
    these_samples_metadata = sample_data$meta,
    regions_bed = GAMBLR.utils::create_bed_data(
        GAMBLR.data::grch37_ashm_regions,
        fix_names = "concat",
        concat_cols = c("gene","region"),sep="-"
    ),
    streamlined = FALSE,
    basic_columns = FALSE
) %>%
    select(
        any_of(c(colnames(sample_data$grch37$maf), maf_columns_to_keep))
    )

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

hg38_ashm <- get_ssm_by_regions(
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
hg38_ashm <- hg38_ashm %>%
    filter(Tumor_Sample_Barcode %in% sample_data$meta$Tumor_Sample_Barcode)

hg38_ashm <- hg38_ashm %>% mutate(Pipeline = "SLMS-3")

hg38_ashm <- left_join(
    hg38_ashm,
    studies
)

sample_data$grch37$ashm <- grch37_ashm
sample_data$hg38$ashm <- hg38_ashm


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

sample_data$grch37$maf <- get_ssm_by_samples(
    these_samples_metadata = get_gambl_metadata() %>%
        filter(sample_id %in% publication_samples),
    basic_columns = FALSE) %>%
    filter(
        Hugo_Symbol %in% all_lymphoma_genes
    ) %>%
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

sample_data$hg38$maf <- get_ssm_by_samples(
    these_samples_metadata = get_gambl_metadata() %>%
        filter(sample_id %in% publication_samples),
    projection = "hg38",
    basic_columns = FALSE) %>%
    filter(
        Hugo_Symbol %in% all_lymphoma_genes
    ) %>%
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


setwd("~/my_dir/repos/GAMBLR.data/")

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

setwd("/projects/rmorin/projects/gambl-repos/gambl-kdreval/")

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
genome_trios_ssm_grch37 <- get_ssm_by_samples(
    these_samples_metadata = trios_meta %>%
        filter(seq_type == "genome"),
    basic_columns = FALSE,
    subset_from_merge = TRUE
) %>%
    filter(Hugo_Symbol %in% all_lymphoma_genes) %>%
    mutate(
        Pipeline = "SLMS-3",
        Study = "Hilton"
    ) %>%
    select(all_of(colnames(sample_data$grch37$maf)))

capture_trios_ssm_grch37 <- get_ssm_by_samples(
    these_samples_metadata = trios_meta %>%
        filter(seq_type == "capture"),
    basic_columns = FALSE,
    subset_from_merge = TRUE
) %>%
    filter(Hugo_Symbol %in% all_lymphoma_genes) %>%
    mutate(
        Pipeline = "SLMS-3",
        Study = "Hilton"
    ) %>%
    select(all_of(colnames(sample_data$grch37$maf)))

trios_ssm_grch37 <- bind_rows(
    genome_trios_ssm_grch37,
    capture_trios_ssm_grch37
)

# trios hg38 ssm
genome_trios_ssm_hg38 <- get_ssm_by_samples(
    these_samples_metadata = trios_meta %>%
        filter(seq_type == "genome"),
    basic_columns = FALSE,
    subset_from_merge = TRUE,
    projection = "hg38"
) %>%
    filter(Hugo_Symbol %in% all_lymphoma_genes) %>%
    mutate(
        Pipeline = "SLMS-3",
        Study = "Hilton"
    ) %>%
    select(all_of(colnames(sample_data$hg38$maf)))

capture_trios_ssm_hg38 <- get_ssm_by_samples(
    these_samples_metadata = trios_meta %>%
        filter(seq_type == "capture"),
    basic_columns = FALSE,
    subset_from_merge = TRUE,
    projection = "hg38"
) %>%
    filter(Hugo_Symbol %in% all_lymphoma_genes) %>%
    mutate(
        Pipeline = "SLMS-3",
        Study = "Hilton"
    ) %>%
    select(all_of(colnames(sample_data$hg38$maf)))

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

trios_ashm_grch37 <- get_ssm_by_regions(
    these_samples_metadata = trios_meta,
    regions_bed = regions_bed,
    streamlined = FALSE,
    basic_columns = FALSE
)

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

trios_ashm_hg38 <- get_ssm_by_regions(
    these_samples_metadata = trios_meta,
    regions_bed = regions_bed,
    projection = "hg38",
    streamlined = FALSE,
    basic_columns = FALSE
)
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

setwd("~/my_dir/repos/GAMBLR.data/")

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

usethis::use_data(
    sample_data,
    overwrite = TRUE,
    compress = "xz"
)

library(data.tree)

tree <- FromListSimple(sample_data)
tree
