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
    "Sex",
    "COO_consensus",
    "LymphGen",
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
    "Reddy_DLBCL" = 28985567
)


# Importing BL data from Thomas et al
# It has more patients and also contains sample ids, not just patient ids
bl_data <- list()

bl_data$meta <- read_xlsx("inst/extdata/studies/BL_Thomas.xlsx", sheet = 1)

bl_data$meta_to_bundle <- data.frame(
    bl_data$meta$`Patient barcode`,
    bl_data$meta$`Genome sample id`,
    bl_data$meta$`Genome sample id`,
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
    select(names(GAMBLR:::maf_header[1:45])) %>%
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
    names(GAMBLR:::maf_header[1:45]),
    colnames(fl_data$ssm_to_bundle)
)

new_cols <- setNames(rep(NA, length(difference)), difference)

fl_data$ssm_to_bundle  <- fl_data$ssm_to_bundle %>%
    mutate(!!! new_cols) %>%
    select(names(GAMBLR:::maf_header[1:45]))


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
    select(names(GAMBLR:::maf_header[1:45])) %>%
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

# Importing DLBCL Reddy data
reddy_data <- list()

# Importing metadata from Reddy et al and updating IDs to be consistent with GAMBL metadata
reddy_meta = readxl::read_excel("inst/extdata/studies/DLBCL_Reddy.xlsx",sheet=1) %>%
  mutate(patient_id=paste0("Reddy_",`Sample  ID`),sample_id=paste0("Reddy_",`Sample  ID`,"T")) %>%
  mutate(Tumor_Sample_Barcode=sample_id) %>%
  dplyr::rename("Sex"="Gender") %>%
  dplyr::rename("COO_consensus"="ABC GCB (RNAseq)") %>%
  mutate(COO_consensus=ifelse(COO_consensus=="Unclassified","UNCLASS",COO_consensus)) %>%
  dplyr::select(sample_id,patient_id,Tumor_Sample_Barcode,Sex,COO_consensus)

setwd("~/GAMBLR/")

reddy_meta_gambl = get_gambl_metadata(seq_type_filter="capture") %>%
  dplyr::filter(cohort == "dlbcl_reddy") %>%
  dplyr::select(sample_id,lymphgen,EBV_status_inf,cohort,pathology) %>%
  rename("LymphGen" = "lymphgen")

reddy_meta_gambl$reference_PMID = pmids$Reddy_DLBCL

reddy_data$meta_to_bundle = left_join(reddy_meta,reddy_meta_gambl)


reddy_meta_gambl$seq_type = "capture"
reddy_meta_gambl$unix_group = "icgc_dart"
reddy_meta_gambl$genome_build = "hg19-reddy"
reddy_meta_gambl$pairing_status = "unmatched"

#warning: this is very slow!
reddy_full_ssm <- get_ssm_by_samples(
  these_samples_metadata = reddy_meta_gambl
)

#restrict to the most inclusive DLBCL gene list
all_lymphoma_genes = read_tsv("inst/extdata/lymphoma_genes_comprehensive.tsv") %>%
  pull(Gene)

reddy_data$grch37$ssm_to_bundle <- dplyr::filter(reddy_full_ssm,Hugo_Symbol %in% all_lymphoma_genes)



# Importing DLBCL cell lines
cell_lines_data <- list()

cell_lines_data$meta <- get_gambl_metadata() %>%
    filter(sample_id %in% c(
        "DOHH-2", "SU-DHL-10", "OCI-Ly10", "OCI-Ly3", "SU-DHL-4"
    )) %>%
    arrange(sample_id)

cell_lines_data$meta_to_bundle <- data.frame(
    cell_lines_data$meta$patient_id,
    cell_lines_data$meta$sample_id,
    cell_lines_data$meta$sample_id,
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
    these_samples_metadata = cell_lines_data$meta
)

cell_lines_data$hg38$ssm_to_bundle <- get_ssm_by_samples(
    these_samples_metadata = cell_lines_data$meta,
    projection = "hg38"
)

cell_lines_data$grch37$cnv_to_bundle <- get_sample_cn_segments(
    sample_list = cell_lines_data$meta$sample_id,
    multiple_samples = TRUE
)

cell_lines_data$hg38$cnv_to_bundle <- get_sample_cn_segments(
    sample_list = cell_lines_data$meta$sample_id,
    multiple_samples = TRUE,
    projection = "hg38",
    with_chr_prefix = TRUE
)

cell_lines_data$grch37$sv_to_bundle <- get_manta_sv_by_samples(
    these_samples_metadata = cell_lines_data$meta,
)

cell_lines_data$hg38$sv_to_bundle <- get_manta_sv_by_samples(
    these_samples_metadata = cell_lines_data$meta,
    projection = "hg38"
)

# Adding the manta SVs for published studies
full_genome_meta = get_gambl_metadata(seq_type_filter = "genome")

bundled_meta = dplyr::filter(
    full_genome_meta,
    sample_id %in% GAMBLR.data::sample_data$meta$sample_id
)

full_sv_to_bundle = get_manta_sv_by_samples(
  these_samples_metadata = bundled_meta,projection="hg38")

annotated_sv_to_bundle = annotate_sv(full_sv_to_bundle,genome_build = "hg38")
annotated_sv_to_bundle= dplyr::filter(annotated_sv_to_bundle,!is.na(partner)) %>%
  mutate(chrom1=paste0("chr",chrom1),chrom2=paste0("chr",chrom2))

#drop all annotation columns to restore original data subset just to the putative driver SVs
annotated_sv_keep = left_join(
    full_sv_to_bundle,
    annotated_sv_to_bundle,
    by=c(
        "CHROM_A"="chrom1",
        "CHROM_B"="chrom2",
        "START_A"="start1",
        "tumour_sample_id")
    ) %>%
    dplyr::filter(!is.na(partner)) %>%
    select(c(1:16))

# Now same for the grch37 projection
full_sv_to_bundle_grch37 <- get_manta_sv_by_samples(
    these_samples_metadata = bundled_meta
)

annotated_sv_to_bundle_grch37 = annotate_sv(
    full_sv_to_bundle_grch37)

annotated_sv_to_bundle_grch37 = dplyr::filter(
    annotated_sv_to_bundle_grch37,
    !is.na(partner)
)

#drop all annotation columns to restore original data subset just to the putative driver SVs
annotated_sv_keep_grch37 = left_join(
    full_sv_to_bundle_grch37,
    annotated_sv_to_bundle_grch37,
    by=c(
        "CHROM_A"="chrom1",
        "CHROM_B"="chrom2",
        "START_A"="start1",
        "tumour_sample_id")
    ) %>%
    dplyr::filter(!is.na(partner)) %>%
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
    select(-COO_consensus, -LymphGen, -EBV_status_inf) %>%
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
    rename("LymphGen" = "lymphgen") %>%
    select(all_of(colnames_for_bundled_meta))

sample_data$meta = bind_rows(sample_data$meta,reddy_data$meta_to_bundle)

sample_data$grch37$maf <- bind_rows(
    fl_data$ssm_to_bundle,
    cell_lines_data$grch37$ssm,
    reddy_data$grch37$ssm_to_bundle
)

sample_data$hg38$maf <- bind_rows(
    bl_data$ssm_to_bundle,
    dlbcl_data$ssm_to_bundle,
    cell_lines_data$hg38$ssm
)

# add column to indicate which pipeline the data were derived from
sample_data$grch37$maf$Pipeline = "SLMS-3"
sample_data$hg38$maf$Pipeline = "SLMS-3"


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

setwd("~/my_dir/repos/GAMBLR.data/")

usethis::use_data(
    sample_data,
    overwrite = TRUE
)

library(data.tree)

tree <- FromListSimple(sample_data)
tree
