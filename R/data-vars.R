#' @title Coding Class
#' 
#' @description Coding classes as a vector of characters.
#' 
#' @details This object is called by [GAMBLR.data::get_coding_ssm] and [GAMBLR.data::assign_cn_to_ssm].
#' 
#' @noRd
#' 
coding_class = c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", 
                 "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation", 
                 "Nonstop_Mutation", "Silent", "Splice_Region", 
                 "Splice_Site", "Targeted_Region", "Translation_Start_Site")

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    ".", ":=", "CHROM_A", "CHROM_B", "CN", "Chromosome", 
    "End_Position", "End_Position.x", "FILTER", "Gene", 
    "HGVSp_Short", "Hugo_Symbol", "ID", "LOH", "LOH_flag", 
    "SCORE", "START_A", "START_B", "Start_Position", "Start_Position.x", 
    "Tumor_Sample_Barcode", "VAF_tumour", "Variant_Classification", 
    "bin", "category", "chrom", "cohort", "colour", "curated", "end", 
    "ensembl_gene_id", "gambl_metadata", "gene", "genome_build", 
    "grch37_ashm_regions", "group", "head", "hg38_ashm_regions", 
    "hot_spot", "hotspot_regions_grch37", "hotspot_regions_hg38", 
    "is_alias", "log.ratio", "mutated", "mutation_count", "n_mut", 
    "name", "pair_status", "pathology", "patient_id", "region", 
    "row_id", "sample_id", "seq_type", "start", "t_alt_count", 
    "tumour_sample_id", "window_end", "window_start"
  ))
}
