#' @title Get GAMBL Metadata.
#'
#' @description Convenience function for loading the sample metadata.
#'
#' @details This bare bones function was developed to retrieve metadata for
#' non-GSC-users. Specify the seq type (`seq_type_filter`) for the samples you
#' want returned as the only argument.
#' It relies on the bundled metadata in this package.
#' Specify `case_set` argument to retreive samples from particular study.
#' Currently supported case_sets are: FL_Dreval (FL samples from Dreval et al),
#' DLBCL_Dreval (DLBCL samples from Dreval et al), FL-DLBCL-study (all samples
#' from Dreval et al), DLBCL_Arthur (all samples from Arthur et al study),
#' DLBCL_Hilton (all samples from Hilton et al DLBCL Trios study),
#' DLBCL_cell_lines (5 DLBCL cell lines), DLBCL_Chapuy (all samples from Chapuy
#' et al study), DLBCL_Schmitz (all samples from Schmitz et al study),
#' DLBCL_Reddy (all samples from Reddy et al study), DLBCL_Thomas (HTMCP DLBCLs
#' from Thomas et al study), BL_Thomas (BL samples from Thomas et al study)
#'
#' @param seq_type_filter Specify the seq type you want to return metadata for.
#' Default is "genome".
#' @param case_set Optionally specify study details to return samples from a
#' particular case set. See function description for supported case sets.
#' @param ... Any additional parameters.
#'
#' @return A data frame with metadata, tailored for user without GSC access.
#'
#' \describe{
#'   \item{compression}{Format of the original data used as input for our analysis pipelines (cram, bam or fastq)}
#'   \item{bam_available}{Whether or not this file was available when last checked.}
#'   \item{patient_id}{The anonymized unique identifier for this patient. For BC samples, this will be Res ID.}
#'   \item{sample_id}{A unique identifier for the sample analyzed.}
#'   \item{seq_type}{The assay type used to produce this data (one of "genome","capture, "mrna", "promethION")}
#'   \item{genome_build}{The name of the genome reference the data were aligned to.}
#'   \item{cohort}{Name for a group of samples that were added together (usually from a single study), often in the format {pathology}_{cohort_descriptor}.}
#'   \item{pathology}{The diagnosis or pathology for the sample}
#'   \item{time_point}{Timing of biopsy in increasing alphabetical order (A = diagnosis, B = first relapse etc)}
#'   \item{ffpe_or_frozen}{Whether the nucleic acids were extracted from a frozen or FFPE sample}
#'   \item{COO_consensus}{Consensus call of COO between different sources.}
#'   \item{DHITsig_consensus}{Consensus call of DHIT signature status between different sources.}
#'   \item{EBV_status_inf}{Inferred EBV status of the tumor}
#'   \item{lymphgen_no_cnv}{LymphGen label using model without CNV}
#'   \item{lymphgen_with_cnv}{LymphGen label using model with CNV}
#'   \item{lymphgen_cnv_noA53}{LymphGen label using model with CNV but excluding A53 class}
#'   \item{lymphgen_wright}{The LymphGen call for this sample from Wright et all (if applicable)}
#'   \item{fl_grade}{Grade of FL samples}
#'   \item{normal_sample_id}{Sample id for normal tissue used in the analysis}
#'   \item{pairing_status}{Matching status of the sample}
#'   \item{lymphgen}{LymphGen label}
#'   \item{molecular_BL}{label of the sample according to the molecular BL classifier}
#'   \item{Tumor_Sample_Barcode}{Duplicate of sample_id for simplifying joins to MAF data frames}
#'   \item{pathology_rank}{Numeric rank for consistent ordering of samples by pathology}
#'   \item{hiv_status}{HIV status of the sample}
#'   \item{age_group}{Adult_BL or Pediatric_BL or Other, specific to the BLGSP study}
#'   \item{sex}{The biological sex of the patient, if available. Allowable options: M, F, NA}
#' }
#'
#' @import dplyr purrr
#'
#' @export
#'
#' @examples
#' #return metadata for genome samples
#' genome_meta = get_gambl_metadata(seq_type_filter = "genome")
#'
#' #return metadata for capture samples
#' capture_meta = get_gambl_metadata(seq_type_filter = "capture")
#'
#' #return metadata for genome and capture samples
#' all_meta = get_gambl_metadata(seq_type_filter = c("genome", "capture"))
#'
get_gambl_metadata = function(
    seq_type_filter = "genome",
    case_set,
    ...
){

    #check if any invalid parameters are provided
    check_excess_params(...)

    message("Using the bundled metadata in GAMBLR.data...")
    metadata <- GAMBLR.data::sample_data$meta %>%
            dplyr::filter(seq_type %in% seq_type_filter)


    if(!missing(case_set)){

        # pre-defined case sets
        if(case_set == "FL_Dreval"){
            metadata <- metadata %>%
                dplyr::filter(cohort == "FL_Dreval", pathology == "FL")
        }else if(case_set == "DLBCL_Dreval"){
            metadata <- metadata %>%
                dplyr::filter(cohort == "FL_Dreval", pathology == "DLBCL")
        }else if(case_set == "FL-DLBCL-study"){
            metadata <- metadata %>%
                dplyr::filter(cohort == "FL_Dreval")
        }else if(case_set == "DLBCL_Arthur"){
            metadata <- metadata %>%
                dplyr::filter(cohort == "DLBCL_Arthur")
        }else if(case_set == "DLBCL_Hilton"){
            metadata <- metadata %>%
                dplyr::filter(cohort == "DLBCL_Hilton")
        }else if(case_set == "DLBCL_cell_lines"){
            metadata <- metadata %>%
                dplyr::filter(cohort == "DLBCL_cell_lines")
        }else if(case_set == "DLBCL_Chapuy"){
            metadata <- metadata %>%
                dplyr::filter(cohort == "dlbcl_chapuy")
        }else if(case_set == "DLBCL_Schmitz"){
            metadata <- metadata %>%
                dplyr::filter(cohort == "dlbcl_schmitz")
        }else if(case_set == "DLBCL_Reddy"){
            metadata <- metadata %>%
                dplyr::filter(cohort == "dlbcl_reddy")
        }else if(case_set == "BL_Thomas"){
            metadata <- metadata %>%
                dplyr::filter(cohort == "BL_Thomas")
        }else if(case_set == "DLBCL_Thomas"){
            metadata <- metadata %>%
                dplyr::filter(cohort == "DLBCL_Thomas")
        }else{
            message(paste("case set", case_set, "not available"))
            return()
        }
    }

    metadata <- metadata %>%
        dplyr::left_join(
            gambl_metadata,
            by = "sample_id",
            suffix = c(".X", ".Y")
        ) %>%
        split.default(gsub('.[XY]', '', names(.))) %>%
        purrr::map_dfc( ~ if (ncol(.x) == 1)
            .x
            else
            dplyr::mutate(.x,!!sym(gsub('.X', '', names(
                .x
            )[1])) := dplyr::coalesce(!!!syms(names(
                .x
            ))))) %>%
        dplyr::select(!contains("."))
    #ensure only unique rows are returned
    return(unique(metadata))
}
