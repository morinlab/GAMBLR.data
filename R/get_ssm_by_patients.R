#' @title Get SSM By Patients.
#'
#' @description Get MAF-format data frame for more than one patient.
#'
#' @details This function returns variants from a set of patients.
#' This function internally calls [GAMBLR.data::get_ssm_by_samples].
#' Thus, the main contents of this function is to wrangle the provided patient IDs,
#' so that the corresponding sample IDs can be provided to the internal call of `get_ssm_by_samples`.
#' This function expects either a vector of patient IDs (`these_patients_ids`)
#' or an already subset metadata table (`these_samples_metadata`).
#'
#' @param these_patient_ids A vector of patient IDs that you want results for.
#' The user can also use a metadata table that has been subset to the patient IDs of interest (see `these_samples_metadata`).
#' @param these_samples_metadata A metadata subset to contain the rows corresponding to the patients of interest.
#' If the vector of patient IDs is missing (`these_patient_ids`), this function will default to all patient IDs in the metadata table given to this parameter.
#' @param projection Obtain variants projected to this reference (one of grch37 or hg38). Default is grch37.
#' @param this_seq_type The seq type you want results for. Default is "genome".
#' @param tool_name Optionally specify which tool to report variant from. The default is slms-3, also supports "publication" to return the exact variants as reported in the original papers.
#' @param this_study Optionally specify first name of the author for the paper
#'      from which the variants should be returned for.
#' This parameter can either be a vector of indexes (integer) or a vector of characters (matching columns in MAF).
#' @param verbose Set to FALSE to minimize the output to console. Default is TRUE. This parameter also dictates the verbosity of any helper function internally called inside the main function.
#' @param ... Any additional parameters.
#'
#' @return A data frame with SSM calls for the selected patients in MAF format.
#'
#' @import dplyr
#'
#' @export
#'
#' @examples
#'
#' # Lets find which patient_id occur more than once in the metadata first
#' my_ids = get_gambl_metadata(seq_type_filter = c("genome","capture")) %>%
#'              dplyr::group_by(patient_id) %>%
#'              dplyr::tally() %>%
#'              dplyr::filter(n>1) %>%
#'              dplyr::pull(patient_id)
#'
#' #now let's get every SSM for all samples from these patients
#' patient_maf = get_ssm_by_patients(these_patient_ids = my_ids)
#' patient_maf %>% dplyr::group_by(Tumor_Sample_Barcode) %>% 
#'                 dplyr::count() %>% head()
#'
get_ssm_by_patients = function(these_patient_ids,
                               these_samples_metadata,
                               projection = "grch37",
                               this_seq_type = "genome",
                               tool_name = "slms-3",
                               this_study,
                               verbose = FALSE,
                               ...) {

  #check if any invalid parameters are provided
  check_excess_params(...)

  #figure out what patients the user wants
  if(missing(these_patient_ids)) {
    if(missing(these_samples_metadata)) {
      stop("You must provide patient IDs (`these_patient_ids`)or a metadata
      table with the patient IDs of interest (`these_samples_metadata`)...")
    }else{
      message("No patient IDs were provided, this function will resort to
      all available patient IDs in the provided metadata.")
    }
  }else{
    if(missing(these_samples_metadata)){
      these_samples_metadata = GAMBLR.data::get_gambl_metadata(seq_type_filter =
                                                                 this_seq_type)
    }
    message("Patient IDs and metadata were provided, this function will resort to all available patient IDs in the provided metadata.")
    these_samples_metadata = these_samples_metadata %>%
      dplyr::filter(patient_id %in% these_patient_ids)
  }

  #run get_ssm_by_samples with these_samples_metadata parameter
  samples_ssm = get_ssm_by_samples(these_samples_metadata = these_samples_metadata,
                                   projection = projection,
                                   this_seq_type = this_seq_type,
                                         tool_name = tool_name,
                                         verbose = verbose,
                                         ...)

  samples_ssm = create_maf_data(samples_ssm,projection)
  # use S3-safe version of dplyr function

  samples_ssm = mutate.genomic_data(samples_ssm,maf_seq_type = this_seq_type)
}
