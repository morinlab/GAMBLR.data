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
#' @param min_read_support Subset returned variants to a set minimum read support. Default is 3.
#' @param basic_columns Return first 45 columns of MAF rather than full details. Default is TRUE.
#' @param maf_cols if basic_columns is set to FALSE, the user can specify what columns to be returned within the MAF. 
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
#' #load packages
#' library(dplyr)
#' 
#' #basic usage, these_patient_ids
#' my_patient = get_ssm_by_patients(these_patient_ids = "DOHH-2")
#'
#' #using a subset metadata tablee to retreive patient SSMs
#' cell_line_meta = GAMBLR.data::sample_data$meta %>% 
#'  dplyr::filter(cohort == "DLBCL_cell_lines")
#'  
#' patient_maf = get_ssm_by_patients(these_samples_metadata = cell_line_meta, 
#'                                   this_seq_type = "genome")
#'
get_ssm_by_patients = function(these_patient_ids,
                               these_samples_metadata,
                               projection = "grch37",
                               this_seq_type = "genome",
                               min_read_support = 3,
                               basic_columns = TRUE,
                               maf_cols = NULL,
                               verbose = FALSE,
                               ...){
  
  #check if any invalid parameters are provided
  check_excess_params(...)
  
  #figure out what patients the user wants
  if(missing(these_patient_ids)){
    if(missing(these_samples_metadata)){
      stop("You must provide either patient IDs (`these_patient_ids`) or a metadata table with the patient IDs of interest (`these_samples_metadata`)...")
    }else{
      message("No patient IDs were provided, this function will resort to all available patient IDs in the provided metadata.")
    }
  }else{
    if(missing(these_samples_metadata)){
      these_samples_metadata = GAMBLR.data::get_gambl_metadata(seq_type_filter = this_seq_type)
    }
    message("Patient IDs and metadata were provided, this function will resort to all available patient IDs in the provided metadata.")
    these_samples_metadata = these_samples_metadata %>%
      dplyr::filter(patient_id %in% these_patient_ids)
  }
  
  #run get_ssm_by_samples with these_samples_metadata parameter
  return(GAMBLR.data::get_ssm_by_samples(these_samples_metadata = these_samples_metadata,
                                         projection = projection,
                                         this_seq_type = this_seq_type,
                                         min_read_support = min_read_support,
                                         basic_columns = basic_columns,
                                         maf_cols = maf_cols,
                                         verbose = verbose,
                                         ...))
}
