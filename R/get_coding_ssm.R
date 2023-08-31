#' @title Get Coding SSMs
#'
#' @description Convenience function for loading coding Simple Somatic Mutations (SSM) from the bundled data [GAMBLR.data::sample_data].
#'
#' @details This "bare bones" function was developed to retrieve coding SSM calls for non-GSC-users.
#' For users with GSC access it is highly recommended to instead call [GAMBLR.results::get_coding_ssm].
#' Effectively retrieve coding SSM calls. Multiple filtering parameters are available for this function.
#' For more information on how to implement the filtering parameters, refer to the parameter descriptions as well as examples in the vignettes.
#' It relies on the bundled sample data in this package. 
#'
#' @param limit_cohort Supply this to restrict mutations to one or more cohorts in a vector.
#' @param exclude_cohort  Supply this to exclude mutations from one or more cohorts in a vector.
#' @param limit_pathology Supply this to restrict mutations to one pathology.
#' @param limit_samples Supply this to restrict mutations to a vector of sample_id (instead of subsetting using the provided metadata)
#' @param these_samples_metadata Supply a metadata table to auto-subset the data to samples in that table before returning
#' @param force_unmatched_samples Optional argument for forcing unmatched samples, using [GAMBLR.data::get_ssm_by_samples].
#' @param projection Reference genome build for the coordinates in the MAF file. The default is hg19 genome build.
#' @param seq_type The seq_type you want back, default is genome.
#' @param basic_columns Set to FALSE to override the default behavior of returning only the first 45 columns of MAF data.
#' @param maf_cols if basic_columns is set to FALSE, the user can specify what columns to be returned within the MAF. This parameter can either be a vector of indexes (integer) or a vector of characters (matching columns in MAF).
#' @param min_read_support Only returns variants with at least this many reads in t_alt_count (for cleaning up augmented MAFs)
#' @param include_silent Logical parameter indicating whether to include silent mutations into coding mutations. Default is TRUE.
#' @param verbose Set to FALSE to minimize the output to console. Default is TRUE. This parameter also dictates the verbosity of any helper function internally called inside the main function.
#' @param ... Any additional parameters.
#'
#' @return A data frame (.maf) with coding SSMs in reference to the selected projection.
#' 
#' @export
#' 
#' @examples
#' #load pacakges
#' library(dplyr)
#' 
#' #return SSMs in reference to GRCh37:
#' ssm_grch37 = get_coding_ssm(seq_type = "genome")
#'
#' #return SSMs in reference to hg38:
#' cell_line_meta = GAMBLR.data::sample_data$meta %>% 
#'   dplyr::filter(cohort == "DLBCL_cell_lines")
#' 
#' ssm_hg38 = get_coding_ssm(projection = "hg38", 
#'                           these_samples_metadata = cell_line_meta, 
#'                           seq_type = "genome")
#' 
get_coding_ssm = function(limit_cohort,
                          exclude_cohort,
                          limit_pathology,
                          limit_samples,
                          these_samples_metadata,
                          force_unmatched_samples,
                          projection = "grch37",
                          seq_type = "genome",
                          basic_columns = TRUE,
                          maf_cols = NULL,
                          min_read_support = 3,
                          include_silent = TRUE,
                          verbose = FALSE,
                          ...){
  
  #warn/notify the user what version of this function they are using
  message("Using the bundled SSM calls (.maf) calls in GAMBLR.data...")
  
  #check seq type
  if(seq_type != "genome"){
    stop("Currently, SSM results are only available for genome samples (in GAMBLR.data). Please run this function with `this_seq_type` set to genome...")
  }
  
  #check if any invalid parameters are provided
  check_excess_params(...)
  
  #get valid projections
  valid_projections = grep("meta", names(GAMBLR.data::sample_data), value = TRUE, invert = TRUE)
  
  #return SSMs based on the selected projection
  if(projection %in% valid_projections){
    muts = GAMBLR.data::sample_data[[projection]]$maf
  }else{
    stop(paste("please provide a valid projection. The following are available:",
               paste(valid_projections,collapse=", ")))
  }
  
  if(!include_silent){
    coding_class = coding_class[coding_class != "Silent"]
  }
  
  #get samples with the dedicated helper function
  all_meta = id_ease(these_samples_metadata = these_samples_metadata,
                     verbose = verbose,
                     this_seq_type = "genome") #currently, only genome samples are bundled in the repo.
  
  #limit cohort
  if(!missing(limit_cohort)){
    all_meta = all_meta %>%
      dplyr::filter(cohort %in% limit_cohort)
  }
  
  #exclude cohort
  if(!missing(exclude_cohort)){
    all_meta = all_meta %>%
      dplyr::filter(!cohort %in% exclude_cohort)
  }
  
  #limit pathology
  if(!missing(limit_pathology)){
    all_meta = all_meta %>%
      dplyr::filter(pathology %in% limit_pathology)
  }
  
  #limit samples
  if(!missing(limit_samples)){
    all_meta = all_meta %>%
      dplyr::filter(sample_id %in% limit_samples)
  }
  
  #pull info for loading .CDS.maf
  sample_ids = pull(all_meta, sample_id)
  
  #drop variants with low read support (default is 3), enforce sample IDs and keep only coding variants
  muts = dplyr::filter(muts, t_alt_count >= min_read_support) %>%
    dplyr::filter(Tumor_Sample_Barcode %in% sample_ids) %>%
    dplyr::filter(Variant_Classification %in% coding_class)
  
  #filter maf on selected sample ids
  muts = muts %>%
    dplyr::filter(Tumor_Sample_Barcode %in% sample_ids)
  
  mutated_samples = length(unique(muts$Tumor_Sample_Barcode))
  message(paste("after linking with metadata, we have mutations from", mutated_samples, "samples"))
  
  #drop rows for these samples so we can swap in the force_unmatched outputs instead
  if(!missing(force_unmatched_samples)){
    muts = muts %>%
      dplyr::filter(!Tumor_Sample_Barcode %in% force_unmatched_samples)
    
    nsamp = length(force_unmatched_samples)
    message(paste("dropping variants from", nsamp, "samples and replacing with force_unmatched outputs"))
    
    #get replacements using get_ssm_by_samples
    fu_muts = GAMBLR.data::get_ssm_by_samples(these_sample_ids = force_unmatched_samples)
    muts = bind_rows(muts, fu_muts)
  }
  
  #subset maf to a specific set of columns (defined in maf_cols)
  if(!is.null(maf_cols) && !basic_columns){
    muts = dplyr::select(muts, all_of(maf_cols))
  }
  
  return(muts)
}
