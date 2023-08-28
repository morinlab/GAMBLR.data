#' @title Get SSM By Sample.
#'
#' @description Get the SSMs (i.e. load MAF) for a single sample.
#'
#' @details This was implemented to allow flexibility because there are some samples that we may want to use a different set of variants than those in the main GAMBL merge.
#' The current use case is to allow a force_unmatched output to be used to replace the SSMs from the merge for samples with known contamination in the normal.
#' This will also be useful to apply a blacklist to individual MAFs when coupled with [GAMBLR::annotate_ssm_blacklist].
#' Is this function not what you are looking for? Try one of the following, similar, functions; [GAMBLR::get_coding_ssm], [GAMBLR::get_coding_ssm_status],
#' [GAMBLR::get_ssm_by_patients], [GAMBLR::get_ssm_by_samples], [GAMBLR::get_ssm_by_region], [GAMBLR::get_ssm_by_regions]
#'
#' @param this_sample_id Required. The sample_id you want the data from.
#' @param this_seq_type Required if not specifying these_samples_metadata. The seq_type of the sample you want data from.
#' @param these_samples_metadata Required if not specifying both this_sample_id and this_seq_type a single row or entire metadata table containing your sample_id.
#' @param tool_name This parameter does not do anything for this version of get_manta_sv. See [GAMBLR.results::get_ssm_by_sample] for more info.
#' @param projection The projection genome build. Supports hg38 and grch37.
#' @param these_genes A vector of genes to subset ssm to.
#' @param augmented This parameter does not do anything for this version of get_manta_sv. See [GAMBLR.results::get_ssm_by_sample] for more info.
#' @param flavour This parameter does not do anything for this version of get_manta_sv. See [GAMBLR.results::get_ssm_by_sample] for more info.
#' @param min_read_support Only returns variants with at least this many reads in t_alt_count (for cleaning up augmented MAFs).
#' @param basic_columns Return first 43 columns of MAF rather than full details. Default is TRUE.
#' @param maf_cols if basic_columns is set to FALSE, the user can specify what columns to be returned within the MAF. This parameter can either be a vector of indexes (integer) or a vector of characters.
#' @param verbose Enable for debugging/noisier output.
#'
#' @return data frame in MAF format.
#'
#' @import dplyr tidyr
#' 
#' @export
#'
#' @examples
#' #Basic usage:
#' this_sample_df = get_ssm_by_sample(this_sample_id = "DOHH-2",
#'                                    this_seq_type = "genome",
#'                                    projection = "grch37")
#'
#' #Get meta
#' dohh2_meta = GAMBLR.data::sample_data$meta %>% dplyr::filter(sample_id == "DOHH-2")
#'
#' #Get SSM with some filters applied:
#' dohh2_ssm = get_ssm_by_sample(this_sample_id = "DOHH-2", 
#'                               these_samples_metadata = dohh2_meta, 
#'                               these_genes = c("MYC", "BCL2"), 
#'                               verbose = FALSE,
#'                               this_seq_type = "genome",
#'                               projection = "grch37")
#'
get_ssm_by_sample = function(this_sample_id,
                             this_seq_type,
                             these_samples_metadata,
                             tool_name = NULL,
                             projection = "grch37",
                             these_genes,
                             augmented = NULL,
                             flavour = NULL,
                             min_read_support = 3,
                             basic_columns = TRUE,
                             maf_cols = NULL,
                             verbose = FALSE){
  
  #warn/notify the user what version of this function they are using
  message("Using the bundled SSM calls (.maf) calls in GAMBLR.data...")
  
  #get invalid parameters for this function
  invalid_params = c("tool_name", "augmented", "flavour")
  
  #check if any such parameters are provided
  for(param in invalid_params){
    if(!is.null(get(param))){
      print(paste("Unsupported parameter supplied. This is only available in GAMBLR.results:", param))
      stop()
    }
  }
  
  #figure out which unix_group this sample belongs to
  if(missing(these_samples_metadata)){
    these_samples_metadata = GAMBLR.data::get_gambl_metadata(seq_type_filter = this_seq_type) %>%
      dplyr::filter(sample_id == this_sample_id)
    
  }else{
    these_samples_metadata = these_samples_metadata %>%
      dplyr::filter(sample_id == this_sample_id)
  }
  
  #get valid projections
  valid_projections = grep("meta", names(GAMBLR.data::sample_data), value = TRUE, invert = TRUE)
  
  #return SSMs based on the selected projection
  if(projection %in% valid_projections){
    sample_ssm = GAMBLR.data::sample_data[[projection]]$maf %>%
      dplyr::filter(Tumor_Sample_Barcode == this_sample_id)
  }else{
    stop(paste("please provide a valid projection. The following are available:",
               paste(valid_projections,collapse=", ")))
  }

  if(min_read_support){
    # drop poorly supported reads but only from augmented MAF
    sample_ssm = dplyr::filter(sample_ssm, t_alt_count >= min_read_support)
  }
  
  if(!missing(these_genes)){
    sample_ssm = sample_ssm %>%
      dplyr::filter(Hugo_Symbol %in% these_genes)
  }
  
  #subset maf to only include first 43 columns (default)
  if(basic_columns){
    sample_ssm = dplyr::select(sample_ssm, c(1:45))
  }
  
  #subset maf to a specific set of columns (defined in maf_cols)
  if(!is.null(maf_cols) && !basic_columns){
    sample_ssm = dplyr::select(sample_ssm, all_of(maf_cols))
  }
  
  return(sample_ssm)
}
