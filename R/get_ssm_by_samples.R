#' @title Get SSM By Samples.
#'
#' @description Get the SSMs (i.e. load MAF) for a single sample or a collection fo samples.
#' 
#' @aliases get_ssm_by_sample
#'
#' @details Retrieve a maf for a specific sample or a set of samples. 
#' Either specify the sample IDs of interest with `these_sample_ids`.
#' Or a metadata table subset to the sample IDs of interest.
#'
#' @param these_sample_ids The sample_id you want the data from.
#' @param these_samples_metadata Required if not specifying this_sample_id.
#' @param this_seq_type Required if not specifying these_samples_metadata. The seq_type of the sample you want data from.
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
#' #Basic usage, using a sample ID
#' dohh2_maf = get_ssm_by_sample(these_sample_ids = "DOHH-2")
#' 
#' #Return a MAF for DLBCL cell line
#' cell_line_meta = GAMBLR.data::sample_data$meta %>% dplyr::filter(cohort == "DLBCL_cell_lines")
#' dlbcl_maf = get_ssm_by_samples(these_samples_metadata = cell_line_meta)
#' 
get_ssm_by_samples <- function(these_sample_ids,
                              these_samples_metadata,
                              this_seq_type = "genome",
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
  
  #get sample IDs
  meta_ids = id_ease(these_sample_ids = these_sample_ids, 
                     these_samples_metadata = these_samples_metadata, 
                     this_seq_type = this_seq_type)
  
  #extract sample IDs
  these_ids = meta_ids$these_samples
  
  #get valid projections
  valid_projections = grep("meta", names(GAMBLR.data::sample_data), value = TRUE, invert = TRUE)
  
  #return SSMs based on the selected projection
  if(projection %in% valid_projections){
    sample_ssm = GAMBLR.data::sample_data[[projection]]$maf %>%
      dplyr::filter(Tumor_Sample_Barcode %in% these_ids)
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
