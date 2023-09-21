#' @title Get SSM By Samples.
#'
#' @description Get the SSMs (i.e. load MAF) for a single sample or a collection fo samples.
#'
#' @details Retrieve a maf for a specific sample or a set of samples. 
#' Either specify the sample IDs of interest with `these_sample_ids`.
#' Or a metadata table subset to the sample IDs of interest.
#'
#' @param these_sample_ids The sample_id you want the data from.
#' @param these_samples_metadata Required if not specifying this_sample_id.
#' @param this_seq_type Required if not specifying these_samples_metadata. The seq_type of the sample you want data from.
#' @param projection The projection genome build. Supports hg38 and grch37.
#' @param these_genes A vector of genes to subset ssm to.
#' @param min_read_support Only returns variants with at least this many reads in t_alt_count (for cleaning up augmented MAFs).
#' @param basic_columns Return first 43 columns of MAF rather than full details. Default is TRUE.
#' @param maf_cols if basic_columns is set to FALSE, the user can specify what columns to be returned within the MAF. This parameter can either be a vector of indexes (integer) or a vector of characters.
#' @param verbose Enable for debugging/noisier output.
#' @param ... Any additional parameters.
#'
#' @return data frame in MAF format.
#'
#' @import dplyr
#' 
#' @export
#'
#' @examples
#' #load packages
#' library(dplyr)
#' 
#' #return a MAF for DLBCL cell line
#' cell_line_meta = GAMBLR.data::sample_data$meta %>% 
#'   dplyr::filter(cohort == "DLBCL_cell_lines")
#'   
#' dlbcl_maf = get_ssm_by_samples(these_samples_metadata = cell_line_meta)
#' 
get_ssm_by_samples <- function(these_sample_ids,
                               these_samples_metadata,
                               this_seq_type = "genome",
                               projection = "grch37",
                               these_genes,
                               min_read_support = 3,
                               basic_columns = TRUE,
                               maf_cols = NULL,
                               verbose = FALSE,
                               ...){
  
  #check seq type
  if(this_seq_type != "genome"){
    stop("Currently, SSM results are only available for genome samples (in GAMBLR.data). Please run this function with `this_seq_type` set to genome...")
  }
  
  #warn/notify the user what version of this function they are using
  message("Using the bundled SSM calls (.maf) calls in GAMBLR.data...")
  
  #check if any invalid parameters are provided
  check_excess_params(...)
  
  #get sample IDs
  meta_ids = id_ease(these_sample_ids = these_sample_ids, 
                     these_samples_metadata = these_samples_metadata, 
                     this_seq_type = "genome",
                     verbose = verbose) #currently, only genome samples are bundled in the repo.
  
  #extract sample IDs
  these_ids = meta_ids$sample_id
  
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

  #drop poorly supported reads
  sample_ssm = dplyr::filter(sample_ssm, t_alt_count >= min_read_support)
  
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


#' @rdname get_ssm_by_samples
#' 
#' @examples
#' #basic usage, using a sample ID
#' dohh2_maf = get_ssm_by_samples(these_sample_ids = "DOHH-2")
#' 
#' @export
#' 
get_ssm_by_sample <- get_ssm_by_samples
