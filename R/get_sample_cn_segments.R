#' @title Get Sample CN Segments.
#'
#' @description Get all segments for a single (or multiple) sample_id(s).
#'
#' @details This function returns CN segments for samples. This works for single sample or multiple samples.
#' Specify the sample IDs you are interested in with `these_sample_ids` (as a vector of characters),
#' Or call this function with `these_samples_metadata` if you already ahve a metadata table subset to the sample IDs of interest.
#' If none off the above parameters are specified, the function will return CN segments for available samples.
#' Note, this. function internally calls [GAMBLR::id_ease] for dealing with sample IDs and metadata tables. 
#' Is this function not what you are looking for? Try one of the following, similar, functions; [GAMBLR::assign_cn_to_ssm], [GAMBLR::get_cn_segments], [GAMBLR::get_cn_states],
#'
#' @param these_sample_ids Optional argument, sample_id (vector of characters) for the sample(s) to retrieve segments for. If not provided, the function will return CN segments for all available sample IDs present in the current metadata.
#' @param these_samples_metadata Optional, provide a metadata (data frame) subset to the sample IDs of interest.
#' @param projection Selected genome projection for returned CN segments. Default is "grch37".
#' @param this_seq_type Seq type for returned CN segments. One of "genome" (default) or "capture".
#' @param with_chr_prefix Set to TRUE to add a chr prefix to chromosome names. Default is FALSE.
#' @param streamlined Return a minimal output rather than full details. Default is FALSE.
#' @param ... Any additional parameters.
#'
#' @return A data frame of segments for a specific or multiple sample ID(s).
#'
#' @import dplyr readr
#' @export
#'
#' @examples
#' #load pacakges
#' library(dplyr)
#' 
#' #get CN segments for one sample
#' dohh2_segs = get_sample_cn_segments(these_sample_ids = "DOHH-2",
#'                                     projection = "hg38", 
#'                                     streamlined = TRUE)
#'
#' #get CN segments for DLBCL cell line
#' cell_line_meta = GAMBLR.data::sample_data$meta %>% 
#'   dplyr::filter(cohort == "DLBCL_cell_lines")
#'   
#' dlbcl_segs = get_sample_cn_segments(these_samples_metadata = cell_line_meta, 
#'                                     streamlined = TRUE)
#'
get_sample_cn_segments = function(these_sample_ids,
                                  these_samples_metadata,
                                  projection = "grch37",
                                  this_seq_type = "genome",
                                  with_chr_prefix = FALSE,
                                  streamlined = FALSE,
                                  ...){
  
  #check seq type
  if(this_seq_type != "genome"){
    stop("Currently, only CN segments available for genome samples (in GAMBLR.data). Please run this function with `this_seq_type` set to genome...")
  }
  
  #get sample IDs
  meta_ids = id_ease(these_sample_ids = these_sample_ids, 
                     these_samples_metadata = these_samples_metadata, 
                     this_seq_type = this_seq_type)
  
  #subset returned list to sample IDs
  these_samples = meta_ids$these_samples
  
  #warn/notify the user what version of this function they are using
  message("Using the bundled CN segments (.seg) calls in GAMBLR.data...")
  
  #check if any invalid parameters are provided
  check_excess_params(...)
  
  #get valid projections
  valid_projections = grep("meta", names(GAMBLR.data::sample_data), value = TRUE, invert = TRUE)
  
  #return CN segments based on the selected projection
  if(projection %in% valid_projections){
    all_segs = GAMBLR.data::sample_data[[projection]]$seg %>%
      dplyr::filter(ID %in% these_samples)
  }else{
    stop(paste("please provide a valid projection. The following are available:",
               paste(valid_projections,collapse=", ")))
  }
  
  #deal with chr prefixes
  if(!with_chr_prefix){
    all_segs = all_segs %>%
      dplyr::mutate(chrom = gsub("chr", "", chrom))
  }else{
    if(!grepl("chr", all_segs$chrom[1])){
      all_segs$chrom = paste0("chr", all_segs$chrom)
    }
  }
  
  if(streamlined){all_segs = dplyr::select(all_segs, ID, CN)}
  
  return(all_segs)
}
