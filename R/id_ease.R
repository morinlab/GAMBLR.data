#' @title ID Ease
#'
#' @aliases id_ease, id ease
#'
#' @description Internal convenience function that standardize the way GAMBLR functions deals with sample IDs (these_sample_ids)
#' and metadata (these_samples_metadata).
#'
#' @details This function can take sample IDs as a vector of characters, or a metadata table in data frame format.
#' If no sample IDs are provided to the function, the function will operate on all gambl sample IDs available for the given seq type.
#' It is highly recommended to run this function with `verbose = TRUE`. 
#' Since this will not only improve the overall logic on how the function operates.
#' But also might help with debugging functions that are internally calling this function.
#' The function also performs sanity checks and notifies the user if any of the requested sample IDs are not found in the metadata.
#' In addition, the function also notifies the dimensions of the returned object, providing further insight to what is returned. 
#' As with all GAMBLR functions, providing a curated metadata table to any GAMBLR function (as opposed to a vector of IDs) is the safest way to ensure you get the expected result.
#' 
#' @param these_samples_metadata An optional data frame with metadata, subset to sample IDs of interest.
#' If not provided will retrieve GAMBL metadata for all available samples.
#' @param these_sample_ids Optional character vector of GAMBL sample IDs.
#' @param this_seq_type The seq type of interest. Default is both genome and exome, with priority for genome when a sample has >1 seq_type. 
#' @param verbose Set to FALSE to limit the information that gets printed to the console. Default is FALSE.
#' 
#' @export
#'
#' @return Metadata (data frame).
#'
#' @examples
#' #load packages
#' library(dplyr)
#' 
#' #give the function nothing (i.e return all sample IDs in the metadata for the default seq type)
#' #return metadata for all samples in the default seq type
#' all_meta = id_ease()
#'
#' #return metadata based on a sample ID
#' sample_meta = id_ease(these_sample_ids = "94-15772_tumorA")
#'
#' #return sample IDs based on an already filtered metadata
#' this_metadata = get_gambl_metadata(seq_type_filter = "genome") %>% 
#'   head(5)
#'
#' these_ids = id_ease(these_samples_metadata = this_metadata)
#'
id_ease = function(these_samples_metadata = NULL,
                   these_sample_ids = NULL,
                   this_seq_type = c("genome", "capture"),
                   verbose = FALSE){
  
  #check for provided metadata, else use GAMBL metadata
  if(is.null(these_samples_metadata)){
    if(verbose){
      message("id_ease: No metadata provided, the helper function will fetch metadata for all gambl samples in the selected seq type...") 
    }
    metadata = GAMBLR.data::get_gambl_metadata(seq_type_filter = this_seq_type)
  }else{
    if(verbose){
      message("id_ease: Metadata is provided and samples of the selected seq type are kept...") 
    }
    metadata = dplyr::filter(these_samples_metadata, seq_type %in% this_seq_type)
    not_seq_type = setdiff(these_samples_metadata$sample_id, metadata$sample_id)
    if(length(not_seq_type) > 0){
      not_seq_type_msg = gettextf("id_ease: WARNING! %i samples in the provided metadata were removed because their seq types are not the same as in the `set_type` argument.",
                                  length(not_seq_type))
      if(verbose){
        max_to_show <- 100
        if( length(not_seq_type) > max_to_show ){
          not_seq_type_msg = gettextf("%s Their first %i IDs are:", not_seq_type_msg, 
                                      max_to_show)
          not_seq_type = head(not_seq_type, max_to_show)
        }else{
          not_seq_type_msg = gettextf("%s Their IDs are:", not_seq_type_msg)
        }
        message(not_seq_type_msg)
        print(not_seq_type)
      }else{
        not_seq_type_msg = gettextf("%s Use `verbose = TRUE` to see their IDs.", not_seq_type_msg)
        message(not_seq_type_msg)
      }
    }
  }
  
  #ensure metadata is subset to specified sample IDs
  if(!is.null(these_sample_ids)){
    if(verbose){
      message("id_ease: Sample IDs are provided, filtering the metadata for selected sample IDs...") 
    }
    metadata = dplyr::filter(metadata, sample_id %in% these_sample_ids)
    
    #check if metadata is empty
    if(nrow(metadata) == 0){
      stop("No samples in the metadata, try a different sample ID...")
    }
    #check the existence of provided sample IDs in the metadata
    not_in_meta = setdiff(these_sample_ids, metadata$sample_id)
    if(length(not_in_meta) > 0){
      message("id_ease: WARNING! The following sample IDs were not found in the metadata:")
      print(not_in_meta)
    }
  }else{
    if(verbose){
      message("id_ease: No sample IDs provided, all sample IDs in the metadata will be kept...")
    }
  }
  if(verbose){
    unique_samples = unique(metadata$sample_id)
    message(paste0("id_ease: Returning metadata for ", length(unique_samples), " samples..." ))
  }
  return(metadata) 
}
