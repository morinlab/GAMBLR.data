#' @title Collate Results
#'
#' @description Bring together collated results for a selection of gambl samples.
#'
#' @details Currently, this function only gathers QC metrics (`mirage_metrics`) as the only collated result.
#' Potentially, in the future, additional collated results can be added by this function as well.  
#'
#' @param sample_table A vector of characters with sample IDs, or a data frame with sample IDs in a column (sample_id). 
#' If provided, this will overwrite any sample subsets provided these_samples_metadata.
#' @param these_samples_metadata A metadata table with sample IDs of interest. 
#' If not provided, the function will get metadata for all available samples. 
#' This parameter is intended to use in combination with `join_with_full_metadata`.
#' @param join_with_full_metadata Set to TRUE to horizontally expand metadata with QC results. 
#' Default is FALSE. If `these_samples_metadata` is provided, collated resutls will be added to this metadata table.
#' If not provided, the function will join collated results with all available metadata in the specified seq_type (`seq_type_filter`).
#' @param seq_type_filter Filtering criteria for `get_gambl_metadata` if `these_samples_metadata` is not provided, default is genomes and captures. 
#' @param ... Any additional parameters.
#'
#' @return A data frame with collated results.
#'
#' @import dplyr
#' 
#' @export
#'
#' @examples
#' #load packages
#' library(dplyr)
#' 
#' #return collated results for all available samples
#' all_collated = collate_results()
#'
#' #return available collated results for a metadata subset
#' fl_collated = collate_results(
#'  these_samples_metadata = get_gambl_metadata(
#'    seq_type_filter = "genome") %>% 
#'    dplyr::filter(pathology == "FL"))
#'
#' #horizontally expand a metadata subset with collated results
#' fl_meta_collated = collate_results(
#'  join_with_full_metadata = TRUE, 
#'  these_samples_metadata = get_gambl_metadata(
#'    seq_type_filter = "genome") %>% 
#'    dplyr::filter(pathology == "FL"))
#'
#' #horizontally expand all available metadata with collated results
#' all_meta_collated = collate_results(join_with_full_metadata = TRUE)
#' 
collate_results = function(sample_table,
                           these_samples_metadata,
                           join_with_full_metadata = FALSE,
                           seq_type_filter = c("genome", "capture"),
                           ...){
  
  #check if any invalid parameters are provided
  check_excess_params(...)
  
  #warn/notify the user what version of this function they are using
  message("Using the bundled collated results in GAMBLR.data...")
  
  if(missing(these_samples_metadata)){
    these_samples_metadata = get_gambl_metadata(seq_type_filter = seq_type_filter)
  }
  
  if(missing(sample_table)){
    sample_table = these_samples_metadata %>% 
      pull(sample_id)
  }else{
    if(is.data.frame(sample_table)){
      sample_table = sample_table$sample_id
    }
  }
  
  #read mirage metrics and subset to the sample IDs (in sample_table) we have QC data for
  collated = GAMBLR.data::mirage_metrics %>%
    dplyr::filter(sample_id %in% sample_table)

  #horizontally expand the provided metadata with QC results
  if(join_with_full_metadata){
    full_table = left_join(these_samples_metadata, collated)
    return(full_table)
  }
  return(collated)
}
