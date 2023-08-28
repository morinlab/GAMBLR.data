#' @title Get GAMBL Metadata.
#'
#' @description Convenience function for loading the bundled metadata, [GAMBLR.data::gambl_metadata].
#'
#' @details This bare bones function was developed to retrieve metadata for non-GSC-users.
#' For users with GSC access it is highly recommended to instead call [GAMBLR.results::get_gambl_metadata].
#' Specify the seq type (`seq_type_filter`) for the samples you want returned as the only argument.
#' It relies on the bundled metadata in this package. 
#'
#' @param seq_type_filter Specify the seq type you want to return metadata for. 
#' Default is both genome and capture (all samples).
#'
#' @return A data frame with metadata, tailored for user without GSC access.
#'
#' @import dplyr
#' 
#' @export
#'
#' @examples
#' #return metadata for genome samples
#' genome_meta = get_gambl_metadata(seq_type_filter = "genome")
#' 
#' #return metadata for capture samples
#' genome_meta = get_gambl_metadata(seq_type_filter = "capture")
#' 
#' #return metadata for genome and capture samples
#' genome_meta = get_gambl_metadata(seq_type_filter = c("genome", "capture"))
#'
get_gambl_metadata = function(seq_type_filter = c("genome", "capture")){
  message("Using the bunlded metadata in GAMBLR.data...")
  return(GAMBLR.data::gambl_metadata %>% 
           dplyr::filter(seq_type %in% seq_type_filter))
}
