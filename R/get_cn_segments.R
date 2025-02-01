## GAMBLR.data
#' Create Segmented Data
#'
#' This function creates segmented data from the given input.
#'
#' @param seg_df A data frame containing the segmented data.
#' @param genome_build A string specifying the genome build ("grch37" or "hg38").
#' @return A data frame with class attributes for segmented data.
#' @export
#' @examples
#' seg_df <- data.frame(...)
#' create_seg_data(seg_df, "grch37")
create_seg_data <- function(seg_df, genome_build) {
  if (!inherits(seg_df, "data.frame")) stop("data must be a data frame")
  if (!genome_build %in% c("grch37", "hg38")) stop("Invalid genome build")
  
  structure(seg_df, 
            class = c("seg_data", class(seg_df)), 
            genome_build = genome_build)
}

#' @title Get CN Segments.
#'
#' @description Retrieve all copy number segments from the GAMBL outputs
#'
#' @details This function merely loads and returns all the seg_data available for a projection (genome build)
#' @param these_samples_metadata User must provide a metadata table to restrict the data to the samples in your table. 
#' The metadata also ensures the proper handling of duplicate sample_id across seq_types and ensures the 
#' seq_type in the metadata faithfully represents the seq_type of the data
#' @param projection Desired genome coordinate system for returned CN segments. Default is "grch37".
#' @param this_seq_type Deprecated.
#' @param ... Additional parameters to be passed to the function.
#'
#' @return A data frame with CN segments for the specified region.
#'
#' @import dplyr
#' @export
#'
#' @examples
#' # Example for the capture samples:
#' 
#' genome_metadata = GAMBLR.data::get_gambl_metadata(seq_type_filter="genome") 
#'                       
#' genome_segments_hg38 = get_cn_segments(
#'                              these_samples_metadata = genome_metadata,
#'                              projection="hg38")
#'
#'
get_cn_segments = function(these_samples_metadata,
                           projection = "grch37",
                           this_seq_type,...){
  #warn/notify the user what version of this function they are using
  message("Using the bundled CN segments (.seg) calls in GAMBLR.data...")

  #check if any invalid parameters are provided
  check_excess_params(...)

  #get valid projections
  valid_projections = grep("meta", names(GAMBLR.data::sample_data), value = TRUE, invert = TRUE)

  metadata = these_samples_metadata

  sample_ids = metadata$sample_id
  #return CN segments based on the selected projection
  if(projection %in% valid_projections){
    all_segs = GAMBLR.data::sample_data[[projection]]$seg %>%
      dplyr::filter(ID %in% sample_ids)
  }else{
    stop(paste("please provide a valid projection. The following are available:",
               paste(valid_projections,collapse=", ")))
  }

  #ensure chr prefixes are there when necessary 
  if(projection=="grch37"){
    if(grepl("chr",all_segs$chrom[1])){
      all_segs = all_segs %>%
        dplyr::mutate(chrom = gsub("chr", "", chrom))
    }
  }else{
    if(!grepl("chr",all_segs$chrom[1])){
      all_segs = all_segs %>%
        dplyr::mutate(chrom = paste0("chr", chrom))
    }
  }

  #return S3 class with CN segments and genome_build 
  all_segs = create_seg_data(all_segs,projection)
  return(all_segs)
}
