#' @title Get CN Segments.
#'
#' @description Retrieve all copy number segments from the GAMBL outputs
#'
#' @details This function merely loads and returns all the seg_data
#' available for a projection (genome build)
#' @param these_samples_metadata User must provide a metadata table to
#' restrict the data to the samples in your table.
#' The metadata also ensures the proper handling of duplicate sample_id
#' across seq_types and ensures the seq_type in the metadata faithfully
#' represents the seq_type of the data
#' @param projection Desired genome coordinate system for returned CN segments.
#' Default is "grch37".
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
#' genome_metadata = get_gambl_metadata(seq_type_filter="genome")
#'
#' genome_segments_hg38 = get_cn_segments(
#'                              these_samples_metadata = genome_metadata,
#'                              projection="hg38")
#'
#'
get_cn_segments = function(these_samples_metadata,
                           projection = "grch37",
                           this_seq_type, ...) {
  #warn/notify the user what version of this function they are using
  message("Using the bundled CN segments (.seg) calls in GAMBLR.data...")

  #check if any invalid parameters are provided
  check_excess_params(...)

  #get valid projections
  valid_projections = grep("meta", names(GAMBLR.data::sample_data),
                           value = TRUE, invert = TRUE)

  metadata = these_samples_metadata

  sample_ids = metadata$sample_id
  #return CN segments based on the selected projection
  if (projection %in% valid_projections) {
    all_segs = GAMBLR.data::sample_data[[projection]]$seg %>%
      dplyr::filter(ID %in% sample_ids)
  }else {
    stop(paste("please provide a valid projection.",
               paste(valid_projections, collapse = ", ")))
  }

  #ensure chr prefixes are there when necessary
  if(projection == "grch37") {
    if(grepl("chr", all_segs$chrom[1])) {
      all_segs = all_segs %>%
        dplyr::mutate(chrom = gsub("chr", "", chrom))
    }
  }else {
    if (!grepl("chr",all_segs$chrom[1])) {
      all_segs = all_segs %>%
        dplyr::mutate(chrom = paste0("chr", chrom))
    }
  }

  #return S3 class with CN segments and genome_build
  all_segs = create_seg_data(all_segs, projection)
  return(all_segs)
}
