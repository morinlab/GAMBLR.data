#' @title Get CN Segments.
#'
#' @description Retrieve all copy number segments from the GAMBL database that overlap with a single genomic coordinate range.
#'
#' @details This function returns CN segments for s specified region.
#' There are multiple ways a region can be specified.
#' For example, the user can provide the full region in a "region" format (chr:start-end) to the `region` parameter.
#' Or, the user can provide chromosome, start and end coordinates individually with `chr`, `qstart`, and `qend` parameters.
#' For more usage examples, refer to the parameter descriptions and examples in the vignettes.
#'
#' @param these_sample_ids Optional, a vector of multiple sample_id (or a single sample ID as a string) that you want results for.
#' @param these_samples_metadata Optional, a metadata table (with sample IDs in a column) to subset the return to.
#' If not provided (and if `these_sample_ids` is not provided), the function will return all samples from the specified seq_type in the metadata.
#' @param region Region formatted like chrX:1234-5678 or X:1234-56789.
#' @param chromosome The chromosome you are restricting to. Required parameter if region is not specified.
#' @param qstart Start coordinate of the range you are restricting to. Required parameter if region is not specified.
#' @param qend End coordinate of the range you are restricting to. Required parameter if region is not specified.
#' @param projection Selected genome projection for returned CN segments. Default is "grch37".
#' @param this_seq_type Seq type for returned CN segments. Default is genome.
#' @param with_chr_prefix Boolean parameter for toggling if chr prefixes should be present in the return, default is FALSE.
#' @param streamlined Return a basic rather than full MAF format. Default is FALSE.
#' @param ... Any additional parameters.
#'
#' @return A data frame with CN segments for the specified region.
#'
#' @import dplyr
#'
#' @export
#'
#' @examples
#' #Example using chromosome, qstart and qend parameters:
#' segments_region_grch37 = get_cn_segments(chromosome = "chr8",
#'                                          qstart = 128723128,
#'                                          qend = 128774067)
#'
#' #Example using the regions parameter:
#' segments_region_hg38 = get_cn_segments(region = "chr8:128,723,128-128,774,067",
#'                                        projection = "hg38",
#'                                        with_chr_prefix = TRUE)
#'
get_cn_segments = function(these_sample_ids = NULL,
                           these_samples_metadata = NULL,
                           region,
                           chromosome,
                           qstart,
                           qend,
                           projection = "grch37",
                           this_seq_type = "genome",
                           with_chr_prefix = FALSE,
                           streamlined = FALSE,
                           ...){

  #warn/notify the user what version of this function they are using
  message("Using the bundled CN segments (.seg) calls in GAMBLR.data...")

  #check if any invalid parameters are provided
  check_excess_params(...)

  #get valid projections
  valid_projections = grep("meta", names(GAMBLR.data::sample_data), value = TRUE, invert = TRUE)

  #get samples with the dedicated helper function
  metadata = id_ease(these_samples_metadata = these_samples_metadata,
                     these_sample_ids = these_sample_ids,
                     this_seq_type = this_seq_type)

  sample_ids = metadata$sample_id

  #return CN segments based on the selected projection
  if(projection %in% valid_projections){
    all_segs = GAMBLR.data::sample_data[[projection]]$seg %>%
      dplyr::filter(ID %in% sample_ids)
  }else{
    stop(paste("please provide a valid projection. The following are available:",
               paste(valid_projections,collapse=", ")))
  }

  #perform wrangling on the region to have it in the correct format.
  if(!missing(region)){
    if(length(region) > 1){
      stop("You are providing more than one region...")
    }
    region = gsub(",", "", region)
    split_chunks = unlist(strsplit(region, ":"))
    chromosome = split_chunks[1]
    startend = unlist(strsplit(split_chunks[2], "-"))
    qstart = startend[1]
    qend = startend[2]
  }else{
    if(missing(chromosome)){
      stop("You have not provided a region, or a region in an acceptable format..")
    }else{
      if(length(chromosome) > 1){
        stop("You are providing more than one region...")
      }
    }
  }

  #deal with chr prefixes for region, based on selected genome projection.
  if(projection == "grch37"){
    if(all(grepl("chr", chromosome))){
      chromosome = gsub("chr", "", chromosome)
    }
  }else{
    if(all(!grepl("chr", chromosome))){
      chromosome = paste0("chr", chromosome)
    }
  }

  #enforce data type for qend and qstart coordiantes.
  qstart = as.numeric(qstart)
  qend = as.numeric(qend)

  all_segs = all_segs %>%
    dplyr::filter((chrom == chromosome & start <= qstart & end >= qend) | (chrom == chromosome & start >= qstart & end <= qend)) %>%
    as.data.frame()

  #deal with chr prefixes
  if(!with_chr_prefix){
    if(all(str_detect(all_segs$chrom, "chr"))){
      all_segs = all_segs %>%
        dplyr::mutate(chrom = gsub("chr", "", chrom))
    }
  }else{
    if(all(!str_detect(all_segs$chrom, "chr"))){
      all_segs = all_segs %>%
        dplyr::mutate(chrom = paste0("chr", chrom))
    }
  }

  #subset to only a few columns with streamlined = TRUE.
  if(streamlined){
    all_segs = dplyr::select(all_segs, ID, CN)
  }

  #return data frame with CN segments
  return(all_segs)
}
