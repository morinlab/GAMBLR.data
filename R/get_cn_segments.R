#' @title Get CN Segments.
#'
#' @description Retrieve all copy number segments from the GAMBL database that overlap with a single genomic coordinate range.
#'
#' @details This function returns CN segments for s specified region.
#' There are multiple ways a region can be specified.
#' For example, the user can provide the full region in a "region" format (chr:start-end) to the `region` parameter.
#' Or, the user can provide chromosome, start and end coordinates individually with `chr`, `start`, and `end` parameters.
#' For more usage examples, refer to the parameter descriptions and examples in the vignettes.
#'
#' @param region Region formatted like chrX:1234-5678 or X:1234-56789.
#' @param chromosome The chromosome you are restricting to. Required parameter if region is not specified.
#' @param qstart Start coordinate of the range you are restricting to. Required parameter if region is not specified.
#' @param qend End coordinate of the range you are restricting to. Required parameter if region is not specified.
#' @param projection Selected genome projection for returned CN segments. Default is "grch37".
#' @param this_seq_type Seq type for returned CN segments. This version of this function currently only supports "genome". See [GAMBLR.results::get_cn_segments] for capture samples.
#' @param with_chr_prefix Boolean parameter for toggling if chr prefixes should be present in the return, default is FALSE.
#' @param streamlined Return a basic rather than full MAF format. Default is FALSE.
#' @param from_flatfile This parameter does not do anything for this version of get_manta_sv. See [GAMBLR.results::get_cn_segments] for more info.
#'
#' @return A data frame with CN segments for the specified region.
#'
#' @import dplyr stringr
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
get_cn_segments = function(region,
                           chromosome,
                           qstart,
                           qend,
                           projection = "grch37",
                           this_seq_type = "genome",
                           with_chr_prefix = FALSE,
                           streamlined = FALSE,
                           from_flatfile = NULL){
  
  #check seq type
  if(this_seq_type != "genome"){
    stop("Currently, only CN segments available for genome samples (in GAMBLR.data). Please run this function with `this_seq_type` set to genome...")
  }
  
  #warn/notify the user what version of this function they are using
  message("Using the bundled CN segments (.seg) calls in GAMBLR.data...")
  
  #get invalid parameters for this function
  invalid_params = c("from_flatfile")
  
  #check if any such parameters are provided
  for(param in invalid_params){
    if(!is.null(get(param))){
      print(paste("Unsupported parameter supplied. This is only available in GAMBLR.results:", param))
      stop()
    }
  }
  
  #get valid projections
  valid_projections = grep("meta", names(GAMBLR.data::sample_data), value = TRUE, invert = TRUE)
  
  #return CN segments based on the selected projection
  if(projection %in% valid_projections){
    all_segs = GAMBLR.data::sample_data[[projection]]$seg
  }else{
    stop(paste("please provide a valid projection. The following are available:",
               paste(valid_projections,collapse=", ")))
  }
  
  #perform wrangling on the region to have it in the correct format.
  if(!missing(region)){
    region = gsub(",", "", region)
    split_chunks = unlist(strsplit(region, ":"))
    chromosome = split_chunks[1]
    startend = unlist(strsplit(split_chunks[2], "-"))
    qstart = startend[1]
    qend = startend[2]
  }
  
  #deal with chr prefixes for region, based on selected genome projection.
  if(projection == "grch37"){
    if(grepl("chr", chromosome)){
      chromosome = gsub("chr", "", chromosome)
    }
  }else{
    if(!grepl("chr", chromosome)){
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