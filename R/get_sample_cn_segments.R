#' @title Get Sample CN Segments.
#'
#' @description Get all segments for a single (or multiple) sample_id(s).
#'
#' @details This function returns CN segments for samples. This works for single sample or multiple samples.
#' For multiple samples, remember to set the Boolean parameter `multiple_samples = TRUE` and give the `sample_lsit` a vector of characters with one sample ID per row.
#' For more information on how this function can be run in different ways, refer to the parameter descriptions, examples and vignettes.
#' Is this function not what you are looking for? Try one of the following, similar, functions; [GAMBLR::assign_cn_to_ssm], [GAMBLR::get_cn_segments], [GAMBLR::get_cn_states],
#'
#' @param this_sample_id Optional argument, single sample_id for the sample to retrieve segments for.
#' @param multiple_samples Set to TRUE to return cn segments for multiple samples specified in `samples_list` parameter. Default is FALSE.
#' @param sample_list Optional vector of type character with one sample per row, required if multiple_samples is set to TRUE.
#' @param from_flatfile This parameter does not do anything for this version of get_manta_sv. See [GAMBLR.results::get_sample_cn_segments] for more info.
#' @param projection Selected genome projection for returned CN segments. Default is "grch37".
#' @param this_seq_type Seq type for returned CN segments. One of "genome" (default) or "capture".
#' @param with_chr_prefix Set to TRUE to add a chr prefix to chromosome names. Default is FALSE.
#' @param streamlined Return a minimal output rather than full details. Default is FALSE.
#'
#' @return A data frame of segments for a specific or multiple sample ID(s).
#'
#' @import dplyr readr
#' 
#' @export
#'
#' @examples
#' #return sample CNN segmetns for one sample
#' my_segs = get_sample_cn_segments(this_sample_id = "DOHH-2", projection = "hg38", streamlined = TRUE)
#'
#' #define some samples
#' my_samples = c("SU-DHL-4", "DOHH-2")
#'
#' #return CN segments for thhese samples
#' these_segs = get_sample_cn_segments(multiple_samples = TRUE,
#'                                     sample_list = my_samples)
#'
get_sample_cn_segments = function(this_sample_id,
                                  multiple_samples = FALSE,
                                  sample_list,
                                  from_flatfile = NULL,
                                  projection = "grch37",
                                  this_seq_type = "genome",
                                  with_chr_prefix = FALSE,
                                  streamlined = FALSE){
  
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
  
  if(!missing(this_sample_id) & !multiple_samples){
    all_segs = dplyr::filter(all_segs, ID %in% this_sample_id)
  }else if(!missing(sample_list)){
    all_segs = dplyr::filter(all_segs, ID %in% sample_list)
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
  
  if(streamlined){
    all_segs = dplyr::select(all_segs, ID, CN)
    }
  
  return(all_segs)
}
