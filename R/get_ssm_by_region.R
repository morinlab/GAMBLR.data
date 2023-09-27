#' @title Get SSM By Region.
#'
#' @description Retrieve all SSMs from the GAMBL database within a single genomic coordinate range.
#'
#' @details This function lets the user specify a region of interest for returning SSM calls within that region.
#' There are multiple ways a region can be specified. For example, the user can provide the full region in a "region" format (chr:start-end) to the `region` parameter.
#' Or, the user can provide chromosome, start and end coordinates individually with `chr`, `start`, and `end` parameters.
#' 
#' @param chromosome The chromosome you are restricting to (with or without a chr prefix).
#' @param qstart Query start coordinate of the range you are restricting to.
#' @param qend Query end coordinate of the range you are restricting to.
#' @param region Region formatted like chrX:1234-5678 instead of specifying chromosome, start and end separately.
#' @param streamlined Return Start_Position and Tumor_Smaple_Barcode as the only two MAF columns. Default is FALSE.
#' @param projection Obtain variants projected to this reference (one of grch37 or hg38).
#' @param seq_type The seq_type you want back, default is genome. Currently only genome is supported.
#' @param min_read_support Only returns variants with at least this many reads in t_alt_count (for cleaning up augmented MAFs).
#' @param verbose Set to FALSE to prevent ANY message to be printed. 
#' In most cases, this parameter should be left to TRUE. 
#' The parameter was added to accommodate for noisy output 
#' when running this function in a loop for retrieving SSM 
#' for multiple regions [GAMBLR.data::get_ssm_by_regions].
#' @param ... Any additional parameters.
#'
#' @return A data frame containing all mutations (MAF) in the specified region.
#'
#' @import dplyr stringr
#' 
#' @export
#'
#' @examples
#' my_mutations = get_ssm_by_region(region = "chr8:128,723,128-128,774,067")
#'
#' #specifying chromosome, start and end individually
#' my_mutations = get_ssm_by_region(chromosome = "8",
#'                                  qstart = 128723128,
#'                                  qend = 128774067)
#'
get_ssm_by_region = function(chromosome,
                             qstart,
                             qend,
                             region = "",
                             streamlined = FALSE,
                             projection = "grch37",
                             seq_type = "genome",
                             min_read_support = 3,
                             verbose = TRUE,
                             ...){
  
  #check seq type
  if(seq_type != "genome"){
    stop("Currently, SSM results are only available for genome samples (in GAMBLR.data). Please run this function with `this_seq_type` set to genome...")
  }
  
  if(verbose){
    #warn/notify the user what version of this function they are using
    message("Using the bundled SSM calls (.maf) calls in GAMBLR.data...")  
  }
  
  #check if any invalid parameters are provided
  check_excess_params(...)
  
  #get valid projections
  valid_projections = grep("meta", names(GAMBLR.data::sample_data), value = TRUE, invert = TRUE)
  
  #return SSMs based on the selected projection
  if(projection %in% valid_projections){
    this_maf = GAMBLR.data::sample_data[[projection]]$maf
  }else{
    stop(paste("please provide a valid projection. The following are available:",
               paste(valid_projections,collapse=", ")))
  }
  
  #drop poorly supported reads
  this_maf = dplyr::filter(this_maf, t_alt_count >= min_read_support)
  
  #split region into chunks (chr, start, end) and deal with chr prefixes based on the selected projection
  if(length(region) > 1){
    stop("You are providing more than one region, please refer to get_ssm_by_regions for multiple regions...")
  }
  
  if(!region == ""){
    region = gsub(",", "", region)
    split_chunks = unlist(strsplit(region, ":"))
    
    if(projection == "grch37"){
      region = stringr::str_replace(region, "chr", "")
    }
    
    chromosome = split_chunks[1]
    startend = unlist(strsplit(split_chunks[2], "-"))
    qstart = as.numeric(startend[1])
    qend = as.numeric(startend[2])
  }else{
    if(projection =="grch37"){
      chromosome = gsub("chr", "", chromosome)
    }
    region = paste0(chromosome, ":", qstart, "-", qend)
  }
  
  if(projection =="grch37"){
    chromosome = gsub("chr", "", chromosome)
  }
  
  #subset the maf to the specified region
  muts_region = dplyr::filter(this_maf, Chromosome == chromosome & Start_Position > qstart & Start_Position < qend)

  if(streamlined){
    muts_region = muts_region %>%
      dplyr::select(Start_Position, Tumor_Sample_Barcode)
  }
  
  return(muts_region)
}
