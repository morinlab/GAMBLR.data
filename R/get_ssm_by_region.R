#' @title Get SSM By Region.
#'
#' @description Retrieve all SSMs from the GAMBL database within a single genomic coordinate range.
#'
#' @details This function lets the user specify a region of interest for returning SSM calls within that region.
#' There are multiple ways a region can be specified. For example, the user can provide the full region in a "region" format (chr:start-end) to the `region` parameter.
#' Or, the user can provide chromosome, start and end coordinates individually with `chr`, `start`, and `end` parameters.
#'
#' @param these_sample_ids Optional, a vector of multiple sample_id (or a single sample ID as a string) that you want results for.
#' @param these_samples_metadata Optional, a metadata table (with sample IDs in a column) to subset the return to.
#' If not provided (and if `these_sample_ids` is not provided), the function will return all samples from the specified seq_type in the metadata.
#' @param maf_data Optional data frame with mutations in MAF format.
#' If user provides a maf, the function trusts that the user has already subset this to samples of interest, correct seq_type.
#' i.e the following parameters are ignored; `these_samples_metadata`, `these_sample_ids`, and `this_seq_type`
#' @param chromosome The chromosome you are restricting to (with or without a chr prefix).
#' @param qstart Query start coordinate of the range you are restricting to.
#' @param qend Query end coordinate of the range you are restricting to.
#' @param region Region formatted like chrX:1234-5678 instead of specifying chromosome, start and end separately.
#' @param streamlined Return Start_Position and Tumor_Smaple_Barcode as the only two MAF columns. Default is FALSE.
#' @param projection Obtain variants projected to this reference (one of grch37 or hg38).
#' @param this_seq_type The seq_type you want back, default is genome.
#' @param tool_name Optionally specify which tool to report variant from. The default is slms-3, also supports "publication" to return the exact variants as reported in the original papers.
#' @param this_study Optionally specify first name of the author for the paper
#'      from which the variants should be returned for.
#' @param verbose Set to FALSE to prevent ANY message to be printed.
#' In most cases, this parameter should be left to TRUE.
#' The parameter was added to accommodate for noisy output
#' when running this function in a loop for retrieving SSM
#' for multiple regions [GAMBLR.data::get_ssm_by_regions].
#' @param ... Any additional parameters.
#'
#' @return A data frame containing all mutations (MAF) in the specified region.
#'
#' @import dplyr
#'
#' @examples
#' my_mutations = get_ssm_by_region(region = "chr8:128,723,128-128,774,067")
#'
#' #specifying chromosome, start and end individually
#' my_mutations = get_ssm_by_region(chromosome = "8",
#'                                  qstart = 128723128,
#'                                  qend = 128774067)
#'
get_ssm_by_region = function(these_sample_ids = NULL,
                             these_samples_metadata = NULL,
                             maf_data,
                             chromosome,
                             qstart,
                             qend,
                             region = "",
                             streamlined = FALSE,
                             projection = "grch37",
                             this_seq_type = "genome",
                             tool_name = "slms-3",
                             this_study,
                             verbose = FALSE,
                             ...){

  if(verbose){
    if(missing(maf_data)){
      #warn/notify the user what version of this function they are using
      message("Using the bundled SSM calls (.maf) calls in GAMBLR.data...")
    }
  }

  #check if any invalid parameters are provided
  check_excess_params(...)

  #get samples with the dedicated helper function
  metadata = id_ease(these_samples_metadata = these_samples_metadata,
                     these_sample_ids = these_sample_ids,
                     verbose = verbose,
                     this_seq_type = this_seq_type)

  sample_ids = metadata$sample_id

  

  # Optionally return variants from a particular study
  if(!missing(this_study)){
    this_maf <- this_maf %>%
      dplyr::filter((!!sym("Study")) == this_study)
  }

  #split region into chunks (chr, start, end) and deal with chr prefixes based on the selected projection
  if(length(region) > 1){
    stop("You are providing more than one region, please refer to get_ssm_by_regions for multiple regions...")
  }

  if(!region == ""){
    region = gsub(",", "", region)
    split_chunks = unlist(strsplit(region, ":"))

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

  if(projection == "grch37"){
    chromosome = gsub("chr", "", chromosome)
  }

  #return SSMs based on the selected projection
  if(missing(maf_data)){
    # Filter by position on-the-fly to avoid wastefully building the same large MAF each time
    this_maf = GAMBLR.data::sample_data[[projection]]$maf %>%
      dplyr::filter(Chromosome == chromosome & Start_Position > qstart & Start_Position < qend) %>%
      dplyr::filter(Tumor_Sample_Barcode %in% sample_ids) %>%
      dplyr::filter((tolower(!!sym("Pipeline")) == tool_name))
    muts_region <- GAMBLR.data::sample_data[[projection]]$ashm %>%
      dplyr::filter(Chromosome == chromosome & Start_Position > qstart & Start_Position < qend) %>%
      dplyr::filter(Tumor_Sample_Barcode %in% sample_ids) %>%
      dplyr::filter((tolower(!!sym("Pipeline")) == tool_name)) %>%
      bind_rows(this_maf, .)
  }else{
    muts_region = dplyr::filter(maf_data, Tumor_Sample_Barcode %in% sample_ids) %>%
      dplyr::filter(Chromosome == chromosome & Start_Position > qstart & Start_Position < qend)
  }
  
  # Handle possible duplicates
  muts_region <- muts_region %>%
    distinct(Tumor_Sample_Barcode, Chromosome, Start_Position, End_Position, .keep_all = TRUE)

  if(streamlined){
    muts_region = muts_region %>%
      dplyr::select(Start_Position, Tumor_Sample_Barcode)
  }
  muts_region = create_maf_data(muts_region,projection)
  # use S3-safe version of dplyr function
  muts_region = mutate.genomic_data(muts_region,maf_seq_type = this_seq_type)
  return(muts_region)
}
