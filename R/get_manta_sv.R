#' @title Get Manta SVs
#'
#' @description Convenience function for retrieving Manta Structural Variants (SVs) from the bundled data [GAMBLR.data::sample_data].
#'
#' @details To obtain SV calls for multiple samples, give `these_sample_ids` a vector of sample IDs. 
#' Alternatively, the user can also provide the `these_samples_metadata` parameter to make use of an already subset metadata table. 
#' In this case, the returned SVs will be restricted to the sample_ids within that data frame. 
#' This function internally calls [GAMBLR.data::id_ease] to streamline sample ID/metadata parameters.
#' This function can also restrict the returned calls to any genomic regions specified within `chromosome`, `qstart`, `qend`,
#' or the complete region specified under `region` (in chr:start-end format), note that chromosome can be either prefixed or not prefixed.
#' Useful filtering parameters are also available, use `min_vaf` to set the minimum tumour VAF for a SV to be returned and `min_score`
#' to set the lowest Manta somatic score for a SV to be returned. `pair_status` can be used to return variants from either matched or unmatched samples.
#' In addition, the user can chose to return all variants, even the ones not passing the filter criteria. To do so, set `pass = FALSE` (default is TRUE).
#'
#' @param these_sample_ids Optional, a vector of multiple sample_id (or a single sample ID as a string) that you want results for.
#' @param these_samples_metadata Optional, a metadata table (with sample IDs in a column) to subset the return to. 
#' If not provided (and if `these_sample_ids` is not provided), the function will return all samples from the specified seq_type in the metadata.
#' @param projection The projection genome build. Default is grch37.
#' @param this_seq_type The this_seq_type you want back, default is genome.
#' @param chromosome Optional, the chromosome you are restricting to (can be prefixed or not prefixed).
#' @param qstart Optional, query start coordinate of the range you are restricting to.
#' @param qend Optional, query end coordinate of the range you are restricting to.
#' @param region Optional, region formatted like chrX:1234-5678 (chromosome can be prefixed or not prefixed) instead of specifying chromosome, start and end separately.
#' @param pairing_status Use to restrict results (if desired) to matched or unmatched results (default is to return all). This parameter takes the filtering condition as a string ("matched" or "unmatched").
#' @param min_vaf The minimum tumour VAF for a SV to be returned. Default is 0.1.
#' @param min_score The lowest Manta somatic score for a SV to be returned. Default is 40.
#' @param pass If TRUE (default) only return SVs that are annotated with PASS in the FILTER column. Set to FALSE to keep all variants, regardless if they PASS the filters.
#' @param verbose Set to FALSE to minimize the output to console. Default is TRUE. This parameter also dictates the verbosity of any helper function internally called inside the main function.
#' @param ... Any additional parameters.
#' 
#' @export
#' 
#' @import dplyr
#' 
#' @examples
#' #load packages
#' library(dplyr)
#' 
#' #lazily get every SV in the table with default quality filters
#' all_sv = get_manta_sv()
#'
#' #get all SVs DLBCL cell line samples
#' cell_line_meta = GAMBLR.data::sample_data$meta %>% 
#'   dplyr::filter(cohort == "DLBCL_cell_lines")
#'   
#' dlbcl_sv = get_manta_sv(these_samples_metadata = cell_line_meta)
#'
#' #get the SVs in a region around MYC
#' myc_locus_sv = get_manta_sv(region = "8:128723128-128774067")
#' 
get_manta_sv = function(these_sample_ids = NULL,
                        these_samples_metadata = NULL,
                        projection = "grch37",
                        this_seq_type = "genome",
                        chromosome,
                        qstart,
                        qend,
                        region,
                        pairing_status,
                        min_vaf = 0.1,
                        min_score = 40,
                        pass = TRUE,
                        verbose = FALSE,
                        ...){
  
  #warn/notify the user what version of this function they are using
  message("Using the bundled Manta SV (.bedpe) calls in GAMBLR.data...")
  
  #check if any invalid parameters are provided
  check_excess_params(...)
  
  #get valid projections
  valid_projections = grep("meta", names(GAMBLR.data::sample_data), value = TRUE, invert = TRUE)
  
  #get samples with the dedicated helper function
  metadata = id_ease(these_samples_metadata = these_samples_metadata,
                     these_sample_ids = these_sample_ids,
                     verbose = verbose,
                     this_seq_type = this_seq_type)
  
  sample_ids = metadata$sample_id
  
  #return manta SV based on the selected projection
  if(projection %in% valid_projections){
    manta_sv = GAMBLR.data::sample_data[[projection]]$bedpe %>% 
      dplyr::filter(tumour_sample_id %in% sample_ids)
  }else{
    stop(paste("please provide a valid projection. The following are available:",
               paste(valid_projections,collapse=", ")))
  }
  
  if(!missing(region)){
    region = gsub(",", "", region)
    split_chunks = unlist(strsplit(region, ":"))
    chromosome = split_chunks[1]
    startend = unlist(strsplit(split_chunks[2], "-"))
    qstart = startend[1]
    qend = startend[2]
  }
  
  manta_sv = manta_sv %>%
    dplyr::filter(VAF_tumour >= min_vaf,
                  SCORE >= min_score)
  
  if(verbose){
    no_manta = setdiff(metadata$sample_id, manta_sv$tumour_sample_id)
    
    if(length(no_manta) > 0){
      message(paste0("No Manta results found for ", length(no_manta), " samples..."))
      print(no_manta)
    }
  }
  
  #deal with chr prefixes based on the selected projection (if return is to be subset to regions...)
  if(!missing(region) || !missing(chromosome)){
    if(projection == "grch37"){
      if(grepl("chr", chromosome)){
        chromosome = gsub("chr", "", chromosome)
      }
    }else if(projection == "hg38"){
      if(!grepl("chr", chromosome)){
        chromosome = paste0("chr", chromosome)
      }
    }
    
    manta_sv = manta_sv %>%
      dplyr::filter((CHROM_A == chromosome & START_A >= qstart & START_A <= qend) | (CHROM_B == chromosome & START_B >= qstart & START_B <= qend))
  }
  
  if(verbose){
    message("\nThe following VCF filters are applied;")
    message(paste0("  Minimum VAF: ", min_vaf))
    message(paste0("  Minimum Score: ", min_score))
    message(paste0("  Only keep variants passing the quality filter: ", pass))
  }
  
  #PASS filter
  if(pass){
    manta_sv = manta_sv %>%
      dplyr::filter(FILTER == "PASS")
  }
  
  #pairing status filter
  if(!missing(pairing_status)){
    if(verbose){
      message(paste0("  Pairing status: ", pairing_status))
    }
    
    manta_sv = manta_sv %>%
      dplyr::filter(pair_status == pairing_status)
  }
  
  #convert to data frame and print some metrics
  manta_sv = as.data.frame(manta_sv)
  
  if(verbose){
    n_variants = nrow(manta_sv)
    unique_samples = unique(manta_sv$tumour_sample_id)
    message(paste0("\nReturning ", n_variants, " variants from ", length(unique_samples), " sample(s)"))
    message("\nDone!")
  }
  
  return(manta_sv)
}
