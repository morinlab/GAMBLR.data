#' @title Assign CN to SSM.
#'
#' @description Annotate mutations with their copy number information.
#'
#' @details This function takes a sample ID with the `this_sample_id` parameter and annotates mutations with copy number information.
#' A variety of parameters are at hand for a customized workflow. For example, the user can specify if only coding mutations are of interest.
#' To do so, set `coding_only = TRUE`. This function internally calls `get_ssm_by_samples` and `get_sample_cn_segments`.
#' This function can also take a vector with genes of interest (`genes`) that the returned data frame will be restricted to.
#'
#' @param this_sample_id Sample ID of the sample you want to annotate.
#' @param genes Genes of interest.
#' @param this_seq_type Specified seq type for returned data.
#' @param projection specified genome projection that returned data is in reference to.
#' @param coding_only Optional. set to TRUE to restrict to only coding variants.
#' @param assume_diploid Optional. If no local seg file is provided, instead of defaulting to a GAMBL sample, this parameter annotates every mutation as copy neutral.
#' @param include_silent Logical parameter indicating whether to include silent mutations into coding mutations. Default is FALSE. This parameter only makes sense if coding only is set to TRUE.
#'
#' @return A list containing a data frame (MAF-like format) with two extra columns:
#' log.ratio is the log ratio from the seg file (NA when no overlap was found)
#' as well as the segmented copy number data with the same copy number information
#' CN is the rounded absolute copy number estimate of the region based on log.ratio (NA when no overlap was found)
#'
#' @rawNamespace import(data.table, except = c("last", "first", "between", "transpose"))
#' @import dplyr readr
#' @export
#'
#' @examples
#' cn_list = assign_cn_to_ssm(this_sample_id = "DOHH-2",
#'                            coding_only = TRUE)
#'
assign_cn_to_ssm = function(this_sample_id,
                            genes,
                            this_seq_type = "genome",
                            projection = "grch37",
                            coding_only = FALSE,
                            assume_diploid = FALSE,
                            include_silent = FALSE,
                            ...){
  
  #warn/notify the user what version of this function they are using
  message("Using the bundled CN segments (.seg) calls in GAMBLR.data...")
  
  #check if any invalid parameters are provided
  check_excess_params(...)
  
  #ensure only one sample ID is provided
  if(length(this_sample_id) > 1){
    stop("This function only supports queries of 1 sample ID at the time...")
  }

  #get maf
  maf_sample = get_ssm_by_sample(these_sample_ids = this_sample_id, 
                                 projection = projection, 
                                 this_seq_type = this_seq_type) %>%
    data.table::as.data.table()
  
  #maf filtering
  #silent mutations
  if(!include_silent){
    coding_class = coding_class[coding_class != "Silent"]
  }
  
  #coding mutations
  if(coding_only){
    maf_sample = dplyr::filter(maf_sample, Variant_Classification %in% coding_class)
  }
  
  #subset to genes of interest
  if(!missing(genes)){
    maf_sample = dplyr::filter(maf_sample, Hugo_Symbol %in% genes)
    if(nrow(maf_sample) == 0){
      stop("No variants left after filtering on the provided genes...")
    }
  }
  
  #get seg
  seg_sample = get_sample_cn_segments(these_sample_ids = this_sample_id, 
                                      projection = projection,  
                                      this_seq_type = this_seq_type)

  #annotate all CN segments as copy number neutral
  if(assume_diploid){
    diploid = dplyr::mutate(maf_sample, CN = 2)
    return(list(maf = diploid))
  }
  
  #wrangle the seg file
  seg_sample = seg_sample %>%
    dplyr::filter(end - start > 100) %>%
    mutate(chrom = gsub("chr", "", chrom)) %>%
    rename(Chromosome = chrom, Start_Position = start, End_Position = end) %>%
    mutate(across(LOH_flag, as.factor)) %>%
    data.table::as.data.table() %>%
    data.table::setkey(Chromosome, Start_Position, End_Position)
  
  #perform an overlap join and add CN columns from the seg file and subset MAF to basic columns (first 45)
  maf_tmp = data.table::foverlaps(maf_sample, seg_sample, type = "any")
  
  #rename and change order of columns to match expected format
  maf_with_segs = subset(maf_tmp, select = -c(ID, Start_Position, End_Position)) %>%
    rename(Start_Position = i.Start_Position, End_Position = i.End_Position) %>%
    dplyr::select(colnames(maf_sample), LOH_flag, log.ratio, CN)

  return(list(maf = maf_with_segs, seg = seg_sample))
}
