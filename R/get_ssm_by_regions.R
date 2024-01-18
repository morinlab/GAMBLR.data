#' @title Get SSM By Regions.
#'
#' @description Efficiently retrieve all mutations across a range of genomic regions.
#'
#' @details This function internally calls [GAMBLR.data::get_ssm_by_region] to retrieve SSM calls for the specified regions.
#'
#' @param these_sample_ids Optional, a vector of multiple sample_id (or a single sample ID as a string) that you want results for.
#' @param these_samples_metadata Optional, a metadata table (with sample IDs in a column) to subset the return to.
#' @param maf_data Optional data frame with mutations in MAF format.
#' If user provides a maf, the function trusts that the user has already subset this to samples of interest, correct seq_type.
#' i.e the following parameters are ignored; `these_samples_metadata`, `these_sample_ids`, and `this_seq_type`
#' If not provided (and if `these_sample_ids` is not provided), the function will return all samples from the specified seq_type in the metadata.
#' @param this_seq_type The this_seq_type you want back, default is genome.
#' @param mode Accepts 2 options of "slms-3" and "publication" to indicate which variant calls to use. Default is "slms-3".
#' @param regions_list A vector of regions in the chr:start-end format to restrict the returned SSM calls to.
#' @param regions_bed A data frame in BED format with the coordinates you want to retrieve (recommended).
#' This parameter can also accept an additional column with region names that will be added to the return if `use_name_column = TRUE`
#' @param streamlined If set to TRUE (default) only 3 columns will be kept in the returned data frame (start, sample_id and region_name).
#' @param use_name_column If your bed-format data frame has a name column (must be named "name") these can be used to name your regions in the returned data frame.
#' @param projection Obtain variants projected to this reference (one of grch37 or hg38), default is grch37.
#' @param min_read_support Only returns variants with at least this many reads in t_alt_count (for cleaning up augmented MAFs).
#' @param verbose Set to TRUE to maximize the output to console. Default is TRUE.
#' This parameter also dictates the verbosity of any helper function internally called inside the main function.
#' @param ... Any additional parameters.
#'
#' @return Returns a data frame of variants in MAF-like format.
#'
#' @import tibble dplyr tidyr
#'
#' @export
#'
#' @examples
#' #basic usage, adding custom names from bundled ashm data frame
#' regions_bed = dplyr::mutate(GAMBLR.data::grch37_ashm_regions,
#'   name = paste(gene, region, sep = "_"))
#'
#' ashm_basic_details = get_ssm_by_regions(regions_bed = regions_bed,
#'                                         streamlined = FALSE,
#'                                         use_name_column = TRUE)
#'
get_ssm_by_regions = function(these_sample_ids = NULL,
                              these_samples_metadata = NULL,
                              maf_data,
                              this_seq_type = "genome",
                              mode = "slms-3",
                              regions_list,
                              regions_bed,
                              streamlined = TRUE,
                              use_name_column = FALSE,
                              projection = "grch37",
                              min_read_support = 3,
                              verbose = FALSE,
                              ...){

  #check if any invalid parameters are provided
  check_excess_params(...)

  bed2region = function(x){
    paste0(x[1], ":", as.numeric(x[2]), "-", as.numeric(x[3]))
  }

  if(missing(regions_list)){
    if(!missing(regions_bed)){
      regions = apply(regions_bed, 1, bed2region)
    }else{
      warning("You must supply either regions_list or regions_bed")
    }
  }else{
    regions = regions_list
  }

  if(verbose){
    print(regions)
  }

  if(missing(maf_data)){
    #warn/notify the user what version of this function they are using
    message("Using the bundled SSM calls (.maf) calls in GAMBLR.data...")

    #get samples with the dedicated helper function
    metadata = id_ease(these_samples_metadata = these_samples_metadata,
                       these_sample_ids = these_sample_ids,
                       verbose = verbose,
                       this_seq_type = this_seq_type)

    region_mafs = lapply(regions, function(x){GAMBLR.data::get_ssm_by_region(region = x,
                                                                             these_samples_metadata = metadata,
                                                                             this_seq_type = this_seq_type,
                                                                             streamlined = streamlined,
                                                                             projection = projection,
                                                                             min_read_support = min_read_support,
                                                                             mode = mode,
                                                                             verbose = FALSE, #force to FALSE, suppressing noisy output
                                                                             ...)})
  }else{
    region_mafs = lapply(regions, function(x){GAMBLR.data::get_ssm_by_region(region = x,
                                                                             maf_data = maf_data,
                                                                             streamlined = streamlined,
                                                                             projection = projection,
                                                                             min_read_support = min_read_support,
                                                                             verbose = FALSE, #force to FALSE, suppressing noisy output
                                                                             ...)})

  }

  #deal with region names
  if(!use_name_column){
    rn = regions
  }else{
    if(missing(regions_bed)){
      stop("use_name_column = TRUE is only available if regions are provided with `regions_bed`")
    }else{
      rn = regions_bed[["name"]]
    }
  }

  #add the region name to the to-be-returned maf
  tibbled_data = tibble(region_mafs, region_name = rn)

  #unnest tibbled data
  unnested_df = tibbled_data %>%
    unnest_longer(region_mafs)

  if(streamlined){
    unlisted_df = mutate(unnested_df, start = region_mafs$Start_Position, sample_id = region_mafs$Tumor_Sample_Barcode) %>%
      dplyr::select(start, sample_id, region_name)
  }else{
      return(bind_rows(region_mafs))
  }
  return(unlisted_df)
}
