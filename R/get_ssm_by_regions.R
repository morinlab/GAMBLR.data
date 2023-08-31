#' @title Get SSM By Regions.
#'
#' @description Efficiently retrieve all mutations across a range of genomic regions.
#'
#' @details This function internally calls [GAMBLR.data::get_ssm_by_region] to retrieve SSM calls for the specified regions.
#'
#' @param regions_vector Either provide a vector of regions in the chr:start-end format. 
#' If this parameter is provided, the function will auto set `use_name_column` to FALSE.
#' @param regions_bed Provide a data frame in BED format with the coordinates you want to retrieve (recommended). 
#' This parameter can also accept an additional column with region names that will be added to the return if `use_name_column = TRUE` 
#' @param streamlined If set to TRUE (default) only 3 columns will be kept in the returned data frame (start, sample_id and region_name). 
#' Set to FALSE to return 6 columns. See `basic_columns` for more information on what columns to include in the return.
#' @param basic_columns Set to TRUE for returning a maf with standarad 45 columns. Default is FALSE. If TRUE, `streamlined` and `use_name_column` will be ignored.
#' @param use_name_column If your bed-format data frame has a name column (must be named "name") these can be used to name your regions. 
#' This parameter. is nott accepted if `regions_vector` is provided.
#' @param projection Obtain variants projected to this reference (one of grch37 or hg38).
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
get_ssm_by_regions = function(regions_vector,
                              regions_bed,
                              streamlined = TRUE,
                              basic_columns = FALSE,
                              use_name_column = FALSE,
                              projection = "grch37",
                              min_read_support = 3,
                              verbose = FALSE,
                              ...){
  
  #warn/notify the user what version of this function they are using
  message("Using the bundled SSM calls (.maf) calls in GAMBLR.data...")
  
  #check if any invalid parameters are provided
  check_excess_params(...)
  
  #streamline parameter combinations
  if(basic_columns){
    streamlined = FALSE
  }
  
  bed2region = function(x){
    paste0(x[1], ":", as.numeric(x[2]), "-", as.numeric(x[3]))
  }
  
  if(missing(regions_vector)){
    if(!missing(regions_bed)){
      regions = apply(regions_bed, 1, bed2region)
    }else{
      warning("You must supply either regions_vector or regions_bed")
    }
  }else{
    use_name_column = FALSE #override behavior of this parameter if a region vector is provided.
    regions = regions_vector
  }

  if(verbose){
    print(regions)  
  }
  
  region_mafs = lapply(regions, function(x){GAMBLR.data::get_ssm_by_region(region = x,
                                                                           streamlined = streamlined,
                                                                           projection = projection, 
                                                                           min_read_support = min_read_support,
                                                                           verbose = FALSE, #force to FALSE, suppressing noisy output
                                                                           ...)})

  #return the standard 45 columns (default)
  if(basic_columns){
    print("bind_rows")
    return(bind_rows(region_mafs))
  }

  #deal with region names
  if(!use_name_column){
    rn = regions
  }else{
    rn = regions_bed[["name"]]
  }
  
  #add the region name to the to-be-returned maf
  tibbled_data = tibble(region_mafs, region_name = rn)
  
  #unnest tibbled data
  unnested_df = tibbled_data %>%
    unnest_longer(region_mafs)
  
  #return only 3 columns, if streamlined or 6 columns if FALSE (and basic_columns = FALSE)
  if(streamlined){
    unlisted_df = mutate(unnested_df, start = region_mafs$Start_Position, sample_id = region_mafs$Tumor_Sample_Barcode) %>%
      dplyr::select(start, sample_id, region_name)
  }else{
    unlisted_df = mutate(unnested_df, Chromosome = region_mafs$Chromosome, End_Position = region_mafs$End_Position, Start_Position = region_mafs$Start_Position, Tumor_Sample_Barcode = region_mafs$Tumor_Sample_Barcode) %>%
      dplyr::select(Chromosome, Start_Position, End_Position, Tumor_Sample_Barcode, region_name)
  }
  return(unlisted_df)
}
