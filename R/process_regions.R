#' @title Process Regions objects.
#'
#' @description INTERNAL FUNCTION to harmonize genomic regions specified as character vectors or data frames.
#'
#' @details INTERNAL FUNCTION to harmonize genomic regions specified as character vectors or data frames.
#'
#' @param regions_list Character vector of genomic regions. If neither regions nor regions_df is specified, will use GAMBLR aSHM regions
#' @param regions_bed Data frame of genomic regions with column names "chrom", "start", "end", "name"
#' @param region_padding Amount to pad the start and end coordinates by. The default is 0 (no padding).
#' @param skip_regions Character vector of genes to drop from GAMBLR aSHM regions.
#' @param only_regions Character vector of genes to include from GAMBLR aSHM regions.
#' @param projection Specify which genome build projection to use. The default is "grch37", also accepts "hg38".
#' @param sort Set to TRUE to force regions_bed to be ordered on chromosome and coordinate
#'
#' @return A list with two objects, regions as a vector and in bed format.
#'
#' @export
#'
#' @examples
#' library(dplyr)
#'
#' regions <- setNames(
#'      c("chr1:10000-15000", "chr1:100000000-100005000"),
#'      c("one_region", "another_region")
#' )
#' process_regions(regions_list = regions)
#'
#' reg_bed = GAMBLR.data::grch37_ashm_regions %>%
#' dplyr::filter(chr_name == "chr17") %>%
#'   mutate(name = region, chrom = chr_name, start = hg19_start, end = hg19_end) %>%
#'   select(chrom, start, end, name)
#'
#' process_regions(regions_bed = reg_bed)
#'
process_regions <- function(regions_list = NULL,
                            regions_bed = NULL,
                            region_padding = 0,
                            skip_regions = NULL,
                            only_regions = NULL,
                            projection = "grch37",
                            sort = FALSE) {

  # Use default ashm region table if no regions are provided
  if (is.null(regions_list)) {
    if (is.null(regions_bed)) {
      message("Using default GAMBLR aSHM regions. ")
      if (projection == "grch37") {
        regions_bed <-  create_bed_data(grch37_ashm_regions,
                                        fix_names="concat",
                                        concat_cols=c("gene","region"),
                                        sep="_")
      } else if(projection=="hg38") {
        regions_bed <-  create_bed_data(hg38_ashm_regions,
                                        fix_names="concat",
                                        concat_cols=c("gene","region"),
                                        sep="_")
      }else{
        stop("unsupported projection!")
      }
      
      if (!is.null(skip_regions)) {
        # drop user-specified regions
        regions_bed <- regions_bed %>%
          dplyr::filter(!gene %in% skip_regions)
      }
      if (!is.null(only_regions)) {
        # keep only user-specified regions
        regions_bed <- regions_bed %>%
          dplyr::filter(gene %in% only_regions)
      }
    }

    required_cols <- c("chrom", "start", "end", "name")
    if (min(required_cols %in% colnames(regions_bed)) == 0) {
      stop("Provided regions_bed lacks required column names. Ensure columns chrom, start, end, and name are present. ")
    }

    # gene column is required for later joins
    if (!"gene" %in% colnames(regions_bed)) {
      regions_bed <- mutate(regions_bed, gene = name)
    }
  } else {
    # Convert character vector of regions to df
    regions_bed <- bind_rows(lapply(regions_list, function(x) {

      chunks <- region_to_chunks(x)
      if(projection=="grch37"){
        chunks$chromosome = gsub("chr","",chunks$chromosome)
      }else if(projection=="hg38" && !any(grepl("chr",chunks$chromosome))){
        chunks$chromosome = paste0("chr",chunks$chromosome)
      }
      df <- data.frame(
        chrom = chunks$chromosome,
        start = as.numeric(chunks$start),
        end = as.numeric(chunks$end)
      )
    }))
    if(sort){
      if(projection=="hg38"){
        chrom_order = c(paste0("chr",c(1:22)),"chrX","chrY")
      }else{
        chrom_order = c(c(1:22),"X","Y")
      }
      
      regions_bed = mutate(regions_bed,
                           chrom=factor(chrom,levels=chrom_order)) %>%
        arrange(chrom,start) %>%
        mutate(chrom = as.character(chrom))
    }
    if (!is.null(names(regions_list))) {
      regions_bed$name <- names(regions_list)
      regions_bed$gene <- names(regions_list)
    } else {
      regions_bed = mutate(regions_bed,name=paste0(chrom,":",start,"-",end))
    }
  }

  # Collapse regions with duplicate names
  if (length(unique(regions_bed$name)) < length(regions_bed$name)) {
    message("Warning: Multiple regions in the provided data frame have the same name. Merging these entries based on min(start) and max(end) per name value. ")
    regions_bed <- regions_bed %>%
      group_by(name) %>%
      mutate(
        start = min(start),
        end = max(end)
      ) %>%
      ungroup() %>%
      distinct()
  }

  regions_list <- unlist(apply(
    regions_bed,
    1,
    function(x) {
      # add specified padding around each region
      paste0(x[1], ":", as.numeric(x[2]) - region_padding, "-", as.numeric(x[3]) + region_padding)
    }
  ))
  names(regions_list) <- regions_bed$name

  return(
    list(
      regions_list = regions_list,
      regions_bed = regions_bed
    )
  )
}
