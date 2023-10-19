#' @title Process Regions objects.
#'
#' @description INTERNAL FUNCTION to harmonize genomic regions specified as character vectors or data frames.
#'
#' @details INTERNAL FUNCTION to harmonize genomic regions specified as character vectors or data frames.
#'
#' @param regions Character vector of genomic regions. If neither regions nor regions_df is specified, will use GAMBLR aSHM regions
#' @param regions_df Data frame of genomic regions with column names "chrom", "start", "end", "name"
#' @param region_padding Amount to pad the start and end coordinates by. [0]
#' @param skip_regions Character vector of genes to drop from GAMBLR aSHM regions.
#' @param only_regions Character vector of genes to include from GAMBLR aSHM regions.
#' @param projection Specify which genome build to use. [grch37]
#'
#' @return Numeric value.
#'
#' @noRd
#'
#' @examples
#' regions <- c("chr18:63074150-63369324")
#' process_regions(regions)
#' 
#' these_regions <- process_regions(
#'  only_regions = c("MYC", "BCL2", "BCL6")
#' )
#' regions_bed <- these_regions$regions_bed
#' regions_vec <- these_regions$regions
#' @export
process_regions <- function(regions_list = NULL,
                            regions_bed = NULL,
                            region_padding = 0,
                            skip_regions = NULL,
                            only_regions = NULL,
                            projection = "grch37") {
  # Use default ashm region table if no regions are provided
  if (is.null(regions_list)) {
    if (is.null(regions_bed)) {
      message("Using default GAMBLR aSHM regions. ")
      if (projection == "grch37") {
        regions_bed <- GAMBLR.data::grch37_ashm_regions %>%
          dplyr::mutate(chr_name = str_remove(chr_name, "chr")) %>% 
          dplyr::mutate(name = str_c(gene, region, sep = "_"))
      } else {
        regions_bed <- GAMBLR.data::hg38_ashm_regions %>%
          dplyr::mutate(name = str_c(gene, region, sep = "_"))
      }
      # Fix column names
      regions_bed <- regions_bed %>%
        rename_with(
          ~ str_remove(.x, "^hg.*_")
        ) %>%
        dplyr::rename(chrom = chr_name)
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
      df <- data.frame(
        chrom = chunks$chromosome,
        start = as.numeric(chunks$start),
        end = as.numeric(chunks$end)
      )
    }))
    if (!is.null(names(regions_list))) {
      regions_bed$name <- names(regions_list)
      regions_bed$gene <- names(regions_list)
    } else {
      message("WARNING: Regions provided as an unnamed character vector. It is strongly recommended to provide a named vector with a unique name per region, or a bed file-formatted data frame. ")
      regions_bed$name <- regions_bed$chrom
      regions_bed$gene <- regions_bed$chrom
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

