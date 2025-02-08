#' @title Get SSM By Samples.
#'
#' @description Get the SSMs (i.e. load MAF) for a single sample or a
#' collection of samples.
#'
#' @details Retrieve a maf for a specific sample or a set of samples.
#' Either specify the sample IDs of interest with `these_sample_ids`.
#' Or a metadata table subset to the sample IDs of interest with
#' `these_samples_metadata`.
#'
#' @param these_sample_ids A vector of one or more sample IDs that you
#' want results for.
#' @param these_samples_metadata Optional, a metadata table (with sample_id
#' column) to auto-subset the data to samples in that table before returning.
#' If not provided and these_sample_ids is also not provided, the function will
#' return SSM for all samples from the specified seq_type in the bundled
#' metadata.
#' @param this_seq_type Default is genome.
#' @param projection The projection genome build. Supports hg38 and grch37.
#' @param tool_name Optionally specify which tool to report variant from.
#' The default is slms-3, also supports "publication" to return the exact
#' variants as reported in the original papers.
#' @param verbose Enable for debugging/noisier output.
#' @param ... Any additional parameters.
#'
#' @return data frame in MAF format.
#'
#' @import dplyr
#'
#' @export
#'
#' @examples
#'
#' #Get genome-wide set of mutations from all DLBCL cell lines
#'
#' # 1. get our metadata for the DLBCL cell lines
#' cell_line_meta = get_gambl_metadata() %>%
#'   dplyr::filter(cohort == "DLBCL_cell_lines")
#'
#' # 2. get the SSMs for the DLBCL cell lines
#' dlbcl_maf = get_ssm_by_samples(these_samples_metadata = cell_line_meta)
#'
#' # 3. have a look:
#' dlbcl_maf %>% dplyr::group_by(Tumor_Sample_Barcode) %>%
#'               dplyr::count()
#'
get_ssm_by_samples <- function(these_sample_ids = NULL,
                               these_samples_metadata = NULL,
                               this_seq_type = "genome",
                               projection = "grch37",
                               tool_name = "slms-3",
                               verbose = FALSE,
                               ...) {

  #warn/notify the user what version of this function they are using
  message("Using the bundled SSM calls (.maf) calls in GAMBLR.data...")

  #check if any invalid parameters are provided
  check_excess_params(...)

  #get samples with the dedicated helper function
  metadata = id_ease(these_samples_metadata = these_samples_metadata,
                     these_sample_ids = these_sample_ids,
                     verbose = verbose,
                     this_seq_type = this_seq_type)

  sample_ids = metadata$sample_id

  #get valid projections
  valid_projections = grep("meta", names(GAMBLR.data::sample_data),
                           value = TRUE, invert = TRUE)

  #return SSMs based on the selected projection
  if(projection %in% valid_projections) {
    sample_ssm = GAMBLR.data::sample_data[[projection]]$maf %>%
      dplyr::filter(Tumor_Sample_Barcode %in% sample_ids) %>%
      dplyr::filter((tolower(!!sym("Pipeline")) == tool_name))
    sample_ssm <- bind_rows(sample_ssm,
      GAMBLR.data::sample_data[[projection]]$ashm %>%
        dplyr::filter(Tumor_Sample_Barcode %in% sample_ids) %>%
        dplyr::filter((tolower(!!sym("Pipeline")) == tool_name))
    )

  }else {
    stop(paste("please provide a valid projection. Available options:",
               paste(valid_projections,collapse=", ")))
  }


  # Handle possible duplicates
  sample_ssm <- sample_ssm %>%
    distinct(Tumor_Sample_Barcode,
             Chromosome,
             Start_Position,
             End_Position,
             .keep_all = TRUE)
  # bundle genome_build with the maf_data
  sample_ssm = create_maf_data(sample_ssm,projection)
  # use S3-safe version of dplyr function
  sample_ssm = mutate.genomic_data(sample_ssm,maf_seq_type = this_seq_type)
  return(sample_ssm)
}
