#' @title Get SSM By Samples.
#'
#' @description Get the SSMs (i.e. load MAF) for a single sample or a collection fo samples.
#'
#' @details Retrieve a maf for a specific sample or a set of samples.
#' Either specify the sample IDs of interest with `these_sample_ids`.
#' Or a metadata table subset to the sample IDs of interest with `these_samples_metadata`.
#'
#' @param this_sample_id A single sample ID you want the data from.
#' @param these_sample_ids A vector of sample IDs that you want results for.
#' @param these_samples_metadata Optional, a metadata table (with sample IDs in a column) to auto-subset the data to samples in that table before returning.
#' If not provided and these_sample_ids is also not provided, the function will return SSM for all samples from the specified seq_type in the bundled metadata.
#' @param this_seq_type Default is genome.
#' @param projection The projection genome build. Supports hg38 and grch37.
#' @param tool_name Optionally specify which tool to report variant from. The default is slms-3, also supports "publication" to return the exact variants as reported in the original papers.
#' @param these_genes A vector of genes to subset ssm to.
#' @param basic_columns Return first 43 columns of MAF rather than full details. Default is TRUE.
#' @param maf_cols if basic_columns is set to FALSE, the user can specify what columns to be returned within the MAF. This parameter can either be a vector of indexes (integer) or a vector of characters.
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
#' #load packages
#' library(dplyr)
#'
#' #return a MAF for DLBCL cell line
#' cell_line_meta = GAMBLR.data::sample_data$meta %>%
#'   dplyr::filter(cohort == "DLBCL_cell_lines")
#'
#' dlbcl_maf = get_ssm_by_samples(these_samples_metadata = cell_line_meta)
#'
get_ssm_by_samples <- function(these_sample_ids = NULL,
                               these_samples_metadata = NULL,
                               this_seq_type = "genome",
                               projection = "grch37",
                               tool_name = "slms-3",
                               this_study,
                               these_genes,
                               basic_columns = TRUE,
                               maf_cols = NULL,
                               verbose = FALSE,
                               ...){

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
  valid_projections = grep("meta", names(GAMBLR.data::sample_data), value = TRUE, invert = TRUE)

  #return SSMs based on the selected projection
  if(projection %in% valid_projections){
    sample_ssm = GAMBLR.data::sample_data[[projection]]$maf %>%
        dplyr::filter(Tumor_Sample_Barcode %in% sample_ids) %>%
        dplyr::filter((tolower(!!sym("Pipeline")) == tool_name))
    sample_ssm <- bind_rows(
        sample_ssm,
        GAMBLR.data::sample_data[[projection]]$ashm %>%
            dplyr::filter(Tumor_Sample_Barcode %in% sample_ids) %>%
            dplyr::filter((tolower(!!sym("Pipeline")) == tool_name))
    )
  }else{
    stop(paste("please provide a valid projection. The following are available:",
               paste(valid_projections,collapse=", ")))
  }

  if(!missing(these_genes)){
    sample_ssm = sample_ssm %>%
      dplyr::filter(Hugo_Symbol %in% these_genes)
  }

  # Optionally return variants from a particular study
  if(!missing(this_study)){
    sample_ssm <- sample_ssm %>%
      dplyr::filter((tolower(!!sym("Study")) == this_study))
  }  

  #subset maf to only include first 43 columns (default)
  if(basic_columns){
    sample_ssm = dplyr::select(sample_ssm, c(1:45))
  }

  #subset maf to a specific set of columns (defined in maf_cols)
  if(!is.null(maf_cols) && !basic_columns){
    sample_ssm = dplyr::select(sample_ssm, all_of(maf_cols))
  }

  # Handle possible duplicates
  sample_ssm <- sample_ssm %>%
    distinct(Tumor_Sample_Barcode, Chromosome, Start_Position, End_Position, .keep_all = TRUE)

  return(sample_ssm)
}


#' @rdname get_ssm_by_samples
#'
#' @examples
#' #basic usage, using a single sample ID
#' dohh2_maf = get_ssm_by_sample(this_sample_id = "DOHH-2")
#'
#' @export
#'
get_ssm_by_sample = function(this_sample_id = NULL,
                             these_samples_metadata = NULL,
                             this_seq_type = "genome",
                             projection = "grch37",
                             these_genes,
                             basic_columns = TRUE,
                             maf_cols = NULL,
                             verbose = FALSE,
                             ...){

  get_ssm_by_samples(these_sample_ids = this_sample_id,
                     these_samples_metadata = these_samples_metadata,
                     this_seq_type = this_seq_type,
                     projection = projection,
                     these_genes = these_genes,
                     basic_columns = basic_columns,
                     maf_cols = maf_cols,
                     verbose = verbose,
                     ...)
}
