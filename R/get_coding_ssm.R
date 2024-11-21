#' @title Get Coding SSMs
#'
#' @description Convenience function for loading coding Simple Somatic Mutations
#'      (SSM) from the bundled data [GAMBLR.data::sample_data].
#'
#' @details This "bare bones" function was developed to retrieve coding SSM
#'      calls for non-GSC-users. Effectively retrieve coding SSM calls. Multiple
#'      filtering parameters are available for this function. For more
#'      information on how to implement the filtering parameters, refer to the
#'      parameter descriptions as well as examples in the vignettes. This
#'      function depends on the bundled sample data in this package.
#'
#' @param limit_cohort Optional, supply this to restrict mutations to one or
#'      more cohorts in a vector.
#' @param exclude_cohort Optional, supply this to exclude mutations from one or
#'      more cohorts in a vector.
#' @param limit_pathology Optional, supply this to restrict mutations to one
#'      pathology.
#' @param limit_samples Optional, supply this to restrict mutations to a vector
#'      of sample_id (instead of sub-setting using the provided metadata).
#' @param these_sample_ids Optional, a vector of multiple sample_id (or a single
#'      sample ID as a string) that you want results for.
#' @param these_samples_metadata Optional, a metadata table (with sample IDs in
#'      a column) to subset the return to. If not provided (and if
#'      `these_sample_ids` is not provided), the function will return all
#'      samples from the specified seq_type in the metadata.
#' @param force_unmatched_samples Optional argument for forcing unmatched
#'      samples, using [GAMBLR.data::get_ssm_by_samples].
#' @param projection Reference genome build for the coordinates in the MAF file.
#'      The default is grch37.
#' @param this_seq_type The this_seq_type you want back, default is genome.
#' @param basic_columns Set to FALSE to override the default behavior of
#'      returning only the first 45 columns of MAF data.
#' @param maf_cols if `basic_columns` is set to FALSE, the user can specify what
#'      columns to be returned within the MAF. This parameter can either be a
#'      vector of indexes (integer) or a vector of characters (matching columns
#'      in MAF).
#' @param min_read_support Only returns variants with at least this many reads
#'      in t_alt_count.
#' @param include_silent Logical parameter indicating whether to include silent
#'      mutations into coding mutations. Default is TRUE.
#' @param verbose Set to FALSE to minimize the output to console. Default is
#'      TRUE. This parameter also dictates the verbosity of any helper function
#'      internally called inside the main function.
#' @param ... Any additional parameters.
#'
#' @return data frame
#'
#' @import dplyr
#'
#' @export
#'
#' @examples
#' #load pacakges
#' library(dplyr)
#'
#' #return SSMs in reference to GRCh37:
#' ssm_grch37 = get_coding_ssm(this_seq_type = "genome")
#'
#' #return SSMs in reference to hg38:
#' cell_line_meta = GAMBLR.data::sample_data$meta %>%
#'   dplyr::filter(cohort == "DLBCL_cell_lines")
#'
#' ssm_hg38 = get_coding_ssm(
#'      projection = "hg38",
#'      these_samples_metadata = cell_line_meta,
#'      this_seq_type = "genome"
#' )
#'
get_coding_ssm = function(
    limit_cohort,
    exclude_cohort,
    limit_pathology,
    limit_samples,
    these_sample_ids = NULL,
    these_samples_metadata = NULL,
    force_unmatched_samples,
    projection = "grch37",
    this_seq_type = "genome",
    basic_columns = TRUE,
    maf_cols = NULL,
    min_read_support = 3,
    include_silent = TRUE,
    verbose = FALSE,
    ...
){

    # Warn/notify the user what version of this function they are using
    message("Using the bundled SSM calls (.maf) calls in GAMBLR.data...")

    #check if any invalid parameters are provided
    check_excess_params(...)

    # Get valid projections
    valid_projections = grep(
        "meta",
        names(GAMBLR.data::sample_data),
        value = TRUE,
        invert = TRUE
    )

    #get samples with the dedicated helper function
    metadata = id_ease(
        these_samples_metadata = these_samples_metadata,
        these_sample_ids = these_sample_ids,
        verbose = verbose,
        this_seq_type = this_seq_type
    )

    sample_ids = metadata$sample_id


    if(!projection %in% valid_projections){
        stop(
            paste(
                "Provide a valid projection. The following are available:",
                paste(
                    valid_projections,
                    collapse = ", "
                )
            )
        )
    }

    #return SSMs based on the selected projection
    muts = GAMBLR.data::sample_data[[projection]]$maf %>% 
        dplyr::filter(Tumor_Sample_Barcode %in% sample_ids)

    if(!include_silent){
        coding_class = coding_class[coding_class != "Silent"]
    }

    #limit cohort
    if(!missing(limit_cohort)){
        metadata = metadata %>%
            dplyr::filter(cohort %in% limit_cohort)
    }

    #exclude cohort
    if(!missing(exclude_cohort)){
        metadata = metadata %>%
            dplyr::filter(!cohort %in% exclude_cohort)
    }

    #limit pathology
    if(!missing(limit_pathology)){
        metadata = metadata %>%
            dplyr::filter(pathology %in% limit_pathology)
    }

    #limit samples
    if(!missing(limit_samples)){
        metadata = metadata %>%
            dplyr::filter(sample_id %in% limit_samples)

        #check if any samples were not found in the metadata
        not_in_meta = setdiff(limit_samples, metadata$sample_id)

        #check so that the sample ID actually exists
        if(length(not_in_meta > 0)){
            message(
                paste0(
                    "The follwoing sample(s) were not found in the metadata: ",
                    not_in_meta
                )
            )
        }
    }

    sample_ids = pull(metadata, sample_id)

    # Drop variants with low read support (default is 3),
    # enforce sample IDs and keep only coding variants
    muts = dplyr::filter(muts, t_alt_count >= min_read_support) %>%
        dplyr::filter(Tumor_Sample_Barcode %in% sample_ids) %>%
        dplyr::filter(Variant_Classification %in% coding_class)

    # Filter maf on selected sample ids
    muts = muts %>%
        dplyr::filter(Tumor_Sample_Barcode %in% sample_ids)

    mutated_samples = length(unique(muts$Tumor_Sample_Barcode))
    message(
        paste(
            "after linking with metadata, we have mutations from",
            mutated_samples,
            "samples"
        )
    )

    # Drop rows for these samples so we can swap in the force_unmatched outputs
    if(!missing(force_unmatched_samples)){
        muts = muts %>%
            dplyr::filter(!Tumor_Sample_Barcode %in% force_unmatched_samples)

        nsamp = length(force_unmatched_samples)
        message(
            paste(
                "Dropping variants from",
                nsamp,
                "samples and replacing with force_unmatched outputs"
            )
        )

        # Get replacements using get_ssm_by_samples
        fu_muts = GAMBLR.data::get_ssm_by_samples(
            these_sample_ids = force_unmatched_samples
        )
        muts = bind_rows(muts, fu_muts)
    }

    # Subset maf to a specific set of columns (defined in maf_cols)
    if(!is.null(maf_cols) && !basic_columns){
        muts = dplyr::select(muts, all_of(maf_cols))
    }

    return(muts)
}
