
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
#' @param these_sample_ids Optional, a vector of multiple sample_id (or a single
#'      sample ID as a string) that you want results for.
#' @param these_samples_metadata Optional, a metadata table (with sample IDs in
#'      a column) to subset the return to. If not provided (and if
#'      `these_sample_ids` is not provided), the function will return all
#'      samples from the specified seq_type in the metadata.
#' @param projection Reference genome build for the coordinates in the MAF file.
#'      The default is grch37.
#' @param this_seq_type The this_seq_type you want back, default is genome.
#' @param min_read_support Only returns variants with at least this many reads
#'      in t_alt_count.
#' @param include_silent Logical parameter indicating whether to include silent
#'      mutations into coding mutations. Default is TRUE.
#' @param verbose Set to FALSE to minimize the output to console. Default is
#'      TRUE. This parameter also dictates the verbosity of any helper function
#'      internally called inside the main function.
#' @param tool_name Optionally specify which tool to report variant from. The
#'      default is slms-3, also supports "publication" to return the exact
#'      variants as reported in the original papers.
#' @param ... Any additional parameters.
#'
#' @return data frame
#'
#' @import dplyr
#'
#' @export
#'
#' @examples
#'
#'  # Get mutations from exome data originally aligned to grch37
#' ssm_exomes_grch37 = get_coding_ssm(projection = "grch37",this_seq_type = "capture")
#' 
#' # Get mutations from genome data, hg38 build
#' ssm_genomes_hg38 = get_coding_ssm(projection = "hg38",this_seq_type = "genome")
#'
#' 
#'
#'
get_coding_ssm = function(
    these_sample_ids = NULL,
    these_samples_metadata = NULL,
    projection = "grch37",
    this_seq_type = "genome",
    tool_name = "slms-3",
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
        dplyr::filter(Tumor_Sample_Barcode %in% sample_ids) %>%
        dplyr::filter((tolower(!!sym("Pipeline")) == tool_name))
    
    if(!include_silent){
        coding_class = coding_class[coding_class != "Silent"]
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
    muts = create_maf_data(muts,projection)
    # use S3-safe version of dplyr function
    muts = mutate.genomic_data(muts,maf_seq_type = this_seq_type)
    return(muts)
}
