#' Create MAF Data
#'
#' This function creates MAF (Mutation Annotation Format) data from the given input.
#'
#' @param maf_df A data frame containing the MAF data.
#' @param genome_build A string specifying the genome build ("grch37" or "hg38").
#' @return A data frame with class attributes for MAF data.
#' @export
create_maf_data <- function(maf_df, genome_build) {
  if (!inherits(maf_df, "data.frame")) stop("data must be a data frame")
  if (!genome_build %in% c("grch37", "hg38")) stop("Invalid genome build")
  
  structure(maf_df,
            class = c("maf_data", "genomic_data", class(maf_df)),  #  "genomic_data" for generic methods
            genome_build = genome_build)
}

#' Get Genome Build
#'
#' This function retrieves the genome build attribute from the data.
#'
#' @param data A data frame with genome build attribute.
#' @return A string specifying the genome build.
#' @export
get_genome_build <- function(data) {
  attr(data, "genome_build")
}

#' Preserve Genomic Attributes
#'
#' This function preserves the genomic attributes and class after dplyr operations.
#'
#' @param new_data A data frame resulting from dplyr operations.
#' @param old_data The original data frame with genomic attributes.
#' @return A data frame with preserved genomic attributes.
#' @export
preserve_genomic_attributes <- function(new_data, old_data) {
  attr(new_data, "genome_build") <- attr(old_data, "genome_build")
  class(new_data) <- class(old_data)
  return(new_data)
}

# S3 methods for genomic_data class
#' @export
mutate.genomic_data <- function(.data, ...) {
  new_data <- dplyr::mutate(as.data.frame(.data), ...)
  preserve_genomic_attributes(new_data, .data)
}
#' @export
filter.genomic_data <- function(.data, ...) {
  new_data <- dplyr::filter(as.data.frame(.data), ...)
  preserve_genomic_attributes(new_data, .data)
}
#' @export
select.genomic_data <- function(.data, ...) {
  new_data <- dplyr::select(as.data.frame(.data), ...)
  preserve_genomic_attributes(new_data, .data)
}
#' @export
rename.genomic_data <- function(.data, ...) {
  new_data <- dplyr::rename(as.data.frame(.data), ...)
  preserve_genomic_attributes(new_data, .data)
}
#' @export
arrange.genomic_data <- function(.data, ...) {
  new_data <- dplyr::arrange(as.data.frame(.data), ...)
  preserve_genomic_attributes(new_data, .data)
}
#' @export
group_by.genomic_data <- function(.data, ..., .add = FALSE) {
  new_data <- dplyr::group_by(as.data.frame(.data), ..., .add = .add)
  preserve_genomic_attributes(new_data, .data)
}
#' @export
ungroup.genomic_data <- function(x, ...) {
  new_data <- dplyr::ungroup(as.data.frame(x), ...)
  preserve_genomic_attributes(new_data, x)
}

#' Bind maf data together
#'
#' @description Combine multiple maf_data objects and retain metadata such as genome_build. This function
#' will not allow you to combine maf_data objects that have different genome_build values. An error will also
#' be thrown if the same sample id is found in more than one of the inputs
#' @param ... All maf_data or seg_data objects to be combined
#'
#' @return data.frame
#' @export
#'
#' @examples
#'
#' merged_maf = bind_genomic_data(maf1,maf2)
#'
#'
bind_genomic_data <- function(...) {
  
  in_list <- list(...)
  if("maf_data" %in% class(in_list[[1]])){
    #MAF format, ID column is Tumor_Sample_Barcode
    id_col = "Tumor_Sample_Barcode"
  }else if("seg_data" %in% class(in_list[[1]])){
    #SEG format, ID column is ID
    id_col = "ID"
  }else{
    stop(paste("unsure how to merge:",class(in_list[[1]])))
  }
  # Ensure all inputs are seg_data or maf_data objects
  if (!all(sapply(in_list, inherits, "maf_data")) & !all(sapply(in_list, inherits, "seg_data"))) {
    stop("All inputs must be maf_data objects or seg_data objects.")
  }
  
  # Extract genome builds
  genome_builds <- unique(sapply(in_list, get_genome_build))
  
  if (length(genome_builds) > 1) {
    stop("Cannot bind seg_data or maf_data objects with different genome builds: ", paste(genome_builds, collapse = ", "))
  }
  
  # Collect unique sample IDs from each dataset
  id_sets <- lapply(in_list, function(df) {
    if (!(id_col %in% colnames(df))) {
      stop("ID column '", id_col, "' not found in input data.")
    }
    unique(df[[id_col]])  # Get unique IDs from each data frame
  })
  
  # Flatten the list and count occurrences of each ID
  all_ids <- unlist(id_sets)
  duplicate_ids <- names(table(all_ids)[table(all_ids) > 1])
  
  # If any ID is found in multiple datasets, throw an error
  if (length(duplicate_ids) > 0) {
    stop("Duplicate IDs found in multiple input data frames: ", paste(duplicate_ids, collapse = ", "))
  }
  
  combined <- bind_rows(in_list)
  attr(combined, "genome_build") <- genome_builds[1]  # Assign the common genome build
  if(!"maf_data" %in% class(combined)){
    class(combined) <- c("maf_data","genomic_data", class(combined))  # Preserve class
  }
  return(combined)
}

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
