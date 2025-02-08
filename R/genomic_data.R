# functions for creating and working with S3 objects

#' Check and set the genome_build/projection
#'
#' This helper function checks the genome build of each genomic data object in
#' \code{genomic_data_list} (using \code{get_genome_build()}) and ensures
#' they are consistent. If all objects share a single, unique genome build,
#' that value is returned. If a user-specified genome build (\code{suggested})
#' is provided, it is compared to the inferred build and must match; otherwise,
#' an error is raised. If the genomic data objects have conflicting genome
#' builds or if no genome build can be inferred and no \code{suggested}
#' value is provided, the function stops with an error.
#'
#' @param genomic_data_list A list of genomic data objects. Each object should
#' have a genome build that can be retrieved by \code{get_genome_build()}.
#' @param suggested An optional character string specifying a genome build
#' (projection) to be used. If provided, it must match the genome build inferred
#' from the data objects.
#'
#' @return A character string representing the genome build to be used.
#' @export
#' @keywords internal
#' @examples
#' # Example 1: When genomic data objects all have the same genome build.
#' # Assuming maf_data and seg_data both have a genome build of "hg38":
#' genomic_data <- list(maf_data = maf_data, seg_data = seg_data)
#' projection <- check_get_projection(genomic_data, suggested = "hg38")
#'
#' # Example 2: When the genomic data objects conflict or no genome build
#' # is available.
#' # This will raise an error:
#' genomic_data <- list(maf_data = maf_data,
#'                      seg_data = seg_data_with_different_build)
#' projection <- check_get_projection(genomic_data, suggested = "hg38")
#' 
check_get_projection <- function(genomic_data_list, suggested) {
  # Extract genome builds from each genomic data object
  builds <- sapply(genomic_data_list, get_genome_build)
  uniq_builds <- unique(builds)
  
  if (length(uniq_builds) == 1) {
    # A single, consistent genome build was inferred.
    if (!missing(suggested) && suggested != uniq_builds) {
      stop("Mismatch between user-specified genome_build and
      the genome_build inferred from objects.")
    }
    return(uniq_builds)
  }
  
  if (length(uniq_builds) > 1) {
    # Conflicting genome builds among the objects.
    stop("Conflicting genome_build values found: ",
    paste(uniq_builds, collapse = ", "))
  }

  # No genome build could be inferred.
  if (missing(suggested)) {
    stop("No projection provided and genome_build
    cannot be inferred from the inputs.")
  }

  return(suggested)
}

## GAMBLR.data
#' Create Segmented Data
#'
#' This function creates segmented data from the given input.
#'
#' @param seg_df A data frame containing the segmented data.
#' @param genome_build Required character vector specifying the genome build
#' currently supported: "grch37" or "hg38".
#' @return A data frame with class attributes for segmented data.
#' @export
#' @examples
#' seg_df <- data.frame(...)
#' create_seg_data(seg_df, "grch37")
create_seg_data <- function(seg_df, genome_build) {
  if (!inherits(seg_df, "data.frame")) stop("data must be a data frame")
  if (!genome_build %in% c("grch37", "hg38")) stop("Invalid genome build")
  structure(seg_df,
            class = c("seg_data", class(seg_df)),
            genome_build = genome_build)
}


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

#' @export
#' @keywords internal
print.maf_data <- function(x, ...) {
  cat("MAF Data Object\n")
  cat("Genome Build:", attr(x, "genome_build"), "\n")
  cat("Showing first 10 rows:\n")
  # Convert to a plain data.frame (if not already) so that printing uses the default
  # data.frame print method rather than printing as a list.
  print(utils::head(as.data.frame(x), 10))
}


#' Get Genome Build
#'
#' This function retrieves the genome build attribute from the data.
#'
#' @param data A data frame with genome build attribute.
#' @return A string specifying the genome build.
#' @export
#' @keywords internal
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
#' @keywords internal
preserve_genomic_attributes <- function(new_data, old_data) {
  # Preserve the genome_build attribute
  attr(new_data, "genome_build") <- attr(old_data, "genome_build")
  
  # Combine the new data’s classes with the genomic classes
  new_data_classes <- class(new_data)
  # Ensure the genomic classes are at the front
  new_classes <- unique(c("maf_data", "genomic_data", new_data_classes))
  class(new_data) <- new_classes
  
  return(new_data)
}

#' Strip Genomic Data Classes
#'
#' This function removes custom classes associated with genomic data objects
#' (by default, "genomic_data", "maf_data", and "bed_data") from the class attribute
#' of an object. This can be useful when you want to revert an S3 object to its
#' underlying data.frame (or data.table) classes without converting the object.
#'
#' @param x An object, such as one of your genomic data objects.
#' @param classes A character vector of class names to remove. The default is
#'        c("genomic_data", "maf_data", "bed_data").
#' @return The object with the specified classes removed.
#' @export
#' @keywords internal
strip_genomic_classes <- function(x, classes = c("genomic_data", "maf_data", "bed_data")) {
  current_classes <- class(x)
  new_classes <- setdiff(current_classes, classes)
  class(x) <- new_classes
  return(x)
}


# S3 methods for genomic_data class
#' @export
#' @keywords internal
mutate.genomic_data <- function(.data, ...) {
  new_data <- dplyr::mutate(as.data.frame(.data), ...)
  preserve_genomic_attributes(new_data, .data)
}
#' @export
#' @keywords internal
filter.genomic_data <- function(.data, ...) {
  new_data <- dplyr::filter(as.data.frame(.data), ...)
  preserve_genomic_attributes(new_data, .data)
}
#' @export
#' @keywords internal
select.genomic_data <- function(.data, ...) {
  new_data <- dplyr::select(as.data.frame(.data), ...)
  preserve_genomic_attributes(new_data, .data)
}
#' @export
#' @keywords internal
rename.genomic_data <- function(.data, ...) {
  new_data <- dplyr::rename(as.data.frame(.data), ...)
  preserve_genomic_attributes(new_data, .data)
}
#' @export
#' @keywords internal
arrange.genomic_data <- function(.data, ...) {
  new_data <- dplyr::arrange(as.data.frame(.data), ...)
  preserve_genomic_attributes(new_data, .data)
}
#' @export
#' @keywords internal
group_by.genomic_data <- function(.data, ..., .add = FALSE) {
  new_data <- dplyr::group_by(as.data.frame(.data), ..., .add = .add)
  preserve_genomic_attributes(new_data, .data)
}
#' @export
ungroup.genomic_data <- function(x, ...) {
  new_data <- dplyr::ungroup(as.data.frame(x), ...)
  preserve_genomic_attributes(new_data, x)
}

#' Bind maf or other genomic data together
#'
#' @description Combine multiple maf_data objects and retain metadata such as genome_build.
#' This function will not allow you to combine maf_data objects that have different genome_build values.
#' An error will also be thrown if the same sample id is found in more than one of the inputs (if check_id is TRUE).
#'
#' @param ... All maf_data or seg_data objects to be combined.
#' @param check_id Logical. If TRUE (the default), the function will check for the presence of the expected ID column
#'        and for duplicate sample IDs across the inputs. Set to FALSE to skip this check.
#'
#' @return data.frame with combined data and preserved genome_build metadata.
#' @export
#'
#' @examples
#'
#' merged_maf = bind_genomic_data(maf1, maf2,check_id=FALSE)
#'
#' @keywords internal
bind_genomic_data <- function(..., check_id = TRUE) {

  in_list <- list(...)
  
  if ("maf_data" %in% class(in_list[[1]])) {
    # MAF format, ID column is Tumor_Sample_Barcode
    id_col <- "Tumor_Sample_Barcode"
  } else if ("seg_data" %in% class(in_list[[1]])) {
    # SEG format, ID column is ID
    id_col <- "ID"
  } else {
    stop(paste("Unsure how to merge:", class(in_list[[1]])))
  }
  
  # Ensure all inputs are either maf_data or seg_data objects
  if (!all(sapply(in_list, inherits, "maf_data")) &&
         !all(sapply(in_list, inherits, "seg_data"))) {
    stop("All inputs must be maf_data objects or seg_data objects.")
  }
  
  # Extract genome builds
  genome_builds <- unique(sapply(in_list, get_genome_build))
  
  if (length(genome_builds) > 1) {
    stop("Cannot bind seg_data or maf_data objects
    with different genome builds: ", 
         paste(genome_builds, collapse = ", "))
  }
  
  # If check_id is TRUE, verify that the expected ID column exists and
  # that IDs are unique.
  if (check_id) {
    # Collect unique sample IDs from each dataset
    id_sets <- lapply(in_list, function(df) {
      if (!(id_col %in% colnames(df))) {
        stop("ID column '", id_col, "' not found in input data.")
      }
      unique(df[[id_col]])
    })
    
    # Flatten the list and count occurrences of each ID
    all_ids <- unlist(id_sets)
    duplicate_ids <- names(table(all_ids)[table(all_ids) > 1])
    
    # If any ID is found in multiple datasets, throw an error
    if (length(duplicate_ids) > 0) {
      stop("Duplicate IDs found in multiple input data frames: ",
      paste(duplicate_ids, collapse = ", "))
    }
  }
  
  combined <- dplyr::bind_rows(in_list)
  attr(combined, "genome_build") <- genome_builds[1]  
# Assign the common genome build
  
  if (!"maf_data" %in% class(combined)) {
    class(combined) <- c("maf_data", "genomic_data", class(combined))
# Preserve class
  }
  
  return(combined)
}



#' Create BED Data
#'
#' This function creates BED (Browser Extensible Data) objects from the given input.
#' It assumes that the BED data should have columns corresponding to chromosome, start,
#' and end. If the second and third columns are not numeric (as expected for start and end),
#' the function will attempt to identify the proper columns by matching column names.
#'
#' In the output, the first three columns will be renamed to "chrom", "start", and "end".
#' If a fourth column exists, it is renamed to "name" (and any additional columns are preserved).
#'
#' Additionally, if a "name" column exists and its values are not unique, the function
#' will warn the user. The user can optionally supply a method to automatically fix the
#' names via the `fix_names` argument:
#'
#'   - If `fix_names = "chrom_start_end"`, the new name will be built as "chrom:start-end".
#'
#'   - If `fix_names = "concat"`, then the columns specified by `concat_cols` (using the
#'     original column names in the input data) will be concatenated to form the new name.
#'     By default, no separator is used, but a separator can be specified via the `sep`
#'     argument.
#'
#' After applying the fix, the function checks if the new names are unique. If they are not,
#' a warning is issued that includes up to five examples of duplicate names and the row numbers
#' where they occur.
#'
#' @param bed_df A data frame containing the BED data.
#' @param genome_build A string specifying the genome build
#' ("grch37" or "hg38").
#'        If NULL, the function will try to infer the genome build
#' from the object name.
#' @param fix_names Either NULL (the default), or one of "chrom_start_end"
#' or "concat".
#'        If not NULL and duplicate names are detected, the function will
#' apply the chosen fix.
#' @param concat_cols When `fix_names = "concat"`, a character vector
#' specifying which columns
#'        from the original data to merge.
#' @param sep The separator to use when concatenating columns if
#' fix_names = "concat".
#'        Defaults to "" (no separator).
#' @return A data frame with class attributes for BED data.
#' 
#' @export
#' 
#' @examples
#' 
#' # get a abed_data object for all aSHM regions
#' ashm_bed = create_bed_data(GAMBLR.data::grch37_ashm_regions,
#'                 fix_names = "concat",
#'                 concat_cols = c("gene","region"),
#'                 sep="-")
#' # the build is automatically inferred if it is in the variable name
#' get_genome_build(ashm_bed)
#' print(ashm_bed)
#' another_bed = create_bed_data(somatic_hypermutation_locations_GRCh37_v_latest,
#'                               fix_names = "concat",
#'                               concat_cols = c("chr_name","hg19_start","hg19_end"))
#' 
#' get_genome_build(another_bed)
#' 
#' # get a bed_data object for all gene regions and combine several columns to make a unique name
#' gene_regions <- create_bed_data(hg38_gene_coordinates,
#'                     fix_names = "concat",
#'                     sep="-",
#'                     concat_cols = c("chromosome","start","end","gene_name"))
#'                     
#' get_genome_build(gene_regions)
#'
create_bed_data <- function(bed_df,
                            genome_build = NULL,
                            fix_names = NULL,
                            concat_cols = NULL,
                            sep = "") {
  # Check that input is a data frame.
  if (!inherits(bed_df, "data.frame")) {
    stop("Input data must be a data frame")
  }
  
  # Capture the original data and column names (before any reordering or renaming)
  orig_df <- bed_df
  orig_names <- names(bed_df)
  
  # If genome_build is not provided, attempt to infer it from the object name.
  if (is.null(genome_build)) {
    object_name <- deparse(substitute(bed_df))
    possible_builds <- character(0)
    
    if (grepl("grch37", object_name, ignore.case = TRUE)) {
      possible_builds <- c(possible_builds, "grch37")
    }
    if (grepl("hg38", object_name, ignore.case = TRUE)) {
      possible_builds <- c(possible_builds, "hg38")
    }
    
    if (length(possible_builds) == 1) {
      genome_build <- possible_builds
    } else if (length(possible_builds) == 0) {
      stop("Could not determine genome build from object name; please supply genome_build argument.")
    } else {
      stop("Ambiguous genome build in object name; please supply genome_build argument explicitly.")
    }
  }
  
  # Validate genome build.
  if (!genome_build %in% c("grch37", "hg38")) {
    stop("Invalid genome build. Please choose either 'grch37' or 'hg38'.")
  }
  
  # Helper function to force column naming for the BED data.
  force_bed_column_names <- function(df) {
    new_names <- names(df)
    # Force first three columns to be "chrom", "start", "end"
    new_names[1:3] <- c("chrom", "start", "end")
    # If there's a fourth column, force it to "name"
    if (ncol(df) >= 4) {
      new_names[4] <- "name"
    }
    names(df) <- new_names
    return(df)
  }
  
  # Check if the first three columns (as supplied) are in the expected form.
  # We expect columns 2 and 3 (start and end) to be numeric.
  if (ncol(bed_df) >= 3 && is.numeric(bed_df[[2]]) && is.numeric(bed_df[[3]])) {
    # The data is assumed to be in the correct order.
    bed_df <- force_bed_column_names(bed_df)
  } else {
    # Attempt to guess the proper columns based on names.
    names_lower <- tolower(names(bed_df))
    
    chrom_idx <- which(names_lower %in% c("chrom", "chromosome"))
    start_idx <- which(names_lower %in% c("start", "start_position", "startpos"))
    end_idx   <- which(names_lower %in% c("end", "end_position", "endpos"))
    
    if (length(chrom_idx) != 1 || length(start_idx) != 1 || length(end_idx) != 1) {
      stop("Columns 2 and 3 (start and end) are not numeric and the chromosome/start/end columns ",
           "cannot be unambiguously identified from the column names.")
    }
    
    # Reorder the data frame so that the candidate columns come first.
    remaining_idx <- setdiff(seq_len(ncol(bed_df)), c(chrom_idx, start_idx, end_idx))
    new_order <- c(chrom_idx, start_idx, end_idx, remaining_idx)
    bed_df <- bed_df[, new_order, drop = FALSE]
    
    # After reordering, check that the new second and third columns are numeric.
    if (!is.numeric(bed_df[[2]]) || !is.numeric(bed_df[[3]])) {
      stop("After reordering based on column names, the start and end columns are not numeric.")
    }
    
    # Force the first three (and optionally the fourth) column names.
    bed_df <- force_bed_column_names(bed_df)
  }
  
  # If a "name" column exists, check that its values are unique.
  if (ncol(bed_df) >= 4) {
    if (anyDuplicated(bed_df[[4]]) > 0) {
      # If no fix is provided, issue a generic warning.
      if (is.null(fix_names)) {
        warning("The values in the 'name' column are not unique.")
      } else {
        # Apply the requested fix.
        if (fix_names == "chrom_start_end") {
          new_names_vec <- paste0(bed_df$chrom, ":", bed_df$start, "-", bed_df$end)
          bed_df[[4]] <- new_names_vec
          if (length(unique(new_names_vec)) != nrow(bed_df)) {
            # Identify duplicate examples.
            dup_idx <- which(duplicated(new_names_vec) | duplicated(new_names_vec, fromLast = TRUE))
            dup_names <- unique(new_names_vec[dup_idx])
            dup_info <- sapply(dup_names, function(nm) {
              rows <- which(new_names_vec == nm)
              paste0(nm, " (rows: ", paste(rows, collapse = ", "), ")")
            })
            warning("The 'chrom_start_end' fix did not result in a unique set of names. Examples: ",
                    paste(dup_info[1:min(5, length(dup_info))], collapse = "; "),
                    ". Please review your data or consider an alternative fix.")
          }
        } else if (fix_names == "concat") {
          if (is.null(concat_cols)) {
            stop("For fix_names = 'concat', you must supply concat_cols indicating which columns to merge.")
          }
          if (!is.character(concat_cols)) {
            stop("For fix_names = 'concat', concat_cols must be a character vector referring to the original column names.")
          }
          if (!all(concat_cols %in% orig_names)) {
            stop("One or more column names specified in concat_cols do not exist in the original data.")
          }
          # Build new names using the original data.
          # Use paste with the specified separator.
          new_names_vec <- do.call(paste, c(orig_df[, concat_cols, drop = FALSE], sep = sep))
          bed_df[[4]] <- new_names_vec
          if (length(unique(new_names_vec)) != nrow(bed_df)) {
            dup_idx <- which(duplicated(new_names_vec) | duplicated(new_names_vec, fromLast = TRUE))
            dup_names <- unique(new_names_vec[dup_idx])
            dup_info <- sapply(dup_names, function(nm) {
              rows <- which(new_names_vec == nm)
              paste0(nm, " (rows: ", paste(rows, collapse = ", "), ")")
            })
            warning("The 'concat' fix did not result in a unique set of names. Examples: ",
                    paste(dup_info[1:min(5, length(dup_info))], collapse = "; "),
                    ". Please review your data or consider an alternative fix.")
          }
        } else {
          stop("Invalid value for fix_names. Use 'chrom_start_end' or 'concat'.")
        }
      }
    }
  }
  # enforce strict matching of chr prefixing
  if(genome_build == "grch37"){
    if(any(grepl("chr",bed_df$chrom))){
      bed_df = mutate(bed_df,chrom = gsub("chr", "", chrom))
    }
  }
  # Create the S3 object with additional class attributes and genome_build attribute.
  structure(bed_df,
            class = c("bed_data", "genomic_data", class(bed_df)),
            genome_build = genome_build)
}

#' @export
print.bed_data <- function(x, ...) {
  cat("BED Data Object\n")
  cat("Genome Build:", attr(x, "genome_build"), "\n")
  cat("Showing first 10 rows:\n")
  # Convert to a plain data.frame (if not already) so that printing uses the default
  # data.frame print method rather than printing as a list.
  print(utils::head(as.data.frame(x), 10))
}


