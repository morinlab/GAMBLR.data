
#' @title Get genes from one or more gene lists
#'
#' @description Retrieve gene names from bundled lymphoma gene lists.
#'
#' @details Complete lists of genes described as significantly mutated in large
#' lymphoma studies have been curated and provided with this package.
#'
#' @param entities Optional vector specifying one or more lymphoma entities
#'      e.g. MCL, DLBCL, BL.
#' @param curated_only Specify FALSE to retrieve all genes or leave default
#'      for the curated subset.
#' @param gene_format Specify what to return as output. Can be one of:
#'      * "symbol" (the default): list of gene symbols
#'      * "ensembl": list of ENSEMBLE IDs
#'      * "data.frame": data frame with column Gene and per-entity gene status
#' @param version Specify which version to return. Currently supported versions
#'      are 0.0 (legacy version from original GAMBLR), 0.1, and _latest. The
#'      latter will always point to the highest numeric version of the genes.
#'
#' @return A character vector of gene symbol or Ensembl IDs or a data frame.
#'
#' @import dplyr
#' @export
#'
#' @examples
#' all_dlbcl_genes <- get_genes(entities = "DLBCL", curated_only = FALSE)
#' all_curated_genes <- get_genes()

get_genes <- function(
        entities = c("DLBCL", "MCL", "BL"),
        curated_only = TRUE,
        gene_format = "symbol",
        version = "_latest"
    ) {

    # We need to handle the legacy version (0.0) separately
    # because it was not pathology-specific
    if (version == "0.0") {
        message(
            paste(
            "You have requested the legacy version of the lymphoma genes.",
            "It is only provided here for backwards compatibility.",
            "We recommend you use the most recent version by keeping the",
            "default value of the argument version. Use at your own risk."
            )
        )

        legacy_lymphoma_genes <- eval(
                parse(
                    text = paste0(
                        "GAMBLR.data::",
                        "lymphoma_genes_lymphoma_genes_v0.0")
                    )
            )

        legacy_lymphoma_genes <- legacy_lymphoma_genes %>%
            dplyr::select(
                Gene, ensembl_gene_id,
                intersect(
                    colnames(.),
                    entities
                )
            )
        if (curated_only) {
            # drop any row where all pathologies have FALSE
            legacy_lymphoma_genes <- legacy_lymphoma_genes %>%
                dplyr::filter(
                    ! if_all(3:ncol(.), ~ . == "FALSE")
                )
        }

        if (gene_format == "symbol") {
            return(legacy_lymphoma_genes$Gene)
        } else if (gene_format == "ensembl") {
            return(legacy_lymphoma_genes$ensembl_gene_id)
        } else if (gene_format == "data.frame") {
            return(legacy_lymphoma_genes)
        } else {
            stop(
                "You requested output format that is not supported."
            )
        }
    }

    #construct file name using entity
    entities <- tolower(entities)

    r_objects <- entities %>%
        paste0(
            "lymphoma_genes_",
            .,
            "_v",
            version
        )

    # check for unsupported gene sets
    all_files <- system.file(
        "extdata",
        package = "GAMBLR.data"
    ) %>%
    list.files(
        recursive = TRUE,
        full.names = TRUE
    )

    all_files <- gsub(
        ".*extdata/",
        "",
        all_files
    )

    all_files <- all_files[grepl("lymphoma_genes", all_files)]

    all_files <- all_files[grepl(version, all_files)]

    available_entities <- gsub("(.*/\\s*(.*$))", "\\2", all_files)

    available_entities <- gsub(".tsv", "", available_entities)

    missing_sets <- setdiff(
        entities,
        available_entities
    )

    if (length(missing_sets) > 0) {
        warning(
            paste(
                "The gene set for the entity",
                missing_sets,
                "is not available and will not be returned."
            )
        )
        r_objects <- r_objects[
                grepl(
                    paste(
                        available_entities,
                        collapse="|"
                    ),
                    r_objects
                )
            ]
        entities <- entities[entities %in% available_entities]
    }

    # Combine all lists into one df
    # Do it in a way that when GAMBLR.data is ot imported, the data objects
    # are still available
    all_entities_data <- list()

    for (i in seq_along(r_objects)) {
        all_entities_data[[i]] <- eval(
            parse(
                text = paste0(
                    "GAMBLR.data::",
                    r_objects[i]
                )
            )
        )
    }

    all_entities_data <- all_entities_data %>%
        # only select necessary columns
        lapply(
            .,
            `[`,
            ,
            c("ensembl_gene_id", "Gene", "curated")
        )

    names(all_entities_data) <- toupper(entities)

    all_entities_data <- all_entities_data %>%
        bind_rows(.id = "entity")

    if (curated_only) {
        # drop any row where curated is FALSE
        all_entities_data <- all_entities_data %>%
            dplyr::filter(
                curated == "TRUE"
            )
    }

    if (gene_format == "symbol") {
        return(all_entities_data$Gene %>% unique %>% sort)
    } else if (gene_format == "ensembl") {
        return(all_entities_data$ensembl_gene_id %>% unique %>% sort)
    } else if (gene_format == "data.frame") {
        return(
            all_entities_data %>%
            select(-curated) %>%
            mutate(is_gene = "TRUE") %>%
            tidyr::pivot_wider(
                names_from = "entity",
                values_from = "is_gene"
            ) %>%
            replace(is.na(.), "FALSE")
        )
    } else {
        stop(
            "You requested output format that is not supported."
        )
    }
}

#' @title Produce colour palettes from your metadata.
#'
#' @description Given a data frame with at least one column, the function will
#' determine whether a colour palette exists and assign the colours to all
#' levels of data in that column.
#'
#' @details This helper function seeks to help you standardize colour mappings
#' within and across projects. It will return either a vector or a list for
#' compatability with ggplot and ComplexHeatmap, respectively.
#'
#' @param this_df Provide a data frame with at least one column. Required.
#' @param check Optionally, whether to perform checks for unsupported values
#'      and return helpful errors on exit (rather than happily returning an
#'      incomplete palette).
#' @param as_list Set to TRUE if you want a named list separating the colours
#'      by the original column names, otherwise all mappings will be in a
#'      single named vector.
#'
#' @return Either a vector or list of Hex codes.
#'
#' @import dplyr
#' @export
#'
#' @examples
#' \dontrun{
#' dplyr::select(
#'  GAMBLR::get_gambl_metadata(),
#'  pathology,
#'  COO_consensus,
#'  EBV_status_inf) %>%
#' get_mapped_colours()
#' }
#'

get_mapped_colours <- function(
        this_df,
        check = FALSE,
        as_list = FALSE
    ) {

    column_names <- colnames(this_df)

    # try to map every column to the colour palette using the name and,
    # if available, user-specified aliases
    mapped_list <- list()
    mapped_vector <- c()
    for (col_name in column_names) {

        unique_values <- unique(this_df[[col_name]])

        df <- dplyr::filter(
            GAMBLR.data::colour_codes,
            name %in% unique_values
        ) %>%
        dplyr::select(
            name,
            colour
        ) %>%
        unique()

        if (any(!unique_values %in% df$name)) {
            message(
                "missing one or more of the values in this set of colours:"
            )
            missing <- unique_values[which(!unique_values %in% df$name)]
            message(
                paste(
                    missing,
                    collapse = ", "
                )
            )

            if (check) {
                stop(
                    paste(
                    "you should correct this issue by modifying, dropping",
                    "or setting the offending values to NA, where applicable"
                    )
                )
            }
        }

        col_vec <- df$colour
        names(col_vec) <- df$name
        mapped_list[[col_name]] <- col_vec
        mapped_vector <- c(
            mapped_vector,
            col_vec
        )
    }

    if (as_list) {
        return(mapped_list)
    } else {
        return(mapped_vector)
    }

}

#' @title Get standardized colours for lymphoid cancers.
#'
#' @description Retrieve and visualize standardized colours
#'      schemes for lymphoid cancers.
#'
#' @details Colours hand picked to represent various common entities and
#'      clinical variables relevant for lymphoid cancers.
#'
#' @param show_available Set to TRUE to see what options are available.
#' @param this_category Optionally supply one of the available categories to
#'      see a subset of the options that are available.
#' @param this_group Optionally supply one of the available groups to see the
#'      palette just for this group.
#' @param as_named_vector Whether to return the colors as named vector.
#' @param drop_alias When FALSE, shows the redundant colours with their aliases.
#' @param legacy_mode When TRUE, will return named list similar to the first implementation of get_gambl_colours
#'
#' @return A data frame or named character vector of colour Hex codes.
#'
#' @import dplyr ggplot2 tibble
#' @export
#'
#' @examples
#' \dontrun{
#' get_colours(show_available = TRUE)
#' # printout shows that "subgroup" is one option to narrow it down,
#' # supply this to the function as this_category:
#' get_colours(show_available = TRUE, this_category = "subgroup")
#' # printout and plot shows several options.
#' # Pick the one you want to visualize it in isolation
#' get_colours(show_available = TRUE, this_group = "LymphGen")
#' # if satisfied, get the result for use with ggplot
#' col_vec <- get_colours(this_group = "LymphGen", as_named_vector = TRUE)
#' ggplot(...) + scale_fill_manual(values = col_vec)
#' }
#'

get_colours <- function(
        show_available = FALSE,
        this_category,
        this_group,
        as_named_vector = FALSE,
        drop_alias = TRUE,
        legacy_mode = FALSE
    ) {

    if (legacy_mode) {
        allcols <- GAMBLR.data::colour_codes$colour
        names(allcols) <- GAMBLR.data::colour_codes$name

        return(allcols)
    }

    if (drop_alias) {
        GAMBLR.data::colour_codes <- dplyr::filter(
            GAMBLR.data::colour_codes,
            is.na(is_alias)
        )
    }

    if (show_available) {
        if (missing(this_category) & missing(this_group)) {

            message(
                paste(
                    "Supply a category using this_category parameter.",
                    "Current options for 'category' are:"
                )
            )

            p <- GAMBLR.data::colour_codes %>%
                group_by(category) %>%
                mutate(
                    n = n()
                ) %>%
                slice_head() %>%
                rename(
                    c("example" = "name")
                )

            print(p)

            message(
                paste(
                    "Supply a group using this_group parameter.",
                    "Current options for 'group' are:"
                )
            )

            p <- GAMBLR.data::colour_codes  %>%
                group_by(group) %>%
                mutate(
                    n = n()
                ) %>%
                slice_head() %>%
                rename(
                    c("example"="name")
                )

            print(p)

        } else if (!missing(this_group)) {
            this_group_df <- GAMBLR.data::colour_codes %>%
                dplyr::filter(
                    group == this_group
                )

            allcols <- this_group_df$colour
            names(allcols) <- this_group_df$name

            p <- this_group_df %>%
                ggplot(
                    aes(
                        x = name,
                        fill = name,
                        y = 1
                    )
                ) +
                geom_col() +
                theme(
                    legend.position = "none"
                ) +
                facet_wrap(
                    ~category,
                    scales = "free_y"
                ) +
                scale_fill_manual(
                    values = allcols
                ) +
                coord_flip()

            print(p)

            return()
        } else if (!missing(this_category)) {
            this_category_df <- GAMBLR.data::colour_codes %>%
                dplyr::filter(
                    category == this_category
                )

            allcols <- this_category_df$colour
            names(allcols) <- this_category_df$name

            p <- this_category_df %>%
                group_by(group) %>%
                mutate(
                    n =n ()
                ) %>%
                slice_head()

            message("Current options for 'group' within this category are:")

            g <- pull(p ,group)
            g <- paste(g, collapse = ",")
            message(g)

            p <- GAMBLR.data::colour_codes %>%
                dplyr::filter(
                    category == this_category
                ) %>%
                ggplot(
                    aes(
                        x = name,
                        fill = name,
                        y = 1
                    )
                ) +
                geom_col() +
                theme(
                    legend.position = "none"
                ) +
                facet_wrap(
                    ~group,
                    scales = "free_y"
                ) +
                scale_fill_manual(
                    values = allcols
                ) +
                coord_flip()

            print(p)
            return()
        }

    }

    if (missing(this_category) & missing(this_group)) {
        stop(
            paste(
                "Provide a category or group via this_category or this_group.",
                "To see what's available run this function with",
                "show_available = TRUE"
            )
        )
    } else if (!missing(this_group)) {
        colour_list <- dplyr::filter(
            GAMBLR.data::colour_codes,
            group == this_group
        ) %>%
        dplyr::select(
            name,
            colour
        ) %>%
        column_to_rownames("name")
    } else if (!missing(this_category)) {
        colour_list <- dplyr::filter(
            GAMBLR.data::colour_codes,
            category == this_category
        ) %>%
        dplyr::select(name, colour) %>%
        column_to_rownames("name")
    }
    if (as_named_vector) {
        #useful for ggplot scale_X_manual
        allcols <- colour_list$colour
        names(allcols) <- rownames(colour_list)
        return(allcols)
    }
    return(colour_list)
}
