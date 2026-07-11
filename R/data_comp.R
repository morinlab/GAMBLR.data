
#' @title Get genes from one or more gene lists
#'
#' @description Retrieve gene names from bundled lymphoma gene lists.
#'
#' @details Gene lists are read from the bundled normalized reference database
#' (`inst/extdata/gambl_reference.db`), which is built from the LLMPP curated
#' lists. Confidence tier and data version are columns in that database, so the
#' previous per-entity, per-version `.rda` objects are no longer required.
#'
#' @param entities Optional vector specifying one or more lymphoma entities
#'      e.g. MCL, DLBCL, BL.
#' @param curated_only When TRUE (default) return the curated set (Tier 1 and
#'      Tier 2) and drop retired Tier 3 genes. Set FALSE to include all tiers.
#'      Ignored when `tier` is supplied.
#' @param tier Optional integer vector of confidence tiers to keep, e.g.
#'      `tier = 1` for high-confidence genes only, or `tier = c(1, 2)`. When
#'      supplied it takes precedence over `curated_only`.
#' @param gene_format Specify what to return as output. Can be one of:
#'      * "symbol" (the default): list of gene symbols
#'      * "ensembl": list of Ensembl IDs
#'      * "data.frame": one row per gene with a `<ENTITY>_Tier` column per
#'        requested entity (the gene's tier there, or NA if absent). Filter
#'        these columns to focus on tiers, e.g.
#'        `dplyr::filter(df, DLBCL_Tier == 1, FL_Tier == 1) |> dplyr::pull(Gene)`.
#' @param version Data version to return. Defaults to "_latest", the current
#'      LLMPP-sourced version in the reference database. Legacy numeric versions
#'      (0.0, 0.1) are no longer bundled.
#'
#' @return A character vector of gene symbols or Ensembl IDs, or a data frame.
#'
#' @import dplyr tidyr
#'
#' @export
#'
#' @examples
#' # high-confidence (Tier 1) DLBCL genes
#' dlbcl_tier1 <- get_genes(entities = "DLBCL", tier = 1)
#'
#' # genes that are Tier 1 in DLBCL, FL and BL simultaneously
#' gene_df <- get_genes(entities = c("DLBCL", "FL", "BL"), gene_format = "data.frame")
#' shared_tier1 <- dplyr::filter(gene_df, DLBCL_Tier == 1, FL_Tier == 1, BL_Tier == 1)
#'
#' all_curated_genes <- get_genes()

get_genes <- function(
        entities = c("DLBCL", "MCL", "BL"),
        curated_only = TRUE,
        tier = NULL,
        gene_format = "symbol",
        version = "_latest"
    ) {

    # Gene lists are sourced from the normalized reference database (built from
    # the LLMPP curated lists) rather than per-entity, per-version .rda objects.
    # "version" and "tier" are columns, so lymphoma_genes_*_v0.1/v0.2/v_latest
    # are no longer needed.
    if (!version %in% c("_latest", "latest")) {
        message(
            paste(
                "Legacy per-version gene lists (e.g. 0.0, 0.1) are no longer",
                "bundled. Data is now sourced from the LLMPP curated lists;",
                "returning the current version."
            )
        )
    }

    # `tier`, when supplied, selects exactly those confidence tiers (e.g.
    # tier = 1 for high-confidence genes only) and takes precedence over
    # `curated_only`. Otherwise curated_only == TRUE keeps the curated set
    # (Tier 1 + Tier 2) and drops retired Tier 3 genes; FALSE returns all tiers.
    keep_tiers <- if (!is.null(tier)) {
        as.integer(tier)
    } else if (curated_only) {
        c(1L, 2L)
    } else {
        c(1L, 2L, 3L)
    }

    con <- gambl_reference_db()
    entities_u <- toupper(entities)

    available <- dplyr::tbl(con, "gene_entity") %>%
        dplyr::distinct(entity) %>%
        dplyr::pull(entity)

    missing_sets <- setdiff(entities_u, available)
    if (length(missing_sets) > 0) {
        warning(
            paste(
                "The gene set for the entity",
                paste(missing_sets, collapse = ", "),
                "is not available and will not be returned."
            )
        )
        entities_u <- intersect(entities_u, available)
    }

    dat <- dplyr::tbl(con, "gene_entity") %>%
        dplyr::filter(entity %in% entities_u, tier %in% keep_tiers) %>%
        dplyr::select(entity, ensembl_gene_id, Gene = gene, tier) %>%
        dplyr::collect()

    if (gene_format == "symbol") {
        return(dat$Gene %>% unique %>% sort)
    } else if (gene_format == "ensembl") {
        return(dat$ensembl_gene_id %>% unique %>% sort)
    } else if (gene_format == "data.frame") {
        # one row per gene; one <ENTITY>_Tier column per requested entity
        # holding the gene's tier there (NA if absent), so callers can do
        # filter(df, DLBCL_Tier == 1, FL_Tier == 1, ...) %>% pull(Gene)
        return(
            dat %>%
                dplyr::distinct(Gene, ensembl_gene_id, entity, tier) %>%
                tidyr::pivot_wider(
                    id_cols = c(Gene, ensembl_gene_id),
                    names_from = "entity",
                    values_from = "tier",
                    names_glue = "{entity}_Tier"
                ) %>%
                dplyr::select(Gene, ensembl_gene_id,
                              dplyr::any_of(paste0(entities_u, "_Tier"))) %>%
                dplyr::arrange(Gene)
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
#' @import dplyr tidyr
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
#' @import dplyr ggplot2 tidyr
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
        colour_codes <- dplyr::filter(
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
