#' @title Get Coding SSM Status.
#'
#' @description Tabulate mutation status (SSM) for a set of genes.
#'
#' @details This function takes a vector of gene symbols and subsets the
#' incoming MAF to specified genes. If no genes are provided, the function will
#' default to all lymphoma genes.  Currently only the grch37 genome build is
#' supported for hotspot annotation and review for this version of the function.
#'
#' @param gene_symbols A vector of gene symbols for which the mutation status
#'      will be tabulated. If not provided, lymphoma genes will be returned
#'      by default.
#' @param these_samples_metadata The metadata for samples of interest to be
#'      included in the returned matrix. Only the column "sample_id" is
#'      required. If not provided, the example metadata is used as default.
#' @param maf_data data frame in maf format. Must be in the grch37 projection.
#' @param include_hotspots Logical parameter indicating whether hotspots object
#'      should also be tabulated. Default is TRUE.
#' @param keep_multihit_hotspot Logical parameter indicating whether to keep the
#'      gene annotation as mutated when the gene has both hot spot and
#'      non-hotspot mutation. Default is FALSE. If set to TRUE, will report the
#'      number of non-hotspot mutations instead of tabulating for just mutation
#'      presence.
#' @param review_hotspots Logical parameter indicating whether hotspots object
#'      should be reviewed to include functionally relevant mutations or rare
#'      lymphoma-related genes. Default is TRUE.
#' @param genes_of_interest A vector of genes for hotspot review. Currently only
#'      FOXO1, MYD88, and CREBBP are supported.
#' @param genome_build Reference genome build for the coordinates in the MAF
#'      file. The default is hg19 genome build.
#'
#' @return A data frame with tabulated mutation status.
#'
#' @import dplyr tidyr
#' @export
#'
#' @examples
#' coding_tabulated_df = get_coding_ssm_status(
#'  maf_data = GAMBLR.data::sample_data$grch37$maf,
#'  gene_symbols = "EGFR"
#' )
#'
#' #all lymphoma genes from bundled NHL gene list
#' coding_tabulated_df = get_coding_ssm_status()
#'
get_coding_ssm_status = function(
        gene_symbols,
        these_samples_metadata,
        maf_data,
        include_hotspots = TRUE,
        keep_multihit_hotspot = FALSE,
        review_hotspots = TRUE,
        genes_of_interest = c("FOXO1", "MYD88", "CREBBP"),
        genome_build = "hg19",
        ...
    ){

    # check if any invalid parameters are provided
    check_excess_params(...)

    # check the projection
    if(!genome_build == "hg19"){
        stop(
            "Currently only grch37 projection (hg19 genome build) is supported."
        )
    }

    if(missing(gene_symbols)){
        message(
            "No gene_symbols provided, defaulting to all lymphoma genes."
        )
        gene_symbols <- GAMBLR.data::lymphoma_genes$Gene
    }

    if(missing(these_samples_metadata)){
        these_samples_metadata <- get_gambl_metadata()
    }

    coding_ssm <- maf_data %>%
        dplyr::filter(
            Variant_Classification %in% c(
                "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del",
                "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation",
                "Nonstop_Mutation", "Silent", "Splice_Region", "Splice_Site",
                "Targeted_Region", "Translation_Start_Site"
            )
        )

    coding <- coding_ssm %>%
        dplyr::filter(
            Hugo_Symbol %in% gene_symbols
        ) %>%
        dplyr::select(Tumor_Sample_Barcode, Hugo_Symbol) %>%
        dplyr::rename(
            "sample_id" = "Tumor_Sample_Barcode",
            "gene" = "Hugo_Symbol"
        ) %>%
        unique() %>%
        dplyr::mutate(mutated = 1)

    samples_table <- dplyr::select(
        these_samples_metadata,
        sample_id
    )
    wide_coding <- pivot_wider(
        coding,
        names_from = "gene",
        values_from = "mutated",
        values_fill = 0
    )
    all_tabulated <- left_join(
        samples_table,
        wide_coding
    )
    all_tabulated <- all_tabulated %>%
        replace(is.na(.), 0)

    # include hotspots if user chooses to do so
    if(include_hotspots){
        # first annotate
        annotated <- GAMBLR.data::annotate_hotspots(
            coding_ssm
        )

        # review for the supported genes
        if(review_hotspots){
            annotated = GAMBLR.helpers::review_hotspots(
                annotated,
                genes_of_interest = genes_of_interest,
                genome_build = genome_build
            )
        }

        message("annotating hotspots")

        hotspots <- annotated %>%
            dplyr::filter(Hugo_Symbol %in% genes_of_interest) %>%
            dplyr::select(Tumor_Sample_Barcode, Hugo_Symbol, hot_spot) %>%
            dplyr::rename(
                "sample_id" = "Tumor_Sample_Barcode",
                "gene" = "Hugo_Symbol"
            ) %>%
            dplyr::mutate(gene = paste0(gene, "HOTSPOT")) %>%
            unique() %>%
            dplyr::mutate(mutated = ifelse(hot_spot == "TRUE", 1, 0)) %>%
            dplyr::filter(mutated == 1) %>%
            dplyr::select(-hot_spot)

        # long to wide hotspots, samples are tabulated with 0 if no hotspot is detected
        wide_hotspots <- pivot_wider(
            hotspots,
            names_from = "gene",
            values_from = "mutated",
            values_fill = 0
        )
        # join with the ssm object
        all_tabulated <- left_join(
            all_tabulated,
            wide_hotspots
        )
        all_tabulated <- all_tabulated %>%
            replace(is.na(.), 0)

        all_tabulated <- all_tabulated %>%
            dplyr::select(where(~ any(. != 0)))

        all_tabulated <- as.data.frame(all_tabulated)
        # make SSM and hotspots non-redundant by giving priority
        # to hotspot feature and setting SSM to 0
        for (hotspot_site in colnames(wide_hotspots)[grepl("HOTSPOT", colnames(wide_hotspots))]){
            message(hotspot_site)
            this_gene <- gsub("HOTSPOT", "", hotspot_site)
            redundant_features <- all_tabulated %>%
                dplyr::select(starts_with(this_gene))

            # if not both the gene and the hotspot are present, go to
            # the next iteration
            if(ncol(redundant_features)!= 2) next
            message("OK")
            # if both gene and it's hotspot are in the matrix, give priority to hotspot feature
            all_tabulated[(all_tabulated[, this_gene] >0 & all_tabulated[, paste0(this_gene, "HOTSPOT")] == 1),][,c(this_gene, paste0(this_gene, "HOTSPOT"))][, this_gene] = 0

            # in case gene has both hotspot and another mutation in the same gene,
            # keep both tabulated as multihits
            if(keep_multihit_hotspot){
                # determine which samples have hot spot and another mutation in same gene
                multihits <- annotated %>%
                    dplyr::filter(Hugo_Symbol == this_gene) %>%
                    group_by(Tumor_Sample_Barcode) %>%
                    dplyr::mutate(n_mut = n()) %>%
                    dplyr::filter(
                        n_mut > 1
                    ) %>%
                    dplyr::distinct(Tumor_Sample_Barcode, n_mut, hot_spot) %>%
                    # account for cases with both hotspot and not hotspot to avoid
                    # double-counting the number of mutations
                    mutate_at(vars(hot_spot), ~replace_na(., "FALSE")) %>%
                    dplyr::mutate(
                        n_mut = ifelse(
                            hot_spot == "TRUE",
                            n_mut - 1,
                            n_mut
                        )
                    ) %>%
                    group_by(Tumor_Sample_Barcode) %>%
                    dplyr::arrange(n_mut) %>%
                    slice_head() %>%
                    ungroup %>%
                    select(-hot_spot)

                # Return the annotation of this gene to mutated in these samples
                all_tabulated <- all_tabulated %>%
                    left_join(
                        .,
                        multihits,
                        by = c("sample_id" = "Tumor_Sample_Barcode")
                    ) %>%
                    dplyr::mutate(
                        {{this_gene}} := ifelse(
                                !is.na(n_mut),
                                n_mut,
                                !!!syms(this_gene)
                            )
                    ) %>%
                    select(- n_mut)
            }

        }

    }
    return(all_tabulated)

}