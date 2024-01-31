#' @title Get ASHM Count Matrix.
#'
#' @description Prepare a matrix with one row per sample and one column per
#' region using a set of hypermutated regions.
#'
#' @details Values are the number of mutations in that patient in the region.
#'
#' @param regions_bed A bed file with one row for each region.
#' @param maf_data Optionally provide a data frame in the MAF format, otherwise
#'      the database will be used.
#' @param these_samples_metadata This is used to complete your matrix. All GAMBL
#'      samples will be used by default. Provide a data frame with at least
#'      sample_id for all samples if you are using non-GAMBL data.
#' @param this_seq_type The seq type to return results for. Only used if no
#'      metadata is provided with these_samples_metadata.
#'
#' @return matrix
#'
#' @import dplyr tibble
#' @export
#'
#' @examples
#' regions_bed <- dplyr::mutate(
#'      GAMBLR.data::grch37_ashm_regions,
#'      name = paste(gene, region, sep = "_")
#' )
#'
#' matrix <- get_ashm_count_matrix(
#'      regions_bed = regions_bed,
#'      this_seq_type = "genome"
#' )
#'
get_ashm_count_matrix = function(
        regions_bed,
        maf_data,
        these_samples_metadata,
        this_seq_type
    ){
    if(missing(this_seq_type)){
        if(missing(these_samples_metadata)){
            stop(
                "Please supply either the this_seq_type or a metadata from which it can be retrieved"
            )
        }
        this_seq_type <- these_samples_metadata %>%
            pull(seq_type) %>%
            unique()
    }

    if(missing(regions_bed)){
        message(
            "Using aSHM regions in grch37 projection as regions_bed"
        )
        regions_bed <- GAMBLR.data::grch37_ashm_regions %>%
            mutate(name = paste(gene, region, sep = "_"))
    }

    ashm_maf <- get_ssm_by_regions(
        regions_bed = regions_bed,
        streamlined = TRUE,
        maf_data = maf_data,
        use_name_column = TRUE
    )

    ashm_counted <- ashm_maf %>%
        group_by(sample_id, region_name) %>%
        tally()

    if(missing(these_samples_metadata)){
        all_meta <- get_gambl_metadata(
            seq_type_filter=this_seq_type
        ) %>%
        dplyr::select(sample_id)
    }else{
        all_meta <- these_samples_metadata %>%
            dplyr::select(sample_id)
    }
    #fill out all combinations so we can get the cases with zero mutations
    eg <- expand_grid(
        sample_id = pull(all_meta, sample_id),
        region_name = unique(ashm_counted$region_name)
    )
    all_counts <- left_join(eg, ashm_counted) %>%
        mutate(n = replace_na(n, 0)) %>%
        unique() #not sure where the duplicates are coming from but its annoying

    all_counts_wide <- pivot_wider(
        all_counts,
        id_cols = sample_id,
        names_from = region_name,
        values_from = n
    ) %>%
        column_to_rownames(var = "sample_id")

    return(all_counts_wide)
}
