#' @title Annotate Hotspots.
#'
#' @description Annotate MAF-like data frome with a hot_spot column indicating
#' recurrent mutations.
#'
#' @details This function takes an already loaded MAF data frame with the
#' `mutation_maf` parameter.
#'
#' @param mutation_maf A data frame in MAF format.
#' @param ... Any other parameter. These parameters will be ignored.
#'
#' @return The same data frame with one additional column "hot_spot".
#'
#' @import dplyr
#' @export
#'
#' @examples
#' my_metadata = get_gambl_metadata()
#' all_coding_ssm = get_coding_ssm(these_samples_metadata = my_metadata,
#'                                 projection = "grch37",
#'                                 this_seq_type = "genome") %>%
#'                    dplyr::filter(Hugo_Symbol %in% c("EZH2",
#'                                  "MEF2B","MYD88","KMT2D")) %>%
#'                    dplyr::arrange(Hugo_Symbol)
#'
#' hot_ssms = annotate_hotspots(all_coding_ssm)
#' hot_ssms %>% dplyr::filter(!is.na(hot_spot)) %>%
#'       dplyr::select(1:5,37,hot_spot)
#'
annotate_hotspots = function(
  mutation_maf,
  ...
) {

  # check if any invalid parameters are provided
  check_excess_params(...)

  filled_coords <- GAMBLR.data::hotspots_annotations
  # just the ssms that match these coordinates!
  hot_ssms <- left_join(
    mutation_maf,
    filled_coords,
    by = c("Chromosome", "Start_Position")
  )
  return(hot_ssms)
}
