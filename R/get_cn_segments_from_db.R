#' @title Query copy-number segments from the mutations database.
#'
#' @description Data-access helper that retrieves copy-number segments from the
#' `seg` table of `gambl_mutations.db`. Backs the GAMBLR.open `get_cn_segments`
#' and `get_sample_cn_segments` accessors.
#'
#' @param projection Genome build / `genome_build` column value. Default "grch37".
#' @param sample_ids Optional character vector of segment `ID`s (sample ids) to keep.
#' @param con Optional DBI connection (defaults to [gambl_mutations_db()]).
#'
#' @return A data frame of seg rows (the `genome_build` helper column is dropped;
#'   chromosome-prefix handling is left to the caller, matching the bundled seg).
#'
#' @export
get_cn_segments_from_db <- function(projection = "grch37",
                                    sample_ids = NULL,
                                    con = NULL) {
  if (is.null(con)) con <- gambl_mutations_db()
  q <- dplyr::tbl(con, "seg") %>%
    dplyr::filter(genome_build == projection)
  if (!is.null(sample_ids)) q <- dplyr::filter(q, ID %in% sample_ids)
  res <- dplyr::collect(q)
  res$genome_build <- NULL
  res
}
