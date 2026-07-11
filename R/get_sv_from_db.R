#' @title Query structural variants (SV) from the mutations database.
#'
#' @description Data-access helper that retrieves Manta SV breakpoints (BEDPE
#' format) from the `bedpe` table of `gambl_mutations.db`, pushing the common
#' filters down to SQL. Backs the GAMBLR.open `get_manta_sv` accessor.
#'
#' @param projection Genome build / `genome_build` column value. Default "grch37".
#' @param sample_ids Optional character vector of `tumour_sample_id` to keep.
#' @param min_vaf Minimum `VAF_tumour`. Default 0.
#' @param min_score Minimum Manta `SCORE`. Default 0.
#' @param pass_only When TRUE, keep only `FILTER == "PASS"` breakpoints.
#' @param region Optional list with `chrom`, `start`, `end`; keeps breakpoints
#'   whose A- or B-end falls in the range (positions compared numerically).
#' @param con Optional DBI connection (defaults to [gambl_mutations_db()]).
#'
#' @return A data frame of BEDPE rows (the `genome_build` helper column is dropped).
#'
#' @export
get_sv_from_db <- function(projection = "grch37",
                           sample_ids = NULL,
                           min_vaf = 0,
                           min_score = 0,
                           pass_only = FALSE,
                           region = NULL,
                           con = NULL) {
  if (is.null(con)) con <- gambl_mutations_db()
  q <- dplyr::tbl(con, "bedpe") %>%
    dplyr::filter(genome_build == projection,
                  VAF_tumour >= min_vaf,
                  SCORE >= min_score)
  if (!is.null(sample_ids)) q <- dplyr::filter(q, tumour_sample_id %in% sample_ids)
  if (pass_only)            q <- dplyr::filter(q, FILTER == "PASS")
  if (!is.null(region)) {
    rc <- region$chrom; rs <- as.numeric(region$start); re <- as.numeric(region$end)
    q <- dplyr::filter(q, (CHROM_A == rc & START_A >= rs & START_A <= re) |
                          (CHROM_B == rc & START_B >= rs & START_B <= re))
  }
  res <- dplyr::collect(q)
  res$genome_build <- NULL
  res
}
