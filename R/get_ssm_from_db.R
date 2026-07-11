#' @title Query SSMs from the mutations database.
#'
#' @description Shared data-access helper that retrieves simple somatic mutations
#' (SSM, MAF format) from `gambl_mutations.db`, pushing the common filters down to
#' indexed SQL. This is the single query path behind the GAMBLR.open SSM
#' accessors (get_coding_ssm, get_ssm_by_samples, get_ssm_by_region(s), ...).
#'
#' @param projection Genome build / `genome_build` column value. Default "grch37".
#' @param sample_ids Optional character vector of `Tumor_Sample_Barcode` to keep.
#' @param tool_name Pipeline to keep (matched case-insensitively). Default
#'   "slms-3"; set NULL to skip the Pipeline filter.
#' @param include_ashm When TRUE, also query the `ashm` table and row-bind it.
#' @param coding_only When TRUE, keep only coding `Variant_Classification`s.
#' @param include_silent When FALSE (and `coding_only`), drop Silent mutations.
#' @param min_read_support Keep only variants with `t_alt_count` >= this value.
#' @param this_study Optional single `Study` to restrict to.
#' @param regions Optional data frame with columns `chrom`, `start`, `end`; rows
#'   are OR-ed (each region is an indexed range scan).
#' @param con Optional DBI connection (defaults to [gambl_mutations_db()]).
#'
#' @return A data frame of MAF rows (the `genome_build` helper column is dropped).
#'
#' @export
get_ssm_from_db <- function(projection = "grch37",
                            sample_ids = NULL,
                            tool_name = "slms-3",
                            include_ashm = FALSE,
                            coding_only = FALSE,
                            include_silent = TRUE,
                            min_read_support = 0,
                            this_study = NULL,
                            regions = NULL,
                            con = NULL) {
  if (is.null(con)) con <- gambl_mutations_db()
  cc <- if (include_silent) coding_class else coding_class[coding_class != "Silent"]
  tn <- if (!is.null(tool_name)) tolower(tool_name) else NULL

  # apply the filters that are common to a single (table, region) query
  base_query <- function(table_name, region = NULL) {
    q <- dplyr::tbl(con, table_name) %>%
      dplyr::filter(genome_build == projection)
    if (!is.null(tn))            q <- dplyr::filter(q, tolower(Pipeline) == tn)
    if (!is.null(this_study))    q <- dplyr::filter(q, Study == this_study)
    if (coding_only)             q <- dplyr::filter(q, Variant_Classification %in% cc)
    if (min_read_support > 0)    q <- dplyr::filter(q, t_alt_count >= min_read_support)
    if (!is.null(sample_ids))    q <- dplyr::filter(q, Tumor_Sample_Barcode %in% sample_ids)
    if (!is.null(region)) {
      rc <- region$chrom; rs <- region$start; re <- region$end
      q <- dplyr::filter(q, Chromosome == rc &
                            Start_Position > rs & Start_Position < re)
    }
    dplyr::collect(q)
  }

  # one query per region (each uses the (genome_build, Chromosome, Start_Position)
  # index), otherwise a single unrestricted query
  gather <- function(table_name) {
    if (is.null(regions) || nrow(regions) == 0) return(base_query(table_name))
    dplyr::bind_rows(lapply(seq_len(nrow(regions)), function(i)
      base_query(table_name, region = list(
        chrom = as.character(regions$chrom[i]),
        start = as.numeric(regions$start[i]),
        end   = as.numeric(regions$end[i])))))
  }

  res <- gather("maf")
  if (include_ashm) res <- dplyr::bind_rows(res, gather("ashm"))
  res$genome_build <- NULL
  res
}
