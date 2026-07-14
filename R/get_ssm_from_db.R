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
#'   "slms-3"; set NULL to skip the Pipeline filter. Resolved via the
#'   `variant_pipeline` join table (see [gambl_mutations_db()]'s schema
#'   docs) rather than `maf`/`ashm`'s own Pipeline column, so a variant
#'   independently called by more than one pipeline is still found by a
#'   query for either one, even though its single maf/ashm row can only
#'   carry one Pipeline value. Falls back to filtering the column directly
#'   against DBs built before variant_pipeline existed.
#' @param include_ashm When TRUE, also query the `ashm` table and row-bind it.
#' @param coding_only When TRUE, keep only coding `Variant_Classification`s.
#' @param include_silent When FALSE (and `coding_only`), drop Silent mutations.
#' @param min_read_support Keep only variants with `t_alt_count` >= this value.
#' @param this_study Optional single study to restrict to, matched against
#'   `sample_study$study` (not a `maf`/`ashm` column -- see
#'   [gambl_mutations_db()]'s schema docs). Resolves to the set of sample_ids
#'   belonging to that study, intersected with `sample_ids` if both are
#'   supplied; composes with `tool_name`/Pipeline exactly as before (e.g.
#'   `this_study="Reddy"` alone returns Reddy's SLMS-3 recall by default,
#'   while `this_study="Reddy", tool_name="publication"` returns Reddy's
#'   as-published rows only).
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

  # this_study resolves to sample_ids via a join against sample_study, done
  # once here rather than as a per-query filter inside base_query() -- reuses
  # the existing sample_ids/Tumor_Sample_Barcode machinery below instead of
  # requiring maf/ashm to carry their own Study column (they don't; see
  # gambl_mutations_db()'s schema docs), and keeps this_study (which
  # samples) and tool_name/Pipeline (which rows for those samples) composable
  # exactly as they were.
  if (!is.null(this_study)) {
    study_sample_ids <- dplyr::tbl(con, "sample_study") %>%
      dplyr::filter(study == this_study) %>%
      dplyr::distinct(sample_id) %>%
      dplyr::pull(sample_id)
    sample_ids <- if (is.null(sample_ids)) study_sample_ids
                  else intersect(sample_ids, study_sample_ids)
  }

  # A variant that's independently produced by more than one pipeline (e.g. a
  # cohort's own published maf and our SLMS-3 recall both calling the same
  # position) only has room for one Pipeline value on its single maf/ashm
  # row -- write_mutations_db() picks one arbitrarily when it deduplicates.
  # variant_pipeline (sample_id/Chromosome/Start_Position/End_Position/
  # Tumor_Seq_Allele2/genome_build, Pipeline) records every pipeline that
  # actually produced each variant, so tool_name is resolved against it via a
  # semi-join instead of the maf/ashm row's own (possibly-not-representative)
  # Pipeline value. Falls back to the old direct-column filter against DBs
  # built before this table existed.
  has_variant_pipeline <- !is.null(tn) && "variant_pipeline" %in% DBI::dbListTables(con)
  variant_key_cols <- c("Tumor_Sample_Barcode", "Chromosome", "Start_Position", "End_Position", "Tumor_Seq_Allele2")

  # apply the filters that are common to a single (table, region) query
  base_query <- function(table_name, region = NULL) {
    q <- dplyr::tbl(con, table_name) %>%
      dplyr::filter(genome_build == projection)
    if (!is.null(tn)) {
      if (has_variant_pipeline) {
        matching_keys <- dplyr::tbl(con, "variant_pipeline") %>%
          dplyr::filter(elem == table_name, genome_build == projection, tolower(Pipeline) == tn) %>%
          dplyr::distinct(dplyr::across(dplyr::all_of(variant_key_cols)))
        q <- dplyr::semi_join(q, matching_keys, by = variant_key_cols)
      } else {
        q <- dplyr::filter(q, tolower(Pipeline) == tn)
      }
    }
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
