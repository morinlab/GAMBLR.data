#' @title Query SSMs from the mutations database.
#'
#' @description Shared data-access helper that retrieves simple somatic mutations
#' (SSM, MAF format) from `gambl_mutations.db`, pushing the common filters down to
#' indexed SQL. This is the single query path behind the GAMBLR.open SSM
#' accessors (get_coding_ssm, get_ssm_by_samples, get_ssm_by_region(s), ...).
#'
#' @param projection Genome build / `genome_build` column value. Default "grch37".
#' @param sample_ids Optional character vector of `Tumor_Sample_Barcode` to keep.
#' @param tool_name Pipeline to keep, matched case-insensitively (the value
#'   you pass is lowercased before comparing; the stored Pipeline column is
#'   itself normalized to lowercase at write time, so this is a plain,
#'   indexed equality check rather than a function-wrapped one). Default
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
#' @param seq_type Optional character vector to restrict rows to one or more
#'   `maf_seq_type` values (`"genome"`/`"capture"`), so a sample_id with both
#'   can be queried separately instead of getting both back merged together.
#'   Default `NULL` (no filter, matching prior behaviour). Rows predating
#'   `maf_seq_type` (e.g. publication-pipeline rows, which don't carry
#'   per-row seq_type provenance) have `NA` here and are excluded whenever
#'   this filter is set, same as any other `NA` under a `%in%` filter. Silently
#'   ignored against a `gambl_mutations.db` built before this column existed.
#' @param regions Optional data frame with columns `chrom`, `start`, `end`;
#'   rows are OR-ed and applied as a single combined query (not one query per
#'   region), so passing hundreds of regions (e.g. one per gene in a panel)
#'   doesn't cost hundreds of separate round trips.
#' @param con Optional DBI connection (defaults to [gambl_mutations_db()]).
#'
#' @return A data frame of MAF rows (the `genome_build`/`mutation_id` helper
#'   columns are dropped).
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
                            seq_type = NULL,
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
  # variant_pipeline (mutation_id, elem, Pipeline) records every pipeline
  # that actually produced each variant, keyed on the surrogate mutation_id
  # surviving rows are assigned at write time (not the 5-column natural key
  # -- see write_mutations_db.R for why), so tool_name is resolved against
  # it via a single-column semi-join instead of the maf/ashm row's own
  # (possibly-not-representative) Pipeline value. Falls back to the old
  # direct-column filter against DBs built before this table existed.
  has_variant_pipeline <- !is.null(tn) && "variant_pipeline" %in% DBI::dbListTables(con)

  # apply the filters that are common to every query against this table,
  # regardless of how many regions (if any) are requested
  base_query <- function(table_name) {
    q <- dplyr::tbl(con, table_name) %>%
      dplyr::filter(genome_build == projection)
    if (!is.null(tn)) {
      if (has_variant_pipeline) {
        matching_ids <- dplyr::tbl(con, "variant_pipeline") %>%
          dplyr::filter(elem == table_name, Pipeline == tn) %>%
          dplyr::distinct(mutation_id)
        q <- dplyr::semi_join(q, matching_ids, by = "mutation_id")
      } else {
        q <- dplyr::filter(q, Pipeline == tn)
      }
    }
    if (coding_only)             q <- dplyr::filter(q, Variant_Classification %in% cc)
    if (min_read_support > 0)    q <- dplyr::filter(q, t_alt_count >= min_read_support)
    if (!is.null(sample_ids))    q <- dplyr::filter(q, Tumor_Sample_Barcode %in% sample_ids)
    if (!is.null(seq_type) && "maf_seq_type" %in% DBI::dbListFields(con, table_name)) {
      q <- dplyr::filter(q, maf_seq_type %in% seq_type)
    }
    q
  }

  # All regions are combined into a single OR'd condition and applied in one
  # query, not one query per region. get_ssm_by_regions() commonly passes a
  # region per gene in a panel (150-250+ for the full lymphoma gene list) --
  # looping here used to mean that many separate round trips, each one
  # redoing the tool_name/coding_only/sample_ids filters and re-running the
  # variant_pipeline semi-join from scratch for every single region. A single
  # combined query also fixes a latent duplicate-row bug the loop had: a
  # variant landing inside two overlapping regions (e.g. neighbouring genes'
  # padded windows) was previously returned once per matching region and
  # bind_rows()'d into two identical rows; SQL's OR doesn't double-count a
  # row just because more than one disjunct matches it.
  gather <- function(table_name) {
    q <- base_query(table_name)
    if (!is.null(regions) && nrow(regions) > 0) {
      region_exprs <- lapply(seq_len(nrow(regions)), function(i) {
        rc <- as.character(regions$chrom[i])
        rs <- as.numeric(regions$start[i])
        re <- as.numeric(regions$end[i])
        rlang::expr(Chromosome == !!rc & Start_Position > !!rs & Start_Position < !!re)
      })
      combined <- Reduce(function(a, b) rlang::expr(!!a | !!b), region_exprs)
      q <- dplyr::filter(q, !!combined)
    }
    dplyr::collect(q)
  }

  res <- gather("maf")
  if (include_ashm) res <- dplyr::bind_rows(res, gather("ashm"))
  res$genome_build <- NULL
  res$mutation_id <- NULL
  res
}
