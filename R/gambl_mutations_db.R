#' @title Connect to the GAMBL mutations database.
#'
#' @description Returns a (cached) read-only DBI connection to `gambl_mutations.db`,
#' the SQLite database holding the large sample-level frames (MAF, aSHM, seg,
#' bedpe) built from the in-memory `sample_data` object assembled by
#' `data-raw/assemble_bundled_data.R`, via the shared `write_mutations_db()`
#' helper (`data-raw/write_mutations_db.R`).
#'
#' @details Unlike the small bundled [gambl_reference_db()], this file is large
#' (several hundred MB) and is NOT shipped inside the package. It is
#' distributed as a GitHub Release asset (see [download_gambl_mutations_db()])
#' and cached locally. The path is resolved in this order:
#' \enumerate{
#'   \item the `db_path` argument;
#'   \item `getOption("GAMBLR.data.mutations_db")`;
#'   \item the `GAMBLR_MUTATIONS_DB` environment variable;
#'   \item the user cache dir, `tools::R_user_dir("GAMBLR.data", "cache")`;
#'   \item `gambl_mutations.db` in the working directory (development).
#' }
#' If none of these exist and `auto_download` is `TRUE` (the default), the
#' file is downloaded automatically into the user cache dir via
#' [download_gambl_mutations_db()]. Set `auto_download = FALSE`, or
#' `options(GAMBLR.data.auto_download = FALSE)` to disable this globally
#' (e.g. offline/CI environments) and get an explicit error instead.
#'
#' # Table reference
#'
#' Both genome builds are stacked into the same table (not split across
#' separate tables); filter on `genome_build` (`"grch37"` / `"hg38"`) to select
#' one. `maf` and `ashm` share the same (curated, 49-column) schema: a reduced
#' column set, not the full upstream flatfile MAF.
#'
#' | Table | Rows (typical) | Grain | Key columns |
#' | --- | --- | --- | --- |
#' | `maf` | ~3.6M | one somatic mutation call | `mutation_id` (surrogate integer key, unique within `maf` across both genome builds -- assigned at write time, after deduplication; `variant_pipeline` references it as a foreign key), `Hugo_Symbol, Chromosome, Start_Position, End_Position, Tumor_Sample_Barcode, Variant_Classification, HGVSp_Short, t_alt_count, n_alt_count, Pipeline (lowercase, e.g. `"slms-3"`/`"publication"` -- normalized at write time so it can be matched with a plain, indexed equality), genome_build` (+ ~35 more MAF-standard columns) |
#' | `ashm` | ~135k | one mutation call in an aSHM region | same schema as `maf` (its own independent `mutation_id` sequence) |
#' | `seg` | ~126k | one copy-number segment | `ID, chrom, start, end, LOH_flag, log.ratio, CN, genome_build` |
#' | `bedpe` | ~900 | one Manta structural-variant breakpoint pair | `CHROM_A, START_A, END_A, CHROM_B, START_B, END_B, manta_name, SCORE, STRAND_A, STRAND_B, tumour_sample_id, normal_sample_id, VAF_tumour, DP, pair_status, FILTER, genome_build` |
#' | `sample_meta` | ~3.3k | one sample | `patient_id, sample_id, Tumor_Sample_Barcode, seq_type, pathology, cohort, study, ...` (its own `genome_build` column records the sample's native alignment build, unrelated to the per-row build stamp used in the other tables) |
#' | `sample_study` | varies | one (sample, study) membership fact -- many-to-many, a sample belonging to N studies is N rows | `sample_id, study, study_id, reference_PMID`. `study_id` is that study's own identifier for the sample where it differs from GAMBL's `sample_id` (e.g. a paper's own case/patient ID); NA where no study-specific ID has been sourced. |
#' | `variant_pipeline` | varies | one (mutation, Pipeline) fact -- many-to-many, a variant independently called by N pipelines is N rows | `mutation_id, elem, Pipeline`. `mutation_id` is a foreign key to `maf.mutation_id` or `ashm.mutation_id` depending on `elem` (`"maf"` or `"ashm"`) -- not a natural key, so this table doesn't repeat the variant's own position/sample columns. |
#' | `build_info` | 10 | key/value | provenance (`source`, `built_at`, `builder`) and expected row counts (`n_maf`, `n_ashm`, `n_seg`, `n_bedpe`, `n_samples`, `n_sample_study`, `n_variant_pipeline`) used by the `data-raw/test_gambl_db.R` regression checks |
#'
#' `maf`/`ashm` do NOT carry a `Study` column. Cohort/study membership is
#' tracked once per sample in `sample_study`, not once per mutation row --
#' a mutation row's `Study` couldn't represent a sample belonging to more
#' than one study anyway, and a per-row tag caused the same sample's
#' mutations to be pulled and duplicated once per cohort-specific code path
#' that claimed it (see `GAMBLR.data::get_ssm_from_db()`'s `this_study`
#' parameter for how to filter by study post-refactor).
#'
#' Similarly, `maf`/`ashm` only ever carry ONE `Pipeline` value per variant,
#' even though the same real mutation can be independently produced by more
#' than one pipeline (e.g. a cohort's own published/curated maf and our
#' separate SLMS-3 recall both calling the same position) -- when that
#' happens, `write_mutations_db()`'s write-time dedup keeps one row
#' arbitrarily (whichever pipeline's data was assembled first), which would
#' otherwise make the variant invisible to a query for the *other* pipeline
#' even though it's still fully present in the table. `variant_pipeline`
#' records every pipeline that actually produced each variant, independent
#' of which one "won" as `maf`/`ashm`'s single representative row; see
#' `GAMBLR.data::get_ssm_from_db()`'s `tool_name` parameter, which resolves
#' against this table rather than `maf`/`ashm`'s own `Pipeline` column.
#'
#' ## Join keys
#' `maf`, `ashm`, and `bedpe` (via `tumour_sample_id`) key to
#' `sample_meta$Tumor_Sample_Barcode` / `sample_meta$sample_id`; `seg$ID` keys
#' to `sample_meta$sample_id`; `sample_study$sample_id` keys to
#' `sample_meta$sample_id`; `variant_pipeline$mutation_id` keys to
#' `maf.mutation_id` or `ashm.mutation_id`, selecting which one via
#' `variant_pipeline$elem`.
#'
#' ## Indexes
#' `maf`/`ashm`: `(genome_build, Chromosome, Start_Position)`,
#' `Tumor_Sample_Barcode`, `Pipeline`, `mutation_id`. No index on
#' `Hugo_Symbol`: gene-restricted queries resolve the gene to a region first
#' (the same logic used to populate these tables) and filter on
#' `(genome_build, Chromosome, Start_Position)` instead -- see
#' `GAMBLR.data::get_ssm_from_db()`. `seg`: `(genome_build, ID)`,
#' `(genome_build, chrom, start)`. `bedpe`: `(tumour_sample_id,
#' genome_build)`, `(genome_build, CHROM_A, START_A)`, `(genome_build,
#' CHROM_B, START_B)` (one per breakpoint end, so a region search that OR's
#' both ends can use a different index per side).
#' `sample_meta`: `sample_id`, `Tumor_Sample_Barcode`. `sample_study`:
#' `sample_id`, `study`. `variant_pipeline`: `mutation_id`,
#' `(elem, Pipeline)`.
#'
#' @param db_path Optional explicit path to the .db file.
#' @param auto_download When no existing copy is found, download it
#'   automatically via [download_gambl_mutations_db()]. Default `TRUE`;
#'   can also be disabled globally with
#'   `options(GAMBLR.data.auto_download = FALSE)`.
#'
#' @return A read-only DBIConnection.
#'
#' @examples
#' \dontrun{
#' con <- gambl_mutations_db()
#' DBI::dbListTables(con)
#'
#' # coding mutations in TP53, grch37, one pipeline (Pipeline is stored
#' # lowercase -- see the Pipeline column note above)
#' dplyr::tbl(con, "maf") |>
#'   dplyr::filter(genome_build == "grch37", Hugo_Symbol == "TP53",
#'                Pipeline == "slms-3") |>
#'   dplyr::collect()
#' }
#' @export
gambl_mutations_db <- function(db_path = NULL, auto_download = TRUE) {
  if (is.null(db_path)) {
    opt <- getOption("GAMBLR.data.mutations_db", default = NULL)
    env <- Sys.getenv("GAMBLR_MUTATIONS_DB", unset = "")
    cache_path <- file.path(tools::R_user_dir("GAMBLR.data", "cache"), "gambl_mutations.db")
    candidates <- c(
      opt,
      if (nzchar(env)) env,
      cache_path,
      "gambl_mutations.db"
    )
    candidates <- candidates[!is.null(candidates)]
    hit <- candidates[file.exists(candidates)]
    if (length(hit) == 0) {
      can_auto_download <- isTRUE(auto_download) && isTRUE(getOption("GAMBLR.data.auto_download", TRUE))
      if (can_auto_download) {
        db_path <- download_gambl_mutations_db(dest_path = cache_path)
      } else {
        stop("gambl_mutations.db not found. Set options(GAMBLR.data.mutations_db = \"/path\"), ",
             "the GAMBLR_MUTATIONS_DB env var, place it in ",
             tools::R_user_dir("GAMBLR.data", "cache"),
             ", or call gambl_mutations_db() with auto_download = TRUE (the default) ",
             "to fetch it automatically.", call. = FALSE)
      }
    } else {
      db_path <- hit[1]
    }
  }
  stopifnot(file.exists(db_path))

  cached <- getOption("gamblr.mutations.con")
  if (!is.null(cached) && DBI::dbIsValid(cached) &&
      identical(normalizePath(cached@dbname), normalizePath(db_path))) {
    return(cached)
  }
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path, flags = RSQLite::SQLITE_RO)
  options(gamblr.mutations.con = con)
  con
}
