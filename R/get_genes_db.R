#' @title Connect to the bundled GAMBLR reference database.
#'
#' @description Returns a (cached) read-only DBI connection to the normalized
#' SQLite reference database built from the LLMPP curated gene lists. This is the
#' shared backing store used by [get_genes()] (and intended for the Lymphopedia
#' wiki / data API), replacing the versioned `lymphoma_genes_*` .rda objects.
#'
#' @details The database is rebuilt from source by `data-raw/build_reference_db.R`
#' and is small (well under 1 MB), so — unlike [gambl_mutations_db()] — it is
#' bundled inside the package (`inst/extdata/gambl_reference.db`).
#'
#' # Table reference
#'
#' | Table | Grain | Key columns |
#' | --- | --- | --- |
#' | `gene_entity` | one row per gene x lymphoma entity | `gene, ensembl_gene_id, entity, tier, ashm, qc, citekey, pmid, mutation_effect, extra` (JSON of entity-specific fields), `source_version` |
#' | `ashm_regions` | one row per gene x aSHM region x genome build | `gene, chrom, start, end, region, regulatory_comment, build, source_version` |
#' | `build_info` | key/value | provenance: `source` (LLMPP resources/curated), `source_version`, `built_at`, `builder` |
#'
#' `gene_entity$entity` takes values like `"DLBCL"`, `"BL"`, `"FL"`, `"PMBL"`,
#' `"MCL"`, `"MZL"`; `tier` is `1` (high-confidence), `2` (lower-confidence), or
#' `3` (retired). `ashm_regions$build` is `"grch37"` or `"hg38"`.
#'
#' @param db_path Optional explicit path to the .db file. Defaults to the copy
#'   installed with the package; falls back to `inst/extdata/gambl_reference.db`
#'   during development.
#'
#' @return A DBIConnection.
#'
#' @examples
#' \dontrun{
#' con <- gambl_reference_db()
#' dplyr::tbl(con, "gene_entity")
#'
#' # high-confidence DLBCL genes
#' dplyr::tbl(con, "gene_entity") |>
#'   dplyr::filter(entity == "DLBCL", tier == 1) |>
#'   dplyr::collect()
#' }
#' @export
gambl_reference_db <- function(db_path = NULL) {
  if (is.null(db_path)) {
    db_path <- system.file("extdata", "gambl_reference.db", package = "GAMBLR.data")
    if (!nzchar(db_path) || !file.exists(db_path)) {
      db_path <- file.path("inst", "extdata", "gambl_reference.db")
    }
  }
  stopifnot(file.exists(db_path))
  cache <- getOption("gamblr.reference.con")
  if (!is.null(cache) && DBI::dbIsValid(cache)) return(cache)
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path, flags = RSQLite::SQLITE_RO)
  options(gamblr.reference.con = con)
  con
}
