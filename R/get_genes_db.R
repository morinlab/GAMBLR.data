#' @title Connect to the bundled GAMBLR reference database.
#'
#' @description Returns a (cached) read-only DBI connection to the normalized
#' SQLite reference database built from the LLMPP curated gene lists. This is the
#' shared backing store used by [get_genes()] (and intended for the Lymphopedia
#' wiki / data API), replacing the versioned `lymphoma_genes_*` .rda objects.
#'
#' @details The database is rebuilt from source by `data-raw/build_reference_db.R`.
#' Tables: `gene_entity` (one row per gene x entity, with tier / aSHM / citekey /
#' ensembl_gene_id / source_version), `ashm_regions`, and `build_info`.
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
