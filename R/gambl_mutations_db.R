#' @title Connect to the GAMBL mutations database.
#'
#' @description Returns a (cached) read-only DBI connection to `gambl_mutations.db`,
#' the SQLite database holding the large sample-level frames (MAF, aSHM, seg,
#' bedpe) built from `sample_data` by `data-raw/build_mutations_db.R`.
#'
#' @details Unlike the small bundled reference database, this file is large
#' (~1 GB) and is NOT shipped inside the package. It is distributed as a cached
#' release asset. The path is resolved in this order:
#' \enumerate{
#'   \item the `db_path` argument;
#'   \item `getOption("GAMBLR.data.mutations_db")`;
#'   \item the `GAMBLR_MUTATIONS_DB` environment variable;
#'   \item the user cache dir, `tools::R_user_dir("GAMBLR.data", "cache")`;
#'   \item `gambl_mutations.db` in the working directory (development).
#' }
#'
#' @param db_path Optional explicit path to the .db file.
#'
#' @return A read-only DBIConnection.
#'
#' @examples
#' \dontrun{
#' con <- gambl_mutations_db()
#' DBI::dbListTables(con)
#' }
#' @export
gambl_mutations_db <- function(db_path = NULL) {
  if (is.null(db_path)) {
    opt <- getOption("GAMBLR.data.mutations_db", default = NULL)
    env <- Sys.getenv("GAMBLR_MUTATIONS_DB", unset = "")
    candidates <- c(
      opt,
      if (nzchar(env)) env,
      file.path(tools::R_user_dir("GAMBLR.data", "cache"), "gambl_mutations.db"),
      "gambl_mutations.db"
    )
    candidates <- candidates[!is.null(candidates)]
    hit <- candidates[file.exists(candidates)]
    if (length(hit) == 0) {
      stop("gambl_mutations.db not found. Set options(GAMBLR.data.mutations_db = \"/path\"), ",
           "the GAMBLR_MUTATIONS_DB env var, or place it in ",
           tools::R_user_dir("GAMBLR.data", "cache"), call. = FALSE)
    }
    db_path <- hit[1]
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
