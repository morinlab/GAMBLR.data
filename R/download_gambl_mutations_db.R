#' @title Download the GAMBL mutations database.
#'
#' @description Downloads `gambl_mutations.db` from this package's GitHub
#' Release assets, if not already present. Normally called automatically by
#' [gambl_mutations_db()] the first time it's needed; call it directly to
#' pre-fetch the file (e.g. before a batch job on a machine with no other
#' network access) or to force a re-download.
#'
#' @details Requires the `piggyback` package (Suggests, not a hard
#' dependency of GAMBLR.data -- most calls to [gambl_mutations_db()] never
#' need to download anything). Downloads are tagged per package version by
#' default (`"data-v<version>"`), so a given install always fetches data
#' built for its own schema, not whatever happens to be the newest release
#' -- see `data-raw/release_mutations_db.R` for the corresponding upload
#' step (admin-only, run manually after a fresh build; bump
#' `DESCRIPTION`'s `Version` and cut a matching release whenever the schema
#' changes, so the two stay in lockstep).
#'
#' @param dest_path Where to save the downloaded file. Defaults to the
#'   standard cache location, `tools::R_user_dir("GAMBLR.data", "cache")`.
#' @param tag GitHub release tag to download from. Defaults to
#'   `paste0("data-v", utils::packageVersion("GAMBLR.data"))`.
#' @param repo GitHub repo in `"owner/repo"` form.
#' @param overwrite Re-download even if `dest_path` already exists. Default
#'   `FALSE`.
#'
#' @return `dest_path`, invisibly.
#'
#' @export
download_gambl_mutations_db <- function(
    dest_path = file.path(tools::R_user_dir("GAMBLR.data", "cache"), "gambl_mutations.db"),
    tag = paste0("data-v", utils::packageVersion("GAMBLR.data")),
    repo = "morinlab/GAMBLR.data",
    overwrite = FALSE) {

  if (file.exists(dest_path) && !overwrite) {
    message("gambl_mutations.db already exists at ", dest_path,
            "; set overwrite = TRUE to re-download.")
    return(invisible(dest_path))
  }

  if (!requireNamespace("piggyback", quietly = TRUE)) {
    stop("Downloading gambl_mutations.db requires the 'piggyback' package. ",
         "Install it with install.packages(\"piggyback\"), or point ",
         "gambl_mutations_db() at an existing copy via the db_path argument, ",
         "options(GAMBLR.data.mutations_db = ...), or the GAMBLR_MUTATIONS_DB ",
         "environment variable.", call. = FALSE)
  }

  dir.create(dirname(dest_path), recursive = TRUE, showWarnings = FALSE)
  message(sprintf(
    "Downloading gambl_mutations.db from %s @ %s into %s -- this is a large file (several hundred MB) and may take a while...",
    repo, tag, dirname(dest_path)
  ))
  piggyback::pb_download(
    "gambl_mutations.db",
    repo = repo,
    tag = tag,
    dest = dirname(dest_path),
    overwrite = TRUE
  )
  if (!file.exists(dest_path)) {
    stop("Download appeared to succeed but ", dest_path, " still doesn't exist -- ",
         "check that the release tag \"", tag, "\" actually has a gambl_mutations.db asset.",
         call. = FALSE)
  }
  message("Done: ", dest_path)
  invisible(dest_path)
}
