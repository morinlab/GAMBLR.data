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
#' Not every package version gets its own data release. When `tag` isn't
#' given and the release for this package version has no
#' `gambl_mutations.db` (or doesn't exist), the newest `data-v*` release at
#' or below this package version that has one is used instead, with a
#' message saying which. A `tag` given explicitly is used as is.
#'
#' @param dest_path Where to save the downloaded file. Defaults to the
#'   standard cache location, `tools::R_user_dir("GAMBLR.data", "cache")`.
#' @param tag GitHub release tag to download from. Defaults to
#'   `paste0("data-v", utils::packageVersion("GAMBLR.data"))`, or, if that
#'   release has no `gambl_mutations.db`, the newest earlier data release
#'   that does (see details).
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
  tag_is_default <- missing(tag)

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

  # default tag: this version's data release, or the newest earlier one with a database (see details)
  if (tag_is_default) {
    tag <- resolve_data_release_tag(utils::packageVersion("GAMBLR.data"), repo = repo)
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

# The data release to download for a package version: "data-v<version>" if that release has a gambl_mutations.db,
# otherwise the newest "data-v*" release at or below `version` that has one. `release_assets` (a data frame of tag
# and file_name, as piggyback::pb_list() returns) defaults to the repo's current assets.
resolve_data_release_tag <- function(version, repo = "morinlab/GAMBLR.data", release_assets = NULL) {
  if (is.null(release_assets)) release_assets <- piggyback::pb_list(repo = repo)
  with_db <- unique(as.character(release_assets$tag[release_assets$file_name == "gambl_mutations.db"]))
  wanted <- paste0("data-v", version)
  if (wanted %in% with_db) return(wanted)
  data_tags <- grep("^data-v[0-9]+(\\.[0-9]+)*$", with_db, value = TRUE)
  versions <- numeric_version(sub("^data-v", "", data_tags))
  earlier <- which(versions <= numeric_version(as.character(version)))
  if (!length(earlier)) {
    stop("No data release with gambl_mutations.db for GAMBLR.data ", version, " or earlier in ", repo,
         " (releases with one: ", if (length(with_db)) paste(with_db, collapse = ", ") else "none", ").", call. = FALSE)
  }
  best <- data_tags[earlier][order(versions[earlier], decreasing = TRUE)[1]]
  message("No gambl_mutations.db released for GAMBLR.data ", version, " (", wanted, "); using the newest earlier data release, ", best, ".")
  best
}
