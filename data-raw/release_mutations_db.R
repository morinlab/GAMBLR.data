# release_mutations_db.R
#
# Admin-only: uploads a freshly-built gambl_mutations.db as a GitHub Release
# asset, so GAMBLR.data::gambl_mutations_db() can auto-download it for new
# installs (see R/download_gambl_mutations_db.R). Run this AFTER
# assemble_bundled_data.R / write_mutations_db() has produced a fresh
# gambl_mutations.db, from the GAMBLR.data repo root. This script only
# uploads -- it does not build the database itself.
#
# Requires:
#  - the piggyback package: install.packages("piggyback")
#  - a GitHub PAT with write access to this repo, either already configured
#    (gitcreds::gitcreds_set()) or set via the GITHUB_PAT/GITHUB_TOKEN
#    environment variable
#
# The release tag defaults to "data-v<DESCRIPTION Version>", matching what
# download_gambl_mutations_db() looks for by default. Keep these in lockstep:
# bump DESCRIPTION's Version whenever the DB schema changes (as it did
# several times in the sample_study/variant_pipeline refactor), and cut a
# matching release here, so a given package version always downloads data
# built for its own schema rather than whatever the newest release happens
# to be.
#
# Usage: Rscript data-raw/release_mutations_db.R [path to .db] [tag]

args <- commandArgs(trailingOnly = TRUE)
db_path <- if (length(args) >= 1) args[[1]] else "gambl_mutations.db"
pkg_version <- read.dcf("DESCRIPTION")[, "Version"]
tag <- if (length(args) >= 2) args[[2]] else paste0("data-v", pkg_version)
repo <- "morinlab/GAMBLR.data"

stopifnot(file.exists(db_path))
if (!requireNamespace("piggyback", quietly = TRUE)) {
  stop("Install piggyback first: install.packages(\"piggyback\")")
}

message(sprintf("Releasing %s to %s @ %s (%.0f MB)",
                db_path, repo, tag, file.info(db_path)$size / 1024^2))

# create the release/tag if it doesn't already exist
existing_releases <- piggyback::pb_releases(repo = repo)
if (!tag %in% existing_releases$tag_name) {
  message("Creating new release tag: ", tag)
  piggyback::pb_new_release(repo = repo, tag = tag)
}

piggyback::pb_upload(db_path, repo = repo, tag = tag)
message(sprintf(
  "Done. Users on GAMBLR.data v%s will now auto-download gambl_mutations.db from this release.",
  pkg_version
))
