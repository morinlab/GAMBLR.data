# Releasing a new `gambl_mutations.db`

This is the maintainer-only process for publishing a freshly-built
`gambl_mutations.db` so that `GAMBLR.data::gambl_mutations_db()` can
auto-download it for new installs. It's not part of the main `README.md`
because it's not something package users ever need to do -- only whoever
is building and shipping a new version of the bundled data.

## Background

`gambl_mutations.db` is too large (several hundred MB) to ship inside the
package itself. Instead it's distributed as a GitHub Release asset and
downloaded on demand:

- `data-raw/assemble_bundled_data.R` (+ `data-raw/write_mutations_db.R`)
  builds the file from GSC data.
- `data-raw/release_mutations_db.R` (this doc) uploads it to a GitHub
  Release.
- `R/download_gambl_mutations_db.R` downloads it back down for any install
  that doesn't already have a copy -- called automatically by
  `R/gambl_mutations_db.R` the first time it's needed.

Releases are tagged `data-v<DESCRIPTION Version>`. This is deliberate:
`download_gambl_mutations_db()` resolves the tag to fetch from
`utils::packageVersion("GAMBLR.data")` on the *installing* machine, so a
given install always gets data built for its own schema -- not whatever
the newest release happens to be. This means **a schema change requires
both a version bump and a new release**, not just a new upload.

## Steps

1. **Build a fresh `gambl_mutations.db`** via `assemble_bundled_data.R` (on
   the GSC). Run `data-raw/test_gambl_db.R` and `data-raw/compare_bundle_changes.R`
   against it first -- don't release something that hasn't passed those.

2. **Bump `DESCRIPTION`'s `Version`.** Required any time the DB schema
   changed since the last release (new/changed tables, columns, or
   indexes) -- if in doubt, bump it. The release tag is derived from this.

3. **Commit and push the version bump.** `release_mutations_db.R` will
   refuse to run otherwise. The reason: the tag it creates is only
   findable by other installs if the version that generated it is what's
   actually committed and installed elsewhere -- an uncommitted or
   unpushed bump would create a release nothing can ever resolve to,
   including your own next reinstall of this checkout.

4. **Install `piggyback`** if you haven't already: `install.packages("piggyback")`.

5. **Make sure you have a GitHub PAT with write access to this repo.**
   Either already configured (`gitcreds::gitcreds_set()`) or set via the
   `GITHUB_PAT`/`GITHUB_TOKEN` environment variable. This is the step most
   likely to trip you up if it's not already set -- `piggyback::pb_upload()`
   fails with an auth error without it, not a helpful one.

6. **Run the release script from the repo root:**
   ```bash
   Rscript data-raw/release_mutations_db.R /path/to/gambl_mutations.db
   ```
   This creates the `data-v<version>` release/tag if it doesn't already
   exist and uploads the file to it. There's no dry-run mode -- double
   check the path argument before running, since this creates a real,
   public release with a several-hundred-MB asset attached.

7. **Verify the download side works** by clearing your local cache
   (`unlink(file.path(tools::R_user_dir("GAMBLR.data", "cache")), recursive = TRUE)`)
   and calling `GAMBLR.data::gambl_mutations_db()` fresh -- it should
   report downloading, then connect successfully.

## Troubleshooting

- **"DESCRIPTION has uncommitted changes"** -- commit and push the version
  bump first (step 3).
- **Auth error from `piggyback::pb_upload()`/`pb_new_release()`** -- your
  GitHub PAT isn't configured, or doesn't have write access to this repo
  (step 5).
- **A user reports a download that 404s or an outdated schema** -- check
  that their installed `packageVersion("GAMBLR.data")` actually matches a
  tag that exists (`piggyback::pb_releases(repo = "morinlab/GAMBLR.data")`);
  usually means either the release for their version was never cut, or
  they installed from a commit whose `DESCRIPTION` version doesn't have a
  matching release yet.
