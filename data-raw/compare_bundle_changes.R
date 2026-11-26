# compare_bundle_changes.R
#
# Compares a PREVIOUS gambl_mutations.db (e.g. the current release, as cached by
# GAMBLR.data::gambl_mutations_db()) against a NEW build at the level of
# individual events (not just row totals), for each of three data types: SNV
# (maf), CNV (seg), SV (bedpe). For each type:
#   - counts, per sample_id, how many events are NEW (gained: present in the
#     new build but not the old) and how many are LOST (missing: present in
#     the old build but not the new)
#   - writes one example gained-event row per affected sample to
#     <outdir>/{snv,cnv,sv}_gained_examples.log, and one example lost-event
#     row per affected sample to <outdir>/{snv,cnv,sv}_lost_examples.log --
#     for quick manual inspection without dumping every gained/lost row. At
#     most one row per sample that has >=1 gained (or lost) event, so at most
#     as many rows as there are samples -- in practice fewer, since not every
#     sample is affected.
#   - also writes the full per-sample gained/lost counts to
#     <outdir>/{snv,cnv,sv}_persample_counts.tsv for further analysis.
#
# "Same event" is determined by an exact match on a type-specific key (sample
# id + genome_build + genomic coordinates), NOT by row position or count --
# two rows only count as "the same" if every key column matches exactly.
# Chromosome names are compared without a "chr" prefix, so a build that only
# changed naming (e.g. grch37 rows once stored as chr1) doesn't show as
# lost/gained.
#
# (Until the sqlite database replaced it, the old side was the bundled
# sample_data.rda object; that object has been retired.)
#
# Usage: Rscript data-raw/compare_bundle_changes.R [new gambl_mutations.db] [previous gambl_mutations.db] [output dir]
#   previous defaults (when missing or "") to the cached release, tools::R_user_dir("GAMBLR.data", "cache")/gambl_mutations.db

suppressMessages({library(DBI); library(dplyr)})

args    <- commandArgs(trailingOnly = TRUE)
db      <- if (length(args) >= 1) args[[1]] else "gambl_mutations.db"
old_db  <- if (length(args) >= 2 && nzchar(args[[2]])) args[[2]] else file.path(tools::R_user_dir("GAMBLR.data", "cache"), "gambl_mutations.db")
outdir  <- if (length(args) >= 3) args[[3]] else "."
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(db), file.exists(old_db))
if (normalizePath(db) == normalizePath(old_db)) stop("the new and previous databases are the same file")
con <- DBI::dbConnect(RSQLite::SQLite(), db, flags = RSQLite::SQLITE_RO)
on.exit(DBI::dbDisconnect(con), add = TRUE)
old_con <- DBI::dbConnect(RSQLite::SQLite(), old_db, flags = RSQLite::SQLITE_RO)
on.exit(DBI::dbDisconnect(old_con), add = TRUE)
message("new: ", db, "\nprevious: ", old_db)
unprefix <- function(df, cols) { for (cc in cols) df[[cc]] <- sub("^chr", "", as.character(df[[cc]])); df }

# Generic per-sample gained/lost comparison for one data type.
#   old_df, new_df : data frames (already build-stamped with genome_build if
#                    the type spans both builds)
#   sample_col     : name of the sample-id column for this type
#   key_cols       : columns (including sample_col and genome_build, if
#                    present) that together define a unique "event" -- an
#                    exact match on all of these means "the same event"
compare_events <- function(old_df, new_df, sample_col, key_cols, label) {
  key_of <- function(df) {
    # ":::" separator avoids ambiguous concatenation, e.g. chrom="1",
    # start=23456789 vs chrom="12",start=3456789 both -> "123456789" if
    # pasted with no separator at all
    do.call(paste, c(as.list(select(df, all_of(key_cols))), sep = ":::"))
  }
  old_df$.key <- key_of(old_df)
  new_df$.key <- key_of(new_df)

  gained   <- new_df %>% filter(!.key %in% old_df$.key)
  lost     <- old_df %>% filter(!.key %in% new_df$.key)
  retained <- old_df %>% filter(.key %in% new_df$.key)

  gained_per_sample   <- gained   %>% count(.data[[sample_col]], name = "n_gained")
  lost_per_sample     <- lost     %>% count(.data[[sample_col]], name = "n_lost")
  retained_per_sample <- retained %>% count(.data[[sample_col]], name = "n_retained")

  # n_retained (present in both old and new) is joined in separately -- it's
  # looked up only for the samples already in per_sample (those with >=1
  # gained or lost event), not full_join'd, since retained_per_sample on its
  # own covers essentially every sample including untouched ones.
  per_sample <- full_join(gained_per_sample, lost_per_sample, by = sample_col) %>%
    left_join(retained_per_sample, by = sample_col) %>%
    mutate(
      n_gained   = coalesce(n_gained, 0L),
      n_lost     = coalesce(n_lost, 0L),
      n_retained = coalesce(n_retained, 0L)
    ) %>%
    arrange(desc(n_gained + n_lost))

  cat(sprintf(
    "\n=== %s: %d gained rows, %d lost rows (%d samples w/ gains, %d samples w/ losses) ===\n",
    label, nrow(gained), nrow(lost),
    sum(per_sample$n_gained > 0), sum(per_sample$n_lost > 0)
  ))
  print(head(per_sample, 10))

  list(per_sample = per_sample, gained = gained, lost = lost)
}

# One example row per affected sample, from either the gained or the lost
# data frame returned by compare_events() -- whichever is passed in.
write_examples <- function(event_df, sample_col, path) {
  if (nrow(event_df) == 0) {
    message("no rows -- skipping ", path)
    return(invisible(NULL))
  }
  examples <- event_df %>%
    group_by(.data[[sample_col]]) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    select(-.key)
  write.table(examples, path, sep = "\t", quote = FALSE, row.names = FALSE)
  message("wrote ", nrow(examples), " example rows to ", path)
}

write_counts <- function(per_sample, path) {
  write.table(per_sample, path, sep = "\t", quote = FALSE, row.names = FALSE)
  message("wrote per-sample counts to ", path)
}

# --- SNV: maf (aSHM-region mutations included; maf and ashm are one table) ---
new_snv <- dbGetQuery(con, "SELECT * FROM maf") %>% unprefix("Chromosome")
old_snv <- dbGetQuery(old_con, "SELECT * FROM maf") %>% unprefix("Chromosome")

snv_key <- c("Tumor_Sample_Barcode", "genome_build", "Chromosome", "Start_Position", "End_Position")
snv_cmp <- compare_events(old_snv, new_snv, "Tumor_Sample_Barcode", snv_key, "SNV (maf, incl. aSHM regions)")
write_examples(snv_cmp$gained, "Tumor_Sample_Barcode", file.path(outdir, "snv_gained_examples.log"))
write_examples(snv_cmp$lost, "Tumor_Sample_Barcode", file.path(outdir, "snv_lost_examples.log"))
write_counts(snv_cmp$per_sample, file.path(outdir, "snv_persample_counts.tsv"))

# --- CNV: seg -----------------------------------------------------------------
new_cnv <- dbGetQuery(con, "SELECT * FROM seg") %>% unprefix("chrom")
old_cnv <- dbGetQuery(old_con, "SELECT * FROM seg") %>% unprefix("chrom")
cnv_key <- c("ID", "genome_build", "chrom", "start", "end")
cnv_cmp <- compare_events(old_cnv, new_cnv, "ID", cnv_key, "CNV (seg)")
write_examples(cnv_cmp$gained, "ID", file.path(outdir, "cnv_gained_examples.log"))
write_examples(cnv_cmp$lost, "ID", file.path(outdir, "cnv_lost_examples.log"))
write_counts(cnv_cmp$per_sample, file.path(outdir, "cnv_persample_counts.tsv"))

# --- SV: bedpe ------------------------------------------------------------------
new_sv <- dbGetQuery(con, "SELECT * FROM bedpe") %>% unprefix(c("CHROM_A", "CHROM_B"))
old_sv <- dbGetQuery(old_con, "SELECT * FROM bedpe") %>% unprefix(c("CHROM_A", "CHROM_B"))
sv_key <- c("tumour_sample_id", "genome_build", "CHROM_A", "START_A", "CHROM_B", "START_B")
sv_cmp <- compare_events(old_sv, new_sv, "tumour_sample_id", sv_key, "SV (bedpe)")
write_examples(sv_cmp$gained, "tumour_sample_id", file.path(outdir, "sv_gained_examples.log"))
write_examples(sv_cmp$lost, "tumour_sample_id", file.path(outdir, "sv_lost_examples.log"))
write_counts(sv_cmp$per_sample, file.path(outdir, "sv_persample_counts.tsv"))

cat("\n=== summary ===\n")
cat(sprintf("SNV: %d samples with gains, %d samples with losses\n",
            sum(snv_cmp$per_sample$n_gained > 0), sum(snv_cmp$per_sample$n_lost > 0)))
cat(sprintf("CNV: %d samples with gains, %d samples with losses\n",
            sum(cnv_cmp$per_sample$n_gained > 0), sum(cnv_cmp$per_sample$n_lost > 0)))
cat(sprintf("SV:  %d samples with gains, %d samples with losses\n",
            sum(sv_cmp$per_sample$n_gained > 0), sum(sv_cmp$per_sample$n_lost > 0)))
