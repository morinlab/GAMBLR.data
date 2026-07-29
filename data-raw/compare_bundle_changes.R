# compare_bundle_changes.R
#
# Compares the OLD bundled sample_data.rda against a NEW gambl_mutations.db
# build at the level of individual events (not just row totals), for each of
# three data types: SNV (maf), CNV (seg), SV (bedpe). For each type:
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
#
# Usage: Rscript data-raw/compare_bundle_changes.R [gambl_mutations.db] [old sample_data.rda] [output dir]

suppressMessages({library(DBI); library(dplyr)})

args    <- commandArgs(trailingOnly = TRUE)
db      <- if (length(args) >= 1) args[[1]] else "gambl_mutations.db"
old_rda <- if (length(args) >= 2) args[[2]] else "data/sample_data.rda"
outdir  <- if (length(args) >= 3) args[[3]] else "."
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

con <- DBI::dbConnect(RSQLite::SQLite(), db)
on.exit(DBI::dbDisconnect(con), add = TRUE)

e <- new.env()
load(old_rda, envir = e)
old_sd <- get("sample_data", envir = e)

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

# --- SNV: maf (which now also absorbs what used to be a separate ashm pull,
# deduplicated against it at write time -- see write_mutations_db.R) --------
# .source records which table/pull each row came from on the OLD side so
# gained/lost examples can still be traced back to the right code path.
# The new side no longer has that distinction to make (maf and ashm were
# merged into one table with no retained provenance column -- deliberately;
# nothing consumed it), so .source is uniformly "maf" there. Not part of the
# comparison key either way: a row that moved from one source to the other
# between builds still counts as retained if its sample/position match.
new_snv <- dbGetQuery(con, "SELECT * FROM maf") %>% mutate(.source = "maf")
old_snv <- bind_rows(
  old_sd$grch37$maf  %>% mutate(genome_build = "grch37", .source = "maf"),
  old_sd$hg38$maf    %>% mutate(genome_build = "hg38", .source = "maf"),
  old_sd$grch37$ashm %>% mutate(genome_build = "grch37", .source = "ashm"),
  old_sd$hg38$ashm   %>% mutate(genome_build = "hg38", .source = "ashm")
)

# Pipeline == "strelka" (Arthur's old, unfiltered raw-flat-file dump) is
# intentionally and entirely gone from the new bundle -- that source file
# isn't read anymore. Confirmed (twice: 08-15460, DO52686) that its
# Tumor_Sample_Barcode values are the raw file's own non-standard IDs, not
# GAMBL's real sample_id, so every one of these rows is expected to show up
# as "lost" regardless of anything else in the pipeline. Excluded here so
# that expected, already-understood disappearance doesn't drown out
# genuinely new findings on every run.
old_snv <- old_snv %>% filter(is.na(Pipeline) | Pipeline != "strelka")

snv_key <- c("Tumor_Sample_Barcode", "genome_build", "Chromosome", "Start_Position", "End_Position")
snv_cmp <- compare_events(old_snv, new_snv, "Tumor_Sample_Barcode", snv_key, "SNV (maf, incl. aSHM-region pull)")
write_examples(snv_cmp$gained, "Tumor_Sample_Barcode", file.path(outdir, "snv_gained_examples.log"))
write_examples(snv_cmp$lost, "Tumor_Sample_Barcode", file.path(outdir, "snv_lost_examples.log"))
write_counts(snv_cmp$per_sample, file.path(outdir, "snv_persample_counts.tsv"))

# --- CNV: seg -----------------------------------------------------------------
new_cnv <- dbGetQuery(con, "SELECT * FROM seg")
old_cnv <- bind_rows(
  old_sd$grch37$seg %>% mutate(genome_build = "grch37"),
  old_sd$hg38$seg   %>% mutate(genome_build = "hg38")
)
cnv_key <- c("ID", "genome_build", "chrom", "start", "end")
cnv_cmp <- compare_events(old_cnv, new_cnv, "ID", cnv_key, "CNV (seg)")
write_examples(cnv_cmp$gained, "ID", file.path(outdir, "cnv_gained_examples.log"))
write_examples(cnv_cmp$lost, "ID", file.path(outdir, "cnv_lost_examples.log"))
write_counts(cnv_cmp$per_sample, file.path(outdir, "cnv_persample_counts.tsv"))

# --- SV: bedpe ------------------------------------------------------------------
new_sv <- dbGetQuery(con, "SELECT * FROM bedpe")
old_sv <- bind_rows(
  old_sd$grch37$bedpe %>% mutate(genome_build = "grch37"),
  old_sd$hg38$bedpe   %>% mutate(genome_build = "hg38")
)
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
