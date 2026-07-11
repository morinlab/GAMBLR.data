# compare_bedpe.R
#
# Compares bedpe (Manta SV) rows between the OLD bundled sample_data.rda and
# the NEW gambl_mutations.db build, to determine whether a row-count
# difference comes from more samples having SVs at all (healthy growth) or
# the same samples yielding more SVs each (possible duplication / behavior
# change) -- e.g. to diagnose why bedpe's regression check in
# test_gambl_db.R doesn't match the old snapshot as closely as maf/ashm/seg.
#
# Usage: Rscript data-raw/compare_bedpe.R [gambl_mutations.db] [old sample_data.rda]

suppressMessages({library(DBI); library(dplyr)})

args    <- commandArgs(trailingOnly = TRUE)
db      <- if (length(args) >= 1) args[[1]] else "gambl_mutations.db"
old_rda <- if (length(args) >= 2) args[[2]] else "data/sample_data.rda"

# --- new build (from the DB just written) -----------------------------------
con <- DBI::dbConnect(RSQLite::SQLite(), db)
new_bedpe <- DBI::dbGetQuery(con, "SELECT tumour_sample_id, genome_build FROM bedpe")
DBI::dbDisconnect(con)

# --- old build (from the bundled .rda used by the optional regression check) -
e <- new.env()
load(old_rda, envir = e)
old_sd <- get("sample_data", envir = e)
old_bedpe <- bind_rows(
  old_sd$grch37$bedpe %>% mutate(genome_build = "grch37"),
  old_sd$hg38$bedpe   %>% mutate(genome_build = "hg38")
)

cat("=== row counts ===\n")
cat("old:", nrow(old_bedpe), " new:", nrow(new_bedpe), "\n\n")

old_samples <- unique(old_bedpe$tumour_sample_id)
new_samples <- unique(new_bedpe$tumour_sample_id)

cat("=== unique samples with >=1 SV ===\n")
cat("old:", length(old_samples), " new:", length(new_samples), "\n\n")

cat("=== avg rows per sample ===\n")
cat("old:", round(nrow(old_bedpe) / length(old_samples), 2),
    " new:", round(nrow(new_bedpe) / length(new_samples), 2), "\n\n")

cat("=== sample-set overlap ===\n")
cat("in BOTH:            ", length(intersect(old_samples, new_samples)), "\n")
cat("ONLY in new (added): ", length(setdiff(new_samples, old_samples)), "\n")
cat("ONLY in old (dropped):", length(setdiff(old_samples, new_samples)), "\n\n")

# For samples present in BOTH builds: did their row count change?
shared <- intersect(old_samples, new_samples)
old_per_sample <- old_bedpe %>% filter(tumour_sample_id %in% shared) %>%
  count(tumour_sample_id, name = "n_old")
new_per_sample <- new_bedpe %>% filter(tumour_sample_id %in% shared) %>%
  count(tumour_sample_id, name = "n_new")
cmp <- full_join(old_per_sample, new_per_sample, by = "tumour_sample_id") %>%
  mutate(n_old = coalesce(n_old, 0), n_new = coalesce(n_new, 0), delta = n_new - n_old)

cat("=== for the", length(shared), "shared samples, did per-sample SV count change? ===\n")
cat("more SVs in new build:", sum(cmp$delta > 0), "\n")
cat("fewer SVs in new build:", sum(cmp$delta < 0), "\n")
cat("unchanged:             ", sum(cmp$delta == 0), "\n")
cat("extra rows from shared samples having MORE (not from new samples):",
    sum(pmax(cmp$delta, 0)), "\n\n")

cat("=== top 10 shared samples by increase ===\n")
print(cmp %>% arrange(desc(delta)) %>% head(10))

# Bottom line: rows added by (new_samples - old_samples) is expected growth.
# Rows added by shared samples having MORE than before is the number that
# would need a real explanation -- e.g. Manta re-run, annotate_sv() filter
# change, duplication.
new_rows_from_new_samples <- nrow(filter(new_bedpe, tumour_sample_id %in% setdiff(new_samples, old_samples)))
dropped_rows <- nrow(filter(old_bedpe, tumour_sample_id %in% setdiff(old_samples, new_samples)))
shared_delta <- sum(cmp$delta)  # net change among samples present in both

cat("\n=== attribution of the", nrow(new_bedpe) - nrow(old_bedpe), "row increase ===\n")
cat("from newly-added samples:         +", new_rows_from_new_samples, "\n")
cat("from dropped samples (removed):   -", dropped_rows, "\n")
cat("net change among shared samples:  ", ifelse(shared_delta >= 0, "+", ""), shared_delta, "\n")
cat("sum (should equal the row increase above):",
    new_rows_from_new_samples - dropped_rows + shared_delta, "\n")
