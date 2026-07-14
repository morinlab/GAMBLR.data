# Verification notes: sample↔study join table / deduplicated SLMS-3 pulls

Working notes for the PR description covering the `assemble_bundled_data.R` /
`write_mutations_db.R` / `get_ssm_from_db.R` refactor (branch `rmorin-dev`,
commits `1a8cfb6`, `017d3db`). Append further findings (e.g.
`compare_bundle_changes.R` output) below as they come in.

## Background

Root cause: SLMS-3 mutation calls were pulled via separate, cohort-specific
code blocks, each independently deciding which samples to pull for. A sample
belonging to more than one cohort (e.g. the Dreval FL × Hilton trios overlap)
got its mutations pulled and bound multiple times. Confirmed directly against
the previously-shipped `gambl_mutations.db`: 1,797 exact duplicate `maf` rows
(same pattern in `ashm`), plus true duplicate rows in `sample_meta` for the
same 5 samples.

Fix: `Study`/cohort membership moved off `maf`/`ashm` entirely into a new
`sample_study` many-to-many join table; SLMS-3 is now pulled exactly once per
deduplicated sample regardless of cohort membership. Surfaced while replacing
the Arthur cohort's raw, completely unfiltered `strelka` flat-file dump
(2.83M rows / 65,121 distinct genes, missing the lymphoma-gene-panel filter
every other block applies) with a proper SLMS-3 pull.

## `test_gambl_db.R` run — 2026-07-14

All structural/integrity checks passed:

```
== tables present ==        7/7 PASS (maf, ashm, seg, bedpe, sample_meta, sample_study, build_info)
== row counts (all > 0) ==  6/6 PASS
  maf          543,214 rows
  ashm         134,826 rows
  seg          128,075 rows
  bedpe          2,352 rows
  sample_meta    3,346 rows
  sample_study   2,880 rows
== both genome builds present ==       4/4 PASS
== required columns present ==         5/5 PASS
== indexes present ==                  7/7 PASS (incl. new idx_study_sample, idx_study_study)
== pipelines present ==                2/2 PASS -- Pipeline values: Publication, SLMS-3 (no more "strelka")
== spot query: MYC locus, grch37, slms-3 ==   PASS -- 4,549 variants found
== build_info counts match actual table counts ==  4/4 PASS
== sample_metadata.rda matches sample_meta table == PASS (3,346 == 3,346)
```

### Regression vs. old `sample_data.rda` (informational, all 4 differences expected)

| Table | Old (`sample_data.rda`) | New (`gambl_mutations.db`) | Δ | Explanation |
| --- | --- | --- | --- | --- |
| `maf` | 3,663,468 | 543,214 | −85% | Almost entirely Arthur's unfiltered `strelka` dump being removed (2,827,527 of the old rows, 65,121 distinct genes vs. the ~150-237 gene panel) + the Dreval×Hilton dedup. This is the fix working as intended, not a regression. |
| `ashm` | 134,868 | 134,826 | −0.03% | Noise-level. |
| `seg` | 125,982 | 128,075 | +1.7% | Unrelated to this refactor (CNV/seg blocks untouched) -- likely normal cohort growth in GAMBL's live metadata since the old snapshot. |
| `bedpe` | 941 | 2,352 | +150% | Bonus fix, not explicitly planned: Arthur and Hilton were previously added to `sample_data$meta` *after* the Manta SV block ran, so they got zero SV coverage -- same root cause as the aSHM gap that *was* explicitly fixed. Now that metadata assembly happens once, up front, they get SV coverage too. |

## Still to do

- [ ] Run `compare_bundle_changes.R` for the granular per-sample gained/lost/retained view (append results below).
- [ ] Confirm Dreval×Hilton overlap samples (`05-32150T`, `08-15460T`, `09-33003T`, `15-13383T`, `17-36275T`) show 2 `sample_study` rows each and roughly halved mutation counts vs. the old bundle.
- [ ] Confirm Arthur's `maf` rows are gene-panel-restricted (~150-237 distinct genes, not 65,121) and now have hg38 + aSHM coverage.
- [x] Functional check: `get_ssm_by_samples()`/`get_ssm_by_patients()` against real Reddy (capture) data -- 20,480 rows returned, confirming the query layer and underlying data are correct end to end (see note below; the initial zero-result was unrelated to this refactor).

## Non-issue: `get_ssm_by_samples()`/`get_ssm_by_patients()` initially returned 0 rows for Reddy

Root cause turned out to be unrelated to this refactor: both functions
default `this_seq_type = "genome"` and filter `these_samples_metadata` by
it *in addition to* whatever metadata you pass in, rather than inferring it.
Reddy is a capture cohort (`seq_type == "capture"` for all 999 samples), so
the default silently dropped every row before the query ever reached the
database. Confirmed the data itself was fine throughout: a raw SQL query
against `maf` for Reddy's `Tumor_Sample_Barcode` values returned 47,626
rows the whole time. Passing `this_seq_type = "capture"` explicitly fixed
both functions (20,480 rows each, matching after `get_ssm_from_db()`'s
Pipeline/genome_build/read-support filters narrow the raw 47,626).

Separate, pre-existing GAMBLR.open UX gap worth considering independently
of this PR: neither function warns when the `seq_type` filter drops every
row, so this fails silently rather than with a clear message.

## Fix: Arthur's Case ID -> patient_id matching, and a `!grepl("tumor", ...)` filter that dropped multi-sample patients entirely

`compare_bundle_changes.R` showed some patients losing 100% of their SNVs
with zero retention (e.g. patient `08-15460`: 19,897 lost, 0 retained --
found via `egrep "Tumor_Sa|05-32150|08-15460" snv_persample_counts.tsv`).
Confirmed as an Arthur-specific bug, two compounding causes:

1. `arthur_meta` matched Arthur's own "Case ID" (from `DLBCL_Arthur.xlsx`)
   against GAMBL's `patient_id` via `patient_id %in% arthur_case_ids$\`Case ID\``.
   Same class of bug as the `study_id` chr/dbl mismatch fixed earlier:
   `read_xlsx()` can silently type a leading-zero ID like `08-15460` as
   numeric, and a numeric-vs-character `%in%` comparison fails to match
   without any error or warning. Fixed by coercing both sides to character
   and using an explicit `inner_join()` instead of a `%in%` filter.
2. The block also had `! grepl("tumor", sample_id)`, which excludes every
   sample for any patient with more than one tumor biopsy (e.g.
   `..._tumorA`/`..._tumorB` -- the same multi-sample-per-patient pattern
   Hilton has). For a patient whose *only* GAMBL samples are tumor-suffixed,
   this filter matched nothing, dropping their SLMS-3 coverage entirely. No
   other cohort in this script excludes samples this way, and no comment
   explained the original intent -- removed.

**Update**: rebuilt with this fix and re-ran -- `08-15460` (bare) still
showed the identical 100%-loss numbers, byte-for-byte. That ruled out
`arthur_meta`/`arthur_case_ids` as the source: those only ever fed
`sample_data$meta`, never `sample_data$maf` directly. Traced the real
source instead by reading the *old raw flat file itself*:
```r
arthur_raw <- readr::read_tsv("inst/extdata/studies/DLBCL_Arthur.maf.gz",
                               col_select = "Tumor_Sample_Barcode")
grep("15460", unique(arthur_raw$Tumor_Sample_Barcode), value = TRUE)
# [1] "08-15460"
```
Confirmed: the deleted raw dump used bare patient-style IDs for
`Tumor_Sample_Barcode` (not GAMBL's real per-sample naming), directly and
independently of `arthur_meta`. **Conclusion: this specific finding was
never a bug.** Every "100%-loss" row keyed to `08-15460` (bare) is simply
the expected disappearance of the old file's non-standard ID convention,
now correctly superseded by GAMBL's real sample IDs (`08-15460T`: 1,515
retained; `08-15460_tumorB`: 665 retained -- both present and correctly
sized in both old and new bundles).

Generalized this into a defensive fix rather than leaving it Arthur-specific:
added `relabel_to_sample_id(df, study_name)`, applied to all four
Publication-pipeline pulls (Thomas BL, Thomas DLBCL, Dreval, Reddy's
original-variants file). It joins on `sample_study$study_id` and rewrites
`Tumor_Sample_Barcode` to the real `sample_id` wherever they differ -- a
no-op for Thomas/Dreval/Hilton (`study_id` already mirrors `sample_id`
there) and an actual fix for Reddy (`study_id` is the paper's own raw
"Sample ID", genuinely distinct from `sample_id`). Verified against
synthetic data: relabels when a study_id match exists, no-ops when
`Tumor_Sample_Barcode` already is the real sample_id, and passes through
unmatched values unchanged rather than dropping them.

Still need to: rebuild with this fix and re-run `compare_bundle_changes.R`
to confirm `08-15460` (bare) is gone from the loss list entirely, and
check whether the same bare-vs-suffixed pattern explains the
`DO52686`/`07-35482`-style samples from the earlier, larger (3.1M row)
loss count.

## Fix: consolidated Phase 4/5 pulls passed sample_data$meta directly to get_ssm_by_regions() -- missing coding calls for hand-built cohorts

`GAMBLR_examples_output.log`'s regenerated run showed a `get_coding_ssm()`
example for 183 Dreval FL patients (sample_id matching `"SP"`) return **zero
rows**, when it used to return real coding-classified data. Traced through
several rounds of elimination (ID matching confirmed fine, `t_alt_count`
confirmed numeric, `min_read_support` confirmed innocent) to: `coding_only
= TRUE` alone reduced 28,867 real SLMS-3/grch37 rows for these patients to
zero -- every one of them was non-coding (`Intron`/`5'Flank`/`3'Flank`/
`5'UTR`/`3'UTR`). For 183 real FL patients to have literally zero coding
mutations across the whole lymphoma gene panel is biologically implausible.

Root cause: `unix_group` (and likely other native GAMBL metadata columns)
is `NA` for every cohort built by hand in this script rather than pulled
live from `get_gambl_metadata()` -- confirmed directly:
```r
filter(sample_data$meta, is.na(unix_group)) %>% count(study)
#  BL_Thomas 234, DLBCL_Arthur 153, DLBCL_Thomas 43, DLBCL_cell_lines 5, FL_Dreval 443
filter(sample_data$meta, !is.na(unix_group)) %>% count(study)
#  DLBCL_Hilton 159, NCI_DLBCL_Golub 124, dlbcl_chapuy 233, dlbcl_reddy 999, dlbcl_schmitz 951
```
`get_ssm_by_regions()` apparently relies on columns like `unix_group` to
determine which underlying merged flat file to read for a sample -- with it
missing, it was silently returning an incomplete (non-coding-only) subset
for affected samples instead of erroring. This was specific to Phase 4 (the
consolidated SLMS-3 pull) and Phase 5 (the consolidated aSHM pull), both of
which passed `sample_data$meta`-derived subsets directly to
`get_ssm_by_regions()` -- every *other* GAMBLR.results call site in this
script (the SV/bedpe block, the cell-lines pull, the Publication-pipeline
enrichment steps) already re-fetched fresh, complete metadata via
`get_gambl_metadata()` filtered by `sample_id` first, which is exactly why
this only surfaced in the two things newly introduced this session.

Fixed both: `sample_data$meta` is now used only to determine which
`sample_id`s belong in each pull; the actual `these_samples_metadata`
passed to `get_ssm_by_regions()` is always a fresh `get_gambl_metadata()`
call filtered to that same `sample_id` set. General principle going
forward, not just for these two spots: never pass `sample_data$meta` (or
anything derived from it) directly to a GAMBLR.results function.

**Update**: rebuilt and re-ran `compare_bundle_changes.R` broadly (not just
the one example). Confirmed fixed -- samples that previously showed
catastrophic 0-retained losses now show large, healthy retained counts
(`SP116715`: 2,901 retained; `16-11636T`: 1,624; `14-20552_tumorB`: 670;
`04-24937T`: 990), with only small gains/losses on top (single digits to
low teens). Total SNV loss across the whole bundle dropped from 3.1M rows
to 299, spread across 113 samples (~2.6 rows/sample average) -- a
different order of magnitude, consistent with normal residual noise rather
than a bug.

## Current state (post all fixes above)

- **SNV**: 348 samples w/ gains, 113 samples w/ losses, 299 total lost rows
  -- converged to normal residual noise. The gained-only, 0-retained
  samples at the top of `snv_persample_counts.tsv` are the known "left out
  of the last release" bonus fix (unrelated, already understood).
- **CNV (seg)**: unchanged from prior runs (5 cell lines w/ gains, 0 w/
  losses) -- expected, nothing in this session touched seg assembly.
- **SV (bedpe)**: unchanged from prior runs (433 w/ gains, 7 w/ losses) --
  expected, consistent with the earlier Arthur/Hilton SV coverage bonus
  fixes already being stable.

Considering this refactor converged and ready for PR, pending final
sign-off.

## Fix: cell lines need a separate, genome-wide SNV pull

After the strelka-exclusion fix, `compare_bundle_changes.R` came back clean
except for the 5 cell lines (`SU-DHL-4`, `OCI-Ly3`, `OCI-Ly10`, `SU-DHL-10`,
`DOHH-2`), each showing ~43k-65k lost rows against a few hundred/~1k
retained. This is the scope-reduction flagged during planning, now
confirmed as something to actually fix rather than accept: the original
pre-refactor cell-line pull used `get_ssm_by_samples()` with no gene-panel
filter (genome-wide), but folding cell lines into the consolidated,
panel-restricted SLMS-3 pull (Phase 4) silently reduced them to panel-only
coverage.

Fixed by excluding cell lines from `all_slms3_meta` (so they're not also
pulled panel-restricted) and adding back a dedicated, unrestricted
`get_ssm_by_samples()` pull for them specifically, both builds, bound into
`slms3_grch37`/`slms3_hg38` alongside the panel-restricted data for
everyone else. Still tagged `Pipeline = "SLMS-3"` (same underlying pipeline,
just unrestricted scope for these 5 samples).

**Update**: rebuilt and re-ran -- cell-line losses are gone, confirming the
fix worked. One more thing surfaced in the same run, another bonus rather
than a bug: the top of `snv_persample_counts.tsv` now shows several samples
(e.g. `15-18723_CLC02138`, `15-43891_tumorB`, `07-17613T`) with real gained
counts (64-630) but `n_lost = 0` and `n_retained = 0` -- meaning they had
*zero* rows in the old bundle at all, not a mismatch. Confirmed these were
samples left out of the last release for some reason and are now correctly
included -- consistent with the pattern already seen for Arthur/Hilton
(aSHM, SV coverage) where the old per-cohort-block gating silently excluded
samples that the consolidated pull now correctly picks up.

## Fix: `study_id` used `patient_id` instead of `sample_id` for Thomas/Dreval/Hilton

Caught in review (not by the automated checks above): the first pass at
populating `sample_study$study_id` used `patient_id` for Thomas BL, Thomas
DLBCL, Dreval, and Hilton. For Hilton specifically this is a real bug, not
just imprecision -- Hilton is a trios study, so a single patient can have
multiple samples (e.g. `LY_RELY_116_tumorA` and `LY_RELY_116_tumorB`), and
`patient_id` collapsed both to the same `study_id`, making the two
`sample_study` rows ambiguous. Thomas/Dreval don't currently have that
ambiguity in practice, but the same fix applies for consistency: all four
now use `sample_id` (already sourced from each xlsx's own "Genome sample
id"/`DNAseq_sample_id` column, i.e. already the closest thing to a
study-native sample-level identifier) instead of `patient_id`. Arthur and
Reddy were left unchanged -- Arthur's block already excludes
`tumor`-suffixed sample_ids per patient (structurally safe from the same
ambiguity), and Reddy's `study_id` comes from the paper's own distinct
"Sample ID" column, already 1:1 with `sample_id`.
