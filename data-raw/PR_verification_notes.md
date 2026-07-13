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
