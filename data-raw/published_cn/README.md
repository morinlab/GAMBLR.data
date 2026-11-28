# Published copy number segments

Inputs for `data-raw/published_cn_segments.R`, which builds
`inst/extdata/published_cn_segments.tsv.gz` (read by
`get_published_cn_segments()`): segments from published analyses of cohorts
that are also in GAMBL, for comparing GAMBL's segmentations with them.

## Raw inputs (`raw/`)

| File | Source | SHA-256 (uncompressed) |
|---|---|---|
| `DLBCL_699_segs.2021-12-15.ccf.seg.gz` | DLBClass (Chapuy et al.), github.com/getzlab/DLBCL-Classifier, `data_tables/additional_gsm_inputs/`, commit 5394f092 (2025-09-24) | 38c421b49575f25cb51ea2b4ba8268031435939ccd42d8416c7ee591ec986d96 |
| `schmitz_cnv_fix.seg.gz` | Schmitz et al. (NEJM 2018) array copy number, as used for the LymphGen validation | 4e0599dcadd761dfa250013bd0116b9ba87b91515198a28342751e30c2de8a9c |

## Sample maps

- `dlbclass_sample_map.tsv`: DLBClass sample IDs (`source_id`) to GAMBL
  `sample_id`, for the 627 samples in GAMBL (214 Chapuy, 413 Schmitz), from the
  linking of DLBClass and GAMBL samples done for the LymphGen/DLBClass
  benchmarks. 72 DLBClass samples aren't in GAMBL.
- `schmitz_array_sample_map.tsv`: Schmitz array names to GAMBL `sample_id`
  (name + "T"), for the 460 of 476 arrays whose sample is in the bundled
  metadata as a dlbcl_schmitz capture sample.

The DLBClass segments of the Schmitz cohort are DLBClass's reprocessing of the
same arrays as `schmitz_array`: two processings of one array dataset.
