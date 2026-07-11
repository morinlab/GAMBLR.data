#' Lightweight sample metadata for the bundled GAMBL samples.
#'
#' The per-sample metadata table (formerly `sample_data$meta`), split into its
#' own object so that loading metadata does not pull the large MAF/seg/bedpe
#' frames of `sample_data` into memory. The heavy mutation frames now live in
#' the separate `gambl_mutations.db` SQLite database (see
#' `data-raw/build_mutations_db.R`). Rebuilt by `data-raw/make_sample_metadata.R`.
#'
#' @format A data frame with 3344 rows and 31 columns, keyed by `sample_id`.
"sample_metadata"
