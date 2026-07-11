# write_mutations_db.R
#
# Reusable build-time helper: given the in-memory `sample_data` list (as
# assembled by data-raw/assemble_bundled_data.R), write the normalized
# `gambl_mutations.db` SQLite database with indexes matching the GAMBLR.open
# accessors. This is the ONGOING path that replaces bundling `sample_data.rda`:
#   assemble_bundled_data.R  ->  sample_data (in memory)  ->  write_mutations_db()
#
# `source()` this file (it defines the function; it does not run anything).
# Requires the DBI and RSQLite packages (Imports of GAMBLR.data).

write_mutations_db <- function(sample_data,
                               out_db = "gambl_mutations.db",
                               builds = c("grch37", "hg38"),
                               source_desc = "assemble_bundled_data.R") {
  stopifnot(is.list(sample_data), !is.null(sample_data$meta))

  if (file.exists(out_db)) file.remove(out_db)
  con <- DBI::dbConnect(RSQLite::SQLite(), out_db)
  on.exit(DBI::dbDisconnect(con))
  # speed up the bulk load
  DBI::dbExecute(con, "PRAGMA journal_mode = OFF")
  DBI::dbExecute(con, "PRAGMA synchronous = OFF")

  # sample metadata (single table, no genome build)
  DBI::dbWriteTable(con, "sample_meta", as.data.frame(sample_data$meta),
                    overwrite = TRUE)

  # genome-build-stamped frames; append per build to keep peak memory low
  write_element <- function(elem) {
    first <- TRUE
    total <- 0L
    for (b in builds) {
      x <- sample_data[[b]][[elem]]
      if (is.null(x) || !nrow(x)) next
      x <- as.data.frame(x)
      x$genome_build <- b
      DBI::dbWriteTable(con, elem, x, append = !first, overwrite = first)
      total <- total + nrow(x)
      first <- FALSE
    }
    total
  }
  counts <- vapply(c("maf", "ashm", "seg", "bedpe"), write_element, integer(1))

  # indexes mirroring how the GAMBLR.open accessors query the data.
  # wrapped in try() so a build missing an optional table/column is non-fatal.
  idx <- c(
    "CREATE INDEX idx_maf_pos     ON maf(genome_build, Chromosome, Start_Position)",
    "CREATE INDEX idx_maf_sample  ON maf(Tumor_Sample_Barcode)",
    "CREATE INDEX idx_maf_gene    ON maf(Hugo_Symbol)",
    "CREATE INDEX idx_maf_pipe    ON maf(Pipeline)",
    "CREATE INDEX idx_maf_study   ON maf(Study)",
    "CREATE INDEX idx_ashm_pos    ON ashm(genome_build, Chromosome, Start_Position)",
    "CREATE INDEX idx_ashm_sample ON ashm(Tumor_Sample_Barcode)",
    "CREATE INDEX idx_seg_sample  ON seg(genome_build, ID)",
    "CREATE INDEX idx_seg_pos     ON seg(genome_build, chrom, start)",
    "CREATE INDEX idx_bedpe_sample ON bedpe(tumour_sample_id)",
    "CREATE INDEX idx_meta_sample  ON sample_meta(sample_id)",
    "CREATE INDEX idx_meta_barcode ON sample_meta(Tumor_Sample_Barcode)"
  )
  for (stmt in idx) try(DBI::dbExecute(con, stmt), silent = TRUE)

  # self-describing provenance / expected-counts table (used by the tests)
  DBI::dbWriteTable(con, "build_info", data.frame(
    key = c("source", "built_at", "builder",
            "n_maf", "n_ashm", "n_seg", "n_bedpe", "n_samples"),
    value = c(source_desc,
              format(Sys.time(), tz = "UTC", usetz = TRUE),
              "write_mutations_db",
              counts[["maf"]], counts[["ashm"]], counts[["seg"]],
              counts[["bedpe"]], nrow(sample_data$meta)),
    stringsAsFactors = FALSE
  ), overwrite = TRUE)

  DBI::dbExecute(con, "VACUUM")
  DBI::dbExecute(con, "ANALYZE")

  message(sprintf("Wrote %s (%.0f MB)  maf=%d ashm=%d seg=%d bedpe=%d samples=%d",
                  out_db, file.info(out_db)$size / 1024^2,
                  counts[["maf"]], counts[["ashm"]], counts[["seg"]],
                  counts[["bedpe"]], nrow(sample_data$meta)))
  invisible(out_db)
}
