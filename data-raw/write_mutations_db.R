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

  # sample metadata (single table, no genome build). reference_PMID is
  # dropped here -- it's single-valued (one row per sample) and can't
  # correctly represent a sample belonging to more than one study (each with
  # its own PMID); sample_study is the source of truth for that now.
  sample_meta_out <- as.data.frame(sample_data$meta)
  sample_meta_out$reference_PMID <- NULL
  DBI::dbWriteTable(con, "sample_meta", sample_meta_out, overwrite = TRUE)

  # sample_study: many-to-many bridge table (sample_id, study) replacing the
  # old inline Study column on maf/ashm rows -- guarded so older sample_data
  # objects that predate this table (e.g. via the legacy .rda adapter) don't
  # hard-crash a build that simply won't have cohort-membership data.
  n_sample_study <- 0L
  if (!is.null(sample_data$sample_study) && nrow(sample_data$sample_study)) {
    DBI::dbWriteTable(con, "sample_study", as.data.frame(sample_data$sample_study),
                      overwrite = TRUE)
    n_sample_study <- nrow(sample_data$sample_study)
  }

  # Natural mutation-call key, reused from GAMBLR.open::get_ssm_by_region()
  # (which already applies the same distinct() at read time) as a write-time
  # safety net: no cohort-specific pull block should ever again be able to
  # introduce duplicate mutation rows for a sample that's pulled more than
  # once, regardless of the reason.
  dedup_keys <- list(
    maf  = c("Tumor_Sample_Barcode", "Chromosome", "Start_Position", "End_Position", "Tumor_Seq_Allele2"),
    ashm = c("Tumor_Sample_Barcode", "Chromosome", "Start_Position", "End_Position", "Tumor_Seq_Allele2")
  )

  # variant_pipeline: many-to-many bridge table (variant key, Pipeline),
  # same pattern as sample_study (sample_id, study). A single Pipeline column
  # on maf/ashm can only hold one value per variant, but the same real
  # mutation can legitimately be produced by more than one pipeline (e.g. a
  # cohort's own published/curated maf and our independent SLMS-3 recall
  # both calling the same position) -- when that happens the dedup below
  # collapses them to one row keyed on whichever pipeline came first in
  # sample_data[[b]][[elem]], silently making the call invisible to a
  # tool_name-filtered query for the other pipeline even though the variant
  # itself is still fully present in the table. Captured here, before dedup,
  # so no pipeline's claim to a variant is ever lost regardless of which row
  # "wins" as maf/ashm's single representative copy.
  variant_pipeline_rows <- list()

  # genome-build-stamped frames; append per build to keep peak memory low
  write_element <- function(elem) {
    first <- TRUE
    total <- 0L
    dk <- dedup_keys[[elem]]
    for (b in builds) {
      x <- sample_data[[b]][[elem]]
      if (is.null(x) || !nrow(x)) next
      x <- as.data.frame(x)
      x$genome_build <- b
      # Normalized to lowercase here, once, at the write boundary -- so
      # every reader (get_ssm_from_db(), raw SQL, this build's own
      # variant_pipeline capture below) can match Pipeline with a plain
      # equality against a plain index, instead of every query needing to
      # wrap the column in LOWER()/tolower() (which a normal b-tree index
      # can't be used to satisfy).
      if ("Pipeline" %in% names(x)) x$Pipeline <- tolower(x$Pipeline)
      if (!is.null(dk) && "Pipeline" %in% names(x)) {
        variant_pipeline_rows[[length(variant_pipeline_rows) + 1]] <<- x %>%
          dplyr::select(dplyr::all_of(dk), genome_build, Pipeline) %>%
          dplyr::distinct() %>%
          dplyr::mutate(elem = elem)
      }
      if (!is.null(dk)) {
        before <- nrow(x)
        x <- dplyr::distinct(x, dplyr::across(dplyr::all_of(dk)), .keep_all = TRUE)
        if (nrow(x) < before) {
          message(sprintf("[write_mutations_db] %s/%s: dropped %d duplicate row(s)",
                          elem, b, before - nrow(x)))
        }
      }
      DBI::dbWriteTable(con, elem, x, append = !first, overwrite = first)
      total <- total + nrow(x)
      first <- FALSE
    }
    total
  }
  counts <- vapply(c("maf", "ashm", "seg", "bedpe"), write_element, integer(1))

  n_variant_pipeline <- 0L
  if (length(variant_pipeline_rows)) {
    variant_pipeline <- dplyr::bind_rows(variant_pipeline_rows)
    DBI::dbWriteTable(con, "variant_pipeline", as.data.frame(variant_pipeline), overwrite = TRUE)
    n_variant_pipeline <- nrow(variant_pipeline)
  }

  # indexes mirroring how the GAMBLR.open accessors query the data.
  # wrapped in try() so a build missing an optional table/column is non-fatal.
  #
  # idx_maf_pipe/idx_vp_pipe are plain indexes on Pipeline -- safe because
  # Pipeline is normalized to lowercase above, at write time, so readers
  # (get_ssm_from_db(), raw SQL) can match it with a plain equality instead
  # of wrapping the column in LOWER()/tolower(), which a normal b-tree index
  # can't be used to satisfy. idx_*_variant_key give the variant_pipeline
  # semi-join (in get_ssm_from_db()'s tool_name resolution) an index to
  # actually use on both sides of the join -- without one, that join has to
  # scan maf/ashm in full for every query with a non-NULL tool_name (the
  # default), which is most of them.
  variant_key_ddl <- "Tumor_Sample_Barcode, Chromosome, Start_Position, End_Position, Tumor_Seq_Allele2"
  idx <- c(
    "CREATE INDEX idx_maf_pos     ON maf(genome_build, Chromosome, Start_Position)",
    "CREATE INDEX idx_maf_sample  ON maf(Tumor_Sample_Barcode)",
    "CREATE INDEX idx_maf_pipe    ON maf(Pipeline)",
    sprintf("CREATE INDEX idx_maf_variant_key ON maf(%s)", variant_key_ddl),
    "CREATE INDEX idx_ashm_pos    ON ashm(genome_build, Chromosome, Start_Position)",
    "CREATE INDEX idx_ashm_sample ON ashm(Tumor_Sample_Barcode)",
    sprintf("CREATE INDEX idx_ashm_variant_key ON ashm(%s)", variant_key_ddl),
    "CREATE INDEX idx_seg_sample  ON seg(genome_build, ID)",
    "CREATE INDEX idx_seg_pos     ON seg(genome_build, chrom, start)",
    "CREATE INDEX idx_bedpe_sample ON bedpe(tumour_sample_id)",
    "CREATE INDEX idx_meta_sample  ON sample_meta(sample_id)",
    "CREATE INDEX idx_meta_barcode ON sample_meta(Tumor_Sample_Barcode)",
    "CREATE INDEX idx_study_sample ON sample_study(sample_id)",
    "CREATE INDEX idx_study_study  ON sample_study(study)",
    "CREATE INDEX idx_vp_sample ON variant_pipeline(Tumor_Sample_Barcode)",
    "CREATE INDEX idx_vp_pipe   ON variant_pipeline(elem, genome_build, Pipeline)",
    sprintf("CREATE INDEX idx_vp_variant_key ON variant_pipeline(%s)", variant_key_ddl)
  )
  for (stmt in idx) try(DBI::dbExecute(con, stmt), silent = TRUE)

  # self-describing provenance / expected-counts table (used by the tests)
  DBI::dbWriteTable(con, "build_info", data.frame(
    key = c("source", "built_at", "builder",
            "n_maf", "n_ashm", "n_seg", "n_bedpe", "n_samples", "n_sample_study", "n_variant_pipeline"),
    value = c(source_desc,
              format(Sys.time(), tz = "UTC", usetz = TRUE),
              "write_mutations_db",
              counts[["maf"]], counts[["ashm"]], counts[["seg"]],
              counts[["bedpe"]], nrow(sample_data$meta), n_sample_study, n_variant_pipeline),
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
