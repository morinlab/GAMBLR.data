# test_gambl_db.R
#
# Simple, dependency-light sanity checks for a freshly built gambl_mutations.db.
# Run after building the DB (on the GSC side) to confirm the output is usable by
# the GAMBLR.open accessors.
#
# Usage:
#   Rscript data-raw/test_gambl_db.R [OUT_DB] [SAMPLE_METADATA_RDA] [SAMPLE_DATA_RDA]
#
# Exits non-zero if any critical check fails. The optional SAMPLE_DATA_RDA
# enables a regression check (DB row counts vs the in-memory frames).

suppressMessages({ library(DBI); library(RSQLite) })

args   <- commandArgs(trailingOnly = TRUE)
db     <- if (length(args) >= 1) args[[1]] else "gambl_mutations.db"
meta   <- if (length(args) >= 2) args[[2]] else "data/sample_metadata.rda"
rda    <- if (length(args) >= 3) args[[3]] else "data/sample_data.rda"

fails <- 0L
check <- function(cond, msg) {
  status <- if (isTRUE(cond)) "PASS" else { fails <<- fails + 1L; "FAIL" }
  cat(sprintf("  [%s] %s\n", status, msg))
}

stopifnot(file.exists(db))
con <- dbConnect(SQLite(), db, flags = SQLITE_RO)
on.exit(dbDisconnect(con))
tbls <- dbListTables(con)
n    <- function(t) dbGetQuery(con, sprintf("SELECT COUNT(*) n FROM %s", t))$n
cols <- function(t) dbListFields(con, t)
builds_in <- function(t) sort(dbGetQuery(con, sprintf("SELECT DISTINCT genome_build g FROM %s", t))$g)

cat("== tables present ==\n")
for (t in c("maf","ashm","seg","bedpe","sample_meta","sample_study","variant_pipeline","build_info"))
  check(t %in% tbls, sprintf("table %s exists", t))

cat("\n== row counts (all > 0) ==\n")
for (t in c("maf","ashm","seg","bedpe","sample_meta","sample_study","variant_pipeline")) {
  cnt <- if (t %in% tbls) n(t) else 0
  check(cnt > 0, sprintf("%-11s has %d rows", t, cnt))
}

cat("\n== both genome builds present ==\n")
for (t in c("maf","ashm","seg","bedpe"))
  if (t %in% tbls)
    check(all(c("grch37","hg38") %in% builds_in(t)),
          sprintf("%-6s has builds: %s", t, paste(builds_in(t), collapse=", ")))

cat("\n== required columns present (what the accessors filter on) ==\n")
req <- list(
  maf   = c("Chromosome","Start_Position","End_Position","Tumor_Sample_Barcode",
            "Pipeline","Hugo_Symbol","Variant_Classification","t_alt_count","genome_build"),
  seg   = c("ID","chrom","start","end","CN","genome_build"),
  bedpe = c("tumour_sample_id","CHROM_A","START_A","CHROM_B","START_B","VAF_tumour","SCORE","FILTER","genome_build"),
  sample_meta = c("sample_id","Tumor_Sample_Barcode","seq_type","study"),
  sample_study = c("sample_id","study","study_id","reference_PMID"),
  variant_pipeline = c("Tumor_Sample_Barcode","Chromosome","Start_Position","End_Position",
                       "Tumor_Seq_Allele2","genome_build","elem","Pipeline")
)
for (t in names(req)) if (t %in% tbls) {
  missing <- setdiff(req[[t]], cols(t))
  check(length(missing) == 0,
        sprintf("%-11s columns present%s", t,
                if (length(missing)) paste0(" (MISSING: ", paste(missing, collapse=", "), ")") else ""))
}

cat("\n== indexes present ==\n")
# NB: idx_maf_gene was asserted here previously but was never actually
# created by write_mutations_db.R (maf has no Hugo_Symbol index -- gene
# queries resolve to a region first, see gambl_mutations_db()'s docs), so
# this check always failed regardless of any other change. Fixed here to
# assert indexes that actually exist.
idx <- dbGetQuery(con, "SELECT name FROM sqlite_master WHERE type='index'")$name
for (i in c("idx_maf_pos","idx_maf_sample","idx_maf_pipe","idx_maf_variant_key",
            "idx_ashm_variant_key","idx_seg_sample",
            "idx_bedpe_sample","idx_study_sample","idx_study_study",
            "idx_vp_sample","idx_vp_pipe","idx_vp_variant_key"))
  check(i %in% idx, sprintf("index %s", i))

cat("\n== pipelines present ==\n")
# Pipeline is normalized to lowercase at write time (see write_mutations_db.R)
# so it can be matched with a plain, indexed equality instead of every query
# needing to wrap the column in LOWER()/tolower().
pipes <- dbGetQuery(con, "SELECT DISTINCT Pipeline p FROM maf")$p
for (p in c("slms-3","publication"))
  check(p %in% pipes, sprintf("Pipeline '%s' present (all: %s)", p, paste(pipes, collapse=", ")))

cat("\n== spot query: MYC locus, grch37, slms-3 returns variants ==\n")
myc <- dbGetQuery(con, "SELECT COUNT(*) n FROM maf WHERE genome_build='grch37'
  AND Chromosome='8' AND Start_Position>128723128 AND Start_Position<128774067
  AND Pipeline='slms-3'")$n
check(myc > 0, sprintf("MYC slms-3 variants found: %d", myc))

cat("\n== build_info counts match actual table counts ==\n")
if ("build_info" %in% tbls) {
  bi <- dbGetQuery(con, "SELECT key, value FROM build_info")
  getbi <- function(k) as.numeric(bi$value[bi$key == k])
  for (t in c("maf","ashm","seg","bedpe"))
    check(isTRUE(getbi(paste0("n_", t)) == n(t)),
          sprintf("build_info n_%s == COUNT(%s)", t, t))
  check(isTRUE(getbi("n_variant_pipeline") == n("variant_pipeline")),
        sprintf("build_info n_variant_pipeline == COUNT(variant_pipeline)"))
}

cat("\n== sample_metadata.rda matches sample_meta table ==\n")
if (file.exists(meta)) {
  e <- new.env(); load(meta, envir = e); sm <- get("sample_metadata", envir = e)
  check(nrow(sm) == n("sample_meta"),
        sprintf("sample_metadata rows (%d) == sample_meta table (%d)", nrow(sm), n("sample_meta")))
} else cat("  [skip] no", meta, "\n")

cat("\n== regression vs sample_data.rda (optional) ==\n")
if (file.exists(rda)) {
  e <- new.env(); load(rda, envir = e); sd <- get("sample_data", envir = e)
  for (t in c("maf","ashm","seg","bedpe")) {
    expect <- sum(vapply(c("grch37","hg38"),
                         function(b) { x <- sd[[b]][[t]]; if (is.null(x)) 0L else nrow(x) }, integer(1)))
    check(expect == n(t), sprintf("%-6s DB rows (%d) == sample_data rows (%d)", t, n(t), expect))
  }
} else cat("  [skip] no", rda, "(fresh GSC build has no legacy .rda — expected)\n")

cat(sprintf("\n==== %s: %d check(s) failed ====\n", if (fails==0) "ALL PASSED" else "FAILURES", fails))
if (fails > 0) quit(status = 1)
