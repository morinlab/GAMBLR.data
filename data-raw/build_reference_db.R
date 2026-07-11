# build_reference_db.R
#
# Proof-of-concept: build a single normalized SQLite reference database from the
# canonical LLMPP curated gene lists (NOT the versioned lymphoma_genes_*.rda
# objects). This one file is intended to become the shared backing store for
# GAMBLR.data accessors, the Lymphopedia wiki, and the Lymphopedia JSON/SQL API.
#
# Usage:
#   Rscript data-raw/build_reference_db.R [LLMPP_CURATED_DIR] [OUT_DB]
#
# Design notes:
#   * The per-entity curated TSVs have heterogeneous schemas. We normalize the
#     common core into one long `gene_entity` table (one row per gene x entity)
#     and stash the entity-specific columns as JSON in `extra` so nothing is lost.
#   * "Which version" becomes a column (source_version), not a filename
#     convention -- this collapses the lymphoma_genes_*_v0.1/v0.2/v_latest soup.

suppressMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(DBI)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
curated_dir <- ifelse(length(args) >= 1, args[[1]],
                      path.expand("~/git/LLMPP/resources/curated"))
out_db <- ifelse(length(args) >= 2, args[[2]],
                 file.path("inst", "extdata", "gambl_reference.db"))

stopifnot(dir.exists(curated_dir))
dir.create(dirname(out_db), recursive = TRUE, showWarnings = FALSE)

# LLMPP source version -> stamped into every row so the data is self-describing
llmpp_version <- tryCatch(
  system(paste("git -C", shQuote(dirname(dirname(curated_dir))),
               "describe --always --dirty"), intern = TRUE),
  error = function(e) NA_character_
)
if (length(llmpp_version) == 0 || is.na(llmpp_version[1])) {
  llmpp_version <- paste0("LLMPP-", format(Sys.Date()))
} else {
  llmpp_version <- paste0("LLMPP@", llmpp_version[1])
}

# entity -> curated file
entity_files <- c(
  DLBCL = "dlbcl_genes.tsv",
  BL    = "bl_genes.tsv",
  FL    = "fl_genes.tsv",
  PMBL  = "pmbl_genes.tsv",
  MCL   = "mcl_genes.tsv",
  MZL   = "mzl_genes.tsv"
)

# pull a column by any of several candidate names (case-insensitive); NA if absent
pick <- function(df, candidates) {
  hit <- names(df)[tolower(names(df)) %in% tolower(candidates)]
  if (length(hit) == 0) return(rep(NA_character_, nrow(df)))
  as.character(df[[hit[1]]])
}

core_cols <- c("gene", "tier", "ashm", "qc", "citekey", "pmid", "mutation_effect")

harmonize_one <- function(entity, file) {
  df <- suppressWarnings(read_tsv(file, col_types = cols(.default = "c"),
                                  name_repair = "unique_quiet"))
  core <- tibble(
    gene            = pick(df, c("Gene", "gene", "Hugo_Symbol")),
    entity          = entity,
    tier            = pick(df, c("Tier")),
    ashm            = pick(df, c("aSHM")),
    qc              = pick(df, c("QC")),
    citekey         = pick(df, c("citekey")),
    pmid            = pick(df, c("PMID")),
    mutation_effect = pick(df, c("MutationEffect"))
  ) |> filter(!is.na(gene), gene != "")
  # keep everything else as JSON so entity-specific fields are not lost
  used <- c("Gene", "gene", "Hugo_Symbol", "Tier", "aSHM", "QC", "citekey",
            "PMID", "MutationEffect")
  extra_df <- df[, !tolower(names(df)) %in% tolower(used), drop = FALSE]
  extra_df <- extra_df[, names(extra_df) != "" & !grepl("^\\.\\.\\.", names(extra_df)),
                       drop = FALSE]
  core$extra <- if (ncol(extra_df) == 0) NA_character_ else
    vapply(seq_len(nrow(core)), function(i)
      toJSON(as.list(extra_df[i, , drop = FALSE]), auto_unbox = TRUE, na = "null"),
      character(1))
  core
}

# symbol -> ensembl map, reused from GAMBLR.data's own gencode_to_symbol object,
# so gene_format = "ensembl" keeps working from the DB. One ensembl id per symbol.
sym2ens <- local({
  e <- new.env(); load(file.path("data", "gencode_to_symbol.rda"), envir = e)
  get("gencode_to_symbol", envir = e) |>
    dplyr::filter(!is.na(hgnc_symbol), hgnc_symbol != "") |>
    dplyr::distinct(hgnc_symbol, .keep_all = TRUE) |>
    dplyr::select(gene = hgnc_symbol, ensembl_gene_id)
})

gene_entity <- bind_rows(
  Map(function(e, f) harmonize_one(e, file.path(curated_dir, f)),
      names(entity_files), entity_files)
) |>
  left_join(sym2ens, by = "gene") |>
  mutate(
    tier = suppressWarnings(as.integer(tier)),
    source_version = llmpp_version
  ) |>
  relocate(ensembl_gene_id, .after = gene) |>
  relocate(source_version, .after = last_col())

# --- second reference table: aSHM regions (both genome builds) ----------------
read_ashm <- function(file, build) {
  df <- suppressWarnings(read_tsv(file, col_types = cols(.default = "c"),
                                  name_repair = "unique_quiet"))
  tibble(
    gene               = pick(df, c("gene")),
    chrom              = pick(df, c("chr_name", "chrom", "chr")),
    start              = as.integer(df[[2]]),
    end                = as.integer(df[[3]]),
    region             = pick(df, c("region")),
    regulatory_comment = pick(df, c("regulatory_comment")),
    build              = build,
    source_version     = llmpp_version
  ) |> filter(!is.na(gene))
}
ashm_dir <- file.path(curated_dir, "aSHM")
ashm_regions <- bind_rows(
  read_ashm(file.path(ashm_dir, "somatic_hypermutation_locations_GRCh37.txt"), "grch37"),
  read_ashm(file.path(ashm_dir, "somatic_hypermutation_locations_hg38.txt"),   "hg38")
)

# --- write the database -------------------------------------------------------
if (file.exists(out_db)) file.remove(out_db)
con <- dbConnect(RSQLite::SQLite(), out_db)
on.exit(dbDisconnect(con))

dbWriteTable(con, "gene_entity", as.data.frame(gene_entity), overwrite = TRUE)
dbWriteTable(con, "ashm_regions", as.data.frame(ashm_regions), overwrite = TRUE)
dbWriteTable(con, "build_info", data.frame(
  key   = c("source", "source_version", "built_at", "builder"),
  value = c("LLMPP resources/curated", llmpp_version,
            format(Sys.time(), tz = "UTC", usetz = TRUE),
            "data-raw/build_reference_db.R")
), overwrite = TRUE)

dbExecute(con, "CREATE INDEX idx_ge_gene ON gene_entity(gene)")
dbExecute(con, "CREATE INDEX idx_ge_entity ON gene_entity(entity)")
dbExecute(con, "CREATE INDEX idx_ge_tier ON gene_entity(tier)")
dbExecute(con, "CREATE INDEX idx_ashm_gene ON ashm_regions(gene)")

cat(sprintf("Wrote %s\n  gene_entity : %d rows (%d genes, %d entities)\n  ashm_regions: %d rows\n  version     : %s\n",
            out_db, nrow(gene_entity), dplyr::n_distinct(gene_entity$gene),
            dplyr::n_distinct(gene_entity$entity), nrow(ashm_regions), llmpp_version))
