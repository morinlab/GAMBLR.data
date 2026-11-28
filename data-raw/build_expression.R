# build_expression.R
#
# Builds the gene expression files GAMBLR.data serves on demand (one gene at a time, cached locally; see
# R/expression_cache.R), from GAMBL's DESeq2 variance-stabilised (VST) matrix:
#
#   gambl_expression.tsv.bgz (+ .tbi)  one line per gene: chrom, start, end (hg38), ensembl_gene_id, hgnc_symbol,
#                                      then one column per RNA sample (VST, 3 decimals); bgzipped and tabix-indexed
#                                      so a gene's line can be read with an HTTP range request
#   gambl_expression_samples.tsv       one row per RNA sample: mrna_sample_id, patient_id, genome_sample_id,
#                                      capture_sample_id, protocol, ffpe_or_frozen, study, vst_zero (what the VST
#                                      makes of a zero count in that sample: its most common value, shared by the
#                                      hundreds of genes with no reads; not its minimum, since batch correction pushes
#                                      some values below it). Used to tell unexpressed genes apart.
#   gambl_expression_genes.tsv         the gene index: ensembl_gene_id, hgnc_symbol, chrom, start, end
#   gambl_expression_build_info.tsv    provenance (source matrix, its md5, counts, date)
#
# Which RNA samples: the non-redundant RNA samples (one per biopsy, by GAMBLR.results::check_gene_expression()'s
# protocol/preservation priorities) whose biopsy has a genome or capture sample in the bundle's sample_meta. RNA and
# DNA sample_ids often differ for the same biopsy, so the link is check_gene_expression()'s biopsy linkage
# (genome_sample_id / capture_sample_id), never sample_id equality. biopsy_id is internal and not written.
# RNA sample_ids that embed internal lab codes (not already public, not a study barcode, not built from the patient
# id) are replaced by the exact sample_id of the linked DNA sample; every other RNA sample_id is kept as is.
#
# Run from the GAMBLR.data repo root after assemble_bundled_data.R (it reads gambl_mutations.db's sample_meta).
# Needs GAMBLR.results (with its config) and Rsamtools. Set GAMBLR_VST_MATRIX to the wide matrix if it isn't at the
# default local-mirror path.

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})

vst_path <- Sys.getenv("GAMBLR_VST_MATRIX",
  unset = file.path(Sys.getenv("HOME"), "gambl_results/icgc_dart/DESeq2-0.0_salmon-1.0/mrna--gambl-icgc-all/vst-matrix.tsv"))
out_dir <- Sys.getenv("GAMBLR_EXPRESSION_OUT", unset = ".")
stopifnot(file.exists(vst_path), file.exists("gambl_mutations.db"))
if (!requireNamespace("Rsamtools", quietly = TRUE)) stop("build_expression.R needs Rsamtools")

# --- which RNA samples: linked to a bundled DNA sample -----------------------------------------------------------
con <- DBI::dbConnect(RSQLite::SQLite(), "gambl_mutations.db", flags = RSQLite::SQLITE_RO)
sample_meta <- DBI::dbReadTable(con, "sample_meta")
DBI::dbDisconnect(con)
# min_corrected_cov = 0, as the bundle uses, so links to low-coverage DNA samples already bundled aren't missed
links <- suppressMessages(suppressWarnings(
  GAMBLR.results::check_gene_expression(show_linkages = TRUE, min_corrected_cov = 0)))
links <- as.data.frame(links)
links$genome_sample_id[!links$genome_sample_id %in% sample_meta$sample_id] <- NA
links$capture_sample_id[!links$capture_sample_id %in% sample_meta$sample_id] <- NA
links <- links[!is.na(links$genome_sample_id) | !is.na(links$capture_sample_id), , drop = FALSE]
dna_id <- dplyr::coalesce(links$genome_sample_id, links$capture_sample_id)
links$study <- sample_meta$study[match(dna_id, sample_meta$sample_id)]

# public RNA ids: already-public ids, study barcodes and ids built from the patient id are kept; others (internal
# lab codes) are replaced by the exact sample_id of the linked DNA sample (genome if any, else capture)
public_ok <- links$mrna_sample_id %in% sample_meta$sample_id |
  grepl("^(BLGSP|HTMCP)-", links$mrna_sample_id) |
  (mapply(startsWith, links$mrna_sample_id, links$patient_id) & !grepl("CLC", links$mrna_sample_id))
links$public_id <- ifelse(public_ok, links$mrna_sample_id, dna_id)
stopifnot(!anyDuplicated(links$public_id))
message(sprintf("[expression] %d RNA samples linked to bundled DNA samples (%d use their DNA sample_id in place of an internal RNA id); by study: %s",
                nrow(links), sum(!public_ok),
                paste(names(table(links$study)), table(links$study), sep = " ", collapse = ", ")))

# --- the matrix: kept samples only, hg38 coordinates, one row per gene --------------------------------------------
header <- names(fread(vst_path, nrows = 0))
missing <- setdiff(links$mrna_sample_id, header)
if (length(missing) > 0) stop(length(missing), " linked RNA samples aren't columns of the VST matrix")
m <- fread(vst_path, select = c("gene_id", "ensembl_gene_id", "hgnc_symbol", links$mrna_sample_id))
# GENCODE annotates each pseudoautosomal gene twice (X, and a "_PAR_Y" copy with the same ensembl_gene_id): keep X
n_par_y <- sum(grepl("_PAR_Y$", m$gene_id))
m <- m[!grepl("_PAR_Y$", m$gene_id)][, gene_id := NULL]
setnames(m, links$mrna_sample_id, links$public_id)
coords <- as.data.table(GAMBLR.data::hg38_gene_coordinates)[, .(ensembl_gene_id, chrom = sub("^chr", "", chromosome),
                                                                 start, end)]
# a gene with several coordinate rows (e.g. pseudoautosomal X and Y copies): keep one, preferring X over Y
coords <- coords[order(ensembl_gene_id, chrom == "Y", chrom, start)][!duplicated(ensembl_gene_id)]
m <- merge(coords, m, by = "ensembl_gene_id")
if (anyDuplicated(m$ensembl_gene_id)) stop("duplicated ensembl_gene_id after the coordinate join")
for (s in links$public_id) set(m, j = s, value = round(m[[s]], 3))
vst_zero <- vapply(links$public_id, function(s) {
  tab <- table(m[[s]])
  as.numeric(names(tab)[which.max(tab)])
}, numeric(1))
setcolorder(m, c("chrom", "start", "end", "ensembl_gene_id", "hgnc_symbol", links$public_id))
setorder(m, chrom, start, end)
message(sprintf("[expression] %d genes x %d samples", nrow(m), length(links$public_id)))

tsv <- file.path(out_dir, "gambl_expression.tsv")
fwrite(m, tsv, sep = "\t", quote = FALSE, na = "NA")
# tabix wants the column header as a comment line
system2("sed", c("-i.bak", "'1s/^/#/'", shQuote(tsv)))
unlink(paste0(tsv, ".bak"))
bgz <- Rsamtools::bgzip(tsv, dest = file.path(out_dir, "gambl_expression.tsv.bgz"), overwrite = TRUE)
Rsamtools::indexTabix(bgz, seq = 1, start = 2, end = 3, comment = "#")
unlink(tsv)

# --- small tables ---------------------------------------------------------------------------------------------------
samples <- data.frame(mrna_sample_id = links$public_id, patient_id = links$patient_id,
                      genome_sample_id = links$genome_sample_id, capture_sample_id = links$capture_sample_id,
                      protocol = links$protocol, ffpe_or_frozen = links$ffpe_or_frozen, study = links$study,
                      vst_zero = unname(vst_zero))
fwrite(samples, file.path(out_dir, "gambl_expression_samples.tsv"), sep = "\t", na = "NA")
fwrite(m[, .(ensembl_gene_id, hgnc_symbol, chrom, start, end)], file.path(out_dir, "gambl_expression_genes.tsv"),
       sep = "\t", na = "NA")
# the build is the expression file and its small tables together
build_id <- unname(tools::md5sum(bgz))
build_id <- paste(substr(build_id, 1, 16),
                  substr(unname(tools::md5sum(file.path(out_dir, "gambl_expression_samples.tsv"))), 1, 8),
                  substr(unname(tools::md5sum(file.path(out_dir, "gambl_expression_genes.tsv"))), 1, 8), sep = "")
fwrite(data.frame(key = c("build_id", "source_matrix", "source_md5", "n_genes", "n_samples", "transformation",
                          "par_y_rows_dropped", "built_at"),
                  value = c(build_id, basename(dirname(vst_path)), unname(tools::md5sum(vst_path)), nrow(m),
                            nrow(samples), "DESeq2 variance-stabilising transformation (vst), 3 decimals",
                            n_par_y, format(Sys.time(), tz = "UTC", usetz = TRUE))),
       file.path(out_dir, "gambl_expression_build_info.tsv"), sep = "\t")
message(sprintf("[expression] wrote %s (%.0f MB) and its index, samples, genes and build info; build_id %s",
                bgz, file.size(bgz) / 1e6, build_id))
