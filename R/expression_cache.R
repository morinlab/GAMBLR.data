#' @title Gene expression, fetched one gene at a time and cached locally
#'
#' @description GAMBL's RNA-seq expression (DESeq2 variance-stabilised
#' values) for the RNA samples linked to the DNA samples in
#' [gambl_mutations_db()], served on demand: the first time a gene is
#' requested, only that gene's line is read from the expression file (a
#' tabix-indexed file, read with an HTTP range request when it is remote) and
#' stored in a local cache; later requests for it are answered from the cache.
#' Nothing else is downloaded except three small tables (samples, genes and
#' build information).
#'
#' @details The source is a folder or URL prefix holding
#' `gambl_expression.tsv.bgz` (+ `.tbi`), `gambl_expression_samples.tsv`,
#' `gambl_expression_genes.tsv` and `gambl_expression_build_info.tsv` (made by
#' `data-raw/build_expression.R`). It is taken from
#' `getOption("GAMBLR.data.expression_source")`, else the
#' `GAMBLR_EXPRESSION_SOURCE` environment variable. The cache is an SQLite file
#' in `tools::R_user_dir("GAMBLR.data", "cache")`; it is cleared
#' automatically when the source's build changes, and by
#' [purge_expression_cache()]. Reading the expression file needs the
#' Bioconductor package Rsamtools.
#'
#' The source can be a local folder, e.g.
#' `options(GAMBLR.data.expression_source = "~/gambl_expression")`. To turn
#' the cache off (every request reads the source; nothing is written to disk),
#' set `options(GAMBLR.data.expression_cache = FALSE)` or the
#' `GAMBLR_EXPRESSION_CACHE` environment variable to `FALSE`.
#'
#' Sample ids: `mrna_sample_id` is the RNA sample's GAMBL sample_id, or, for
#' the few whose id is an internal code, the sample_id of the linked DNA
#' sample; `genome_sample_id` / `capture_sample_id` link each RNA sample to the
#' DNA samples of the same biopsy (an RNA sample_id often differs from them).
#'
#' @param ensembl_gene_ids Ensembl gene ids (without version), e.g.
#'   "ENSG00000172936".
#'
#' @return `get_expression_values()`: a data frame with ensembl_gene_id,
#'   hgnc_symbol, mrna_sample_id and expression, one row per gene and sample.
#'
#' @export
#' @examples
#' \dontrun{
#' options(GAMBLR.data.expression_source = "https://example.org/gambl_expression")
#' myd88 <- get_expression_values("ENSG00000172936")
#' expression_cache_info()
#' }
get_expression_values <- function(ensembl_gene_ids) {
  ids <- unique(sub("\\..*$", "", as.character(unlist(ensembl_gene_ids))))
  ids <- ids[nzchar(ids)]
  if (length(ids) == 0) stop("give at least one ensembl_gene_id", call. = FALSE)
  genes <- expression_genes()
  unknown <- setdiff(ids, genes$ensembl_gene_id)
  if (length(unknown) > 0) {
    stop("not in the expression data: ", paste(utils::head(unknown, 10), collapse = ", "), call. = FALSE)
  }
  if (!expression_cache_enabled()) {
    out <- read_expression_genes(ids, genes)
    if (is.null(out)) out <- data.frame(ensembl_gene_id = character(), mrna_sample_id = character(), expression = numeric())
    out$hgnc_symbol <- genes$hgnc_symbol[match(out$ensembl_gene_id, genes$ensembl_gene_id)]
    return(out[, c("ensembl_gene_id", "hgnc_symbol", "mrna_sample_id", "expression")])
  }
  con <- expression_cache_con()
  on.exit(DBI::dbDisconnect(con))
  cached <- DBI::dbGetQuery(con, "SELECT ensembl_gene_id FROM genes_cached")$ensembl_gene_id
  todo <- setdiff(ids, cached)
  if (length(todo) > 0) {
    rows <- read_expression_genes(todo, genes)
    DBI::dbWithTransaction(con, {
      if (!is.null(rows)) DBI::dbAppendTable(con, "expr", rows)
      DBI::dbAppendTable(con, "genes_cached", data.frame(ensembl_gene_id = todo))
    })
  }
  q <- paste0("SELECT ensembl_gene_id, mrna_sample_id, expression FROM expr WHERE ensembl_gene_id IN (",
              paste(DBI::dbQuoteString(con, ids), collapse = ","), ")")
  out <- DBI::dbGetQuery(con, q)
  out$hgnc_symbol <- genes$hgnc_symbol[match(out$ensembl_gene_id, genes$ensembl_gene_id)]
  out[, c("ensembl_gene_id", "hgnc_symbol", "mrna_sample_id", "expression")]
}

#' @rdname get_expression_values
#' @return `expression_samples()`: the RNA samples (mrna_sample_id, patient_id,
#'   genome_sample_id, capture_sample_id, protocol, ffpe_or_frozen, study and,
#'   in newer builds, vst_zero: the sample's value for a zero count).
#' @export
expression_samples <- function() expression_small_table("samples")

#' @rdname get_expression_values
#' @return `expression_genes()`: the genes (ensembl_gene_id, hgnc_symbol,
#'   chrom, start, end; hg38).
#' @export
expression_genes <- function() expression_small_table("genes")

#' @rdname get_expression_values
#' @return `expression_cache_info()`: a list with the cache file, whether the
#'   cache is on, its size, the build it holds, the source and how many genes
#'   are cached.
#' @export
expression_cache_info <- function() {
  path <- expression_cache_path()
  if (!file.exists(path)) return(list(cache = path, enabled = expression_cache_enabled(), size_mb = 0, genes_cached = 0L,
                                      build_id = NA, source = expression_source(FALSE)))
  con <- DBI::dbConnect(RSQLite::SQLite(), path)
  on.exit(DBI::dbDisconnect(con))
  build <- if (DBI::dbExistsTable(con, "build")) DBI::dbReadTable(con, "build") else data.frame(key = character(), value = character())
  list(cache = path, enabled = expression_cache_enabled(), size_mb = round(file.size(path) / 1e6, 2),
       genes_cached = if (DBI::dbExistsTable(con, "genes_cached")) DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM genes_cached")$n else 0L,
       build_id = build$value[build$key == "build_id"][1], source = expression_source(FALSE))
}

#' @rdname get_expression_values
#' @param genes For `purge_expression_cache()`: ensembl_gene_ids to drop from
#'   the cache, or NULL (default) to delete the whole cache.
#' @return `purge_expression_cache()`: NULL, invisibly.
#' @export
purge_expression_cache <- function(genes = NULL) {
  path <- expression_cache_path()
  if (!file.exists(path)) return(invisible(NULL))
  if (is.null(genes)) {
    unlink(path)
    message("Deleted the expression cache (", path, ")")
    return(invisible(NULL))
  }
  con <- DBI::dbConnect(RSQLite::SQLite(), path)
  on.exit(DBI::dbDisconnect(con))
  ids <- DBI::dbQuoteString(con, sub("\\..*$", "", as.character(genes)))
  DBI::dbExecute(con, paste0("DELETE FROM expr WHERE ensembl_gene_id IN (", paste(ids, collapse = ","), ")"))
  DBI::dbExecute(con, paste0("DELETE FROM genes_cached WHERE ensembl_gene_id IN (", paste(ids, collapse = ","), ")"))
  invisible(NULL)
}

# --- internals --------------------------------------------------------------------------------------------------------

expression_source <- function(required = TRUE) {
  src <- getOption("GAMBLR.data.expression_source", Sys.getenv("GAMBLR_EXPRESSION_SOURCE", unset = ""))
  if (!nzchar(src) && required) {
    stop("No expression source is set. Set options(GAMBLR.data.expression_source = \"<folder or URL>\") or the ",
         "GAMBLR_EXPRESSION_SOURCE environment variable to where the gambl_expression files are.", call. = FALSE)
  }
  sub("/+$", "", src)
}

expression_cache_enabled <- function() {
  opt <- getOption("GAMBLR.data.expression_cache", Sys.getenv("GAMBLR_EXPRESSION_CACHE", unset = "TRUE"))
  !toupper(as.character(opt)[1]) %in% c("FALSE", "F", "0", "NO")
}

is_remote <- function(src) grepl("^(https?|s3|ftp)://", src)

expression_file <- function(name) paste0(expression_source(), "/", name)

expression_cache_path <- function() {
  file.path(tools::R_user_dir("GAMBLR.data", "cache"), "expression_cache.sqlite")
}

#' @rdname get_expression_values
#' @return `expression_build_info()`: the source's build information as a
#'   named character vector (build_id, source_matrix, n_genes, n_samples, ...).
#' @export
expression_build_info <- function() {
  x <- expression_build_table()
  stats::setNames(as.character(x$value), x$key)
}

# The source's build information (read once per session and source, so a rebuilt source invalidates the cache)
expression_build_table <- function() {
  if (is.null(.expr_state$build_info) || !identical(.expr_state$source, expression_source())) {
    f <- expression_file("gambl_expression_build_info.tsv")
    .expr_state$build_info <- utils::read.delim(if (is_remote(f)) url(f) else f, stringsAsFactors = FALSE)
    .expr_state$source <- expression_source()
  }
  .expr_state$build_info
}
.expr_state <- new.env(parent = emptyenv())

# The cache, cleared when it holds another build than the source's
expression_cache_con <- function() {
  path <- expression_cache_path()
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  build_id <- expression_build_info()[["build_id"]]
  con <- DBI::dbConnect(RSQLite::SQLite(), path)
  if (DBI::dbExistsTable(con, "build")) {
    have <- DBI::dbReadTable(con, "build")
    if (!identical(have$value[have$key == "build_id"], build_id)) {
      DBI::dbDisconnect(con)
      unlink(path)
      message("The expression source has a new build; cleared the cached expression")
      con <- DBI::dbConnect(RSQLite::SQLite(), path)
    }
  }
  if (!DBI::dbExistsTable(con, "build")) {
    DBI::dbWriteTable(con, "build", data.frame(key = "build_id", value = build_id))
    DBI::dbExecute(con, "CREATE TABLE expr (ensembl_gene_id TEXT, mrna_sample_id TEXT, expression REAL)")
    DBI::dbExecute(con, "CREATE INDEX idx_expr_gene ON expr(ensembl_gene_id)")
    DBI::dbExecute(con, "CREATE TABLE genes_cached (ensembl_gene_id TEXT PRIMARY KEY)")
  }
  con
}

# The samples or genes table, kept in the cache after the first read (read from the source each time with the cache off)
expression_small_table <- function(which) {
  read_it <- function() {
    f <- expression_file(paste0("gambl_expression_", which, ".tsv"))
    x <- utils::read.delim(if (is_remote(f)) url(f) else f, stringsAsFactors = FALSE, na.strings = c("NA", ""),
                           colClasses = "character")
    if (which == "genes") { x$start <- as.numeric(x$start); x$end <- as.numeric(x$end) }
    if ("vst_zero" %in% names(x)) x$vst_zero <- as.numeric(x$vst_zero)
    x
  }
  if (!expression_cache_enabled()) return(read_it())
  con <- expression_cache_con()
  on.exit(DBI::dbDisconnect(con))
  if (!DBI::dbExistsTable(con, which)) DBI::dbWriteTable(con, which, read_it())
  DBI::dbReadTable(con, which)
}

# Read the genes' lines from the tabix-indexed file: one row per gene and sample (NULL if none)
read_expression_genes <- function(ids, genes) {
  if (!requireNamespace("Rsamtools", quietly = TRUE)) {
    stop("Fetching expression needs the Bioconductor package Rsamtools: BiocManager::install(\"Rsamtools\")", call. = FALSE)
  }
  f <- expression_file("gambl_expression.tsv.bgz")
  tf <- Rsamtools::TabixFile(f, index = paste0(f, ".tbi"))
  header <- sub("^#", "", strsplit(Rsamtools::headerTabix(tf)$header[1], "\t")[[1]])
  samples <- header[-(1:5)]
  g <- genes[match(ids, genes$ensembl_gene_id), , drop = FALSE]
  if (length(ids) > 5) message("Fetching expression for ", length(ids), " genes",
                                if (expression_cache_enabled()) " (cached afterwards)" else "")
  regions <- GenomicRanges::GRanges(g$chrom, IRanges::IRanges(g$start, g$end))
  lines <- unlist(Rsamtools::scanTabix(tf, param = regions), use.names = FALSE)
  fields <- strsplit(unique(lines), "\t", fixed = TRUE)
  fields <- fields[vapply(fields, function(x) x[4] %in% ids, logical(1))]
  do.call(rbind, lapply(fields, function(x) {
    data.frame(ensembl_gene_id = x[4], mrna_sample_id = samples, expression = suppressWarnings(as.numeric(x[-(1:5)])))
  }))
}
