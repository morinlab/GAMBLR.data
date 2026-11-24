#' @title Copy number segments from published analyses
#'
#' @description Segments from published analyses of cohorts that are also in
#' GAMBL, under GAMBL sample_ids, for comparing GAMBL's own segmentations with
#' them. Backs `get_cn_segments(pipeline = ...)` in GAMBLR.open and
#' GAMBLR.results. See [published_cn_sources] for the sources.
#'
#' @details Built by `data-raw/published_cn_segments.R` from the files in
#' `data-raw/published_cn/` (sources, checksums and the sample maps are in the
#' README there) and shipped as `inst/extdata/published_cn_segments.tsv.gz`.
#' Coordinates are as in the source files (grch37); hg38 coordinates were
#' lifted over segment end by segment end, and about 0.5% of segments that
#' couldn't be lifted are missing from hg38. The table is read once per
#' session.
#'
#' @param source Optional source name(s), from [published_cn_sources]. Default
#' NULL: all.
#' @param projection Genome build: "grch37" (default) or "hg38" (chromosomes
#' then have the "chr" prefix).
#' @param sample_ids Optional GAMBL sample_ids to keep.
#'
#' @return A data frame with source, ID (GAMBL sample_id), chrom, start, end,
#' LOH_flag, CN, log.ratio, source_id (the sample's ID in the source),
#' purity and ploidy (where the source has them).
#'
#' @export
#'
#' @examples
#' published_cn_sources()
#' segs = get_published_cn_segments("dlbclass", sample_ids = "DLBCL-RICOVER_299-Tumor")
#' head(segs)
get_published_cn_segments = function(source = NULL, projection = "grch37", sample_ids = NULL){
  if(!projection %in% c("grch37", "hg38")) stop("projection must be \"grch37\" or \"hg38\"")
  segs = .published_cn_table()
  if(!is.null(source)){
    unknown = setdiff(source, published_cn_sources()$source)
    if(length(unknown) > 0){
      stop("unknown source(s): ", paste(unknown, collapse = ", "), "; see published_cn_sources()")
    }
  }
  keep = segs$genome_build == projection &
    (if(is.null(source)) TRUE else segs$source %in% source) &
    (if(is.null(sample_ids)) TRUE else segs$ID %in% sample_ids)
  out = segs[keep, setdiff(colnames(segs), "genome_build"), drop = FALSE]
  rownames(out) = NULL
  out
}

#' @title Sources of published copy number segments
#'
#' @description The sources available through [get_published_cn_segments],
#' each usable as a `pipeline` name in `get_cn_segments()`.
#'
#' @return A data frame with source, seq_type (of the GAMBL samples they are
#' attached to), the number of samples, and a description.
#'
#' @export
#'
#' @examples
#' published_cn_sources()
published_cn_sources = function(){
  segs = .published_cn_table()
  n = tapply(segs$ID[segs$genome_build == "grch37"], segs$source[segs$genome_build == "grch37"],
             function(x) length(unique(x)))
  data.frame(
    source = c("dlbclass", "schmitz_array"),
    seq_type = c("capture", "capture"),
    n_samples = as.integer(n[c("dlbclass", "schmitz_array")]),
    description = c(
      paste("DLBClass (Chapuy et al.): the Broad's exome segmentation (Chapuy) and reprocessed",
            "arrays (Schmitz); log.ratio = Segment_Mean, CN = ABSOLUTE purity/ploidy-corrected",
            "copy number"),
      paste("Schmitz et al. array copy number; CN = estimated raw copy number,",
            "log.ratio = log2(CN / 2)")),
    stringsAsFactors = FALSE)
}

.published_cn_cache = new.env(parent = emptyenv())

.published_cn_table = function(){
  if(is.null(.published_cn_cache$segs)){
    path = system.file("extdata", "published_cn_segments.tsv.gz", package = "GAMBLR.data")
    if(path == "") stop("published_cn_segments.tsv.gz is missing from GAMBLR.data")
    .published_cn_cache$segs = utils::read.delim(
      gzfile(path), stringsAsFactors = FALSE, na.strings = "NA",
      colClasses = c(source = "character", genome_build = "character", ID = "character",
                     chrom = "character", start = "numeric", end = "numeric",
                     LOH_flag = "integer", CN = "numeric", log.ratio = "numeric",
                     source_id = "character", purity = "numeric", ploidy = "numeric"))
  }
  .published_cn_cache$segs
}
