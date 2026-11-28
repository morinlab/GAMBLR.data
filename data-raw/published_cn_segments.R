# published_cn_segments.R
#
# Builds inst/extdata/published_cn_segments.tsv.gz: copy number segments from
# published analyses of cohorts that are also in GAMBL, under GAMBL sample_ids,
# in grch37 and hg38, for comparison with GAMBL's own segmentations. Read by
# get_published_cn_segments(). Inputs are all in data-raw/published_cn/ (see the
# README there); run from the repository root:
#   Rscript data-raw/published_cn_segments.R
# Needs data.table and rtracklayer (build time only).
#
# Sources (`source` column):
#   dlbclass       DLBClass (Chapuy et al.) segments, DLBCL_699_segs.2021-12-15.ccf.seg:
#                  the Broad's exome segmentation (Chapuy cohort) and reprocessed
#                  arrays (Schmitz cohort), with ABSOLUTE purity and ploidy.
#                  log.ratio = Segment_Mean (log2(tau / 2), what the DLBClass GSM
#                  is built from). CN = ABSOLUTE's purity/ploidy-corrected copy
#                  number for every segment:
#                    (2^log.ratio * D - 2 * (1 - purity)) / purity,
#                    D = purity * ploidy + 2 * (1 - purity)
#                  (floored at 0); on clonal altered segments it is within 0.5 of
#                  ABSOLUTE's allelic total (NA + NB) for 88% of segments.
#                  LOH_flag = 1 where ABSOLUTE calls an altered segment (IS_SCNA 1)
#                  with no minor allele (NA 0), 0 for other altered/unaltered
#                  calls, NA where there are no allelic calls (all Schmitz
#                  samples).
#   schmitz_array  Schmitz et al. array copy number (schmitz_cnv_fix.seg):
#                  CN = Estimated_Raw_Copy, log.ratio = log2(CN / 2) (CN floored
#                  at 0.05 so complete losses get a finite log ratio).
#
# Only samples with a GAMBL sample_id (published_cn/*_sample_map.tsv) are kept.
# Coordinates are as in the source files (grch37). hg38 coordinates come from
# lifting each segment's start and end separately (inst/extdata chain; an end
# that doesn't map, e.g. at a chromosome tip, is stepped inward by up to 5 Mb
# until it does), keeping segments whose ends map to the same chromosome in
# order; the rest (about 0.5%) are dropped (counted below).

suppressMessages({
  library(data.table)
  library(rtracklayer)
  library(GenomicRanges)
})

dir_in = "data-raw/published_cn"
out_file = "inst/extdata/published_cn_segments.tsv.gz"
columns = c("source", "genome_build", "ID", "chrom", "start", "end", "LOH_flag", "CN",
            "log.ratio", "source_id", "purity", "ploidy")

# ---- dlbclass ----
d = fread(file.path(dir_in, "raw", "DLBCL_699_segs.2021-12-15.ccf.seg.gz"), na.strings = "",
          colClasses = "character")
setnames(d, make.names(names(d), unique = TRUE))  # "NA" -> "NA."
for(k in c("Start", "End", "Segment_Mean", "purity", "ploidy", "NA.", "NB", "IS_SCNA")){
  set(d, j = k, value = as.numeric(d[[k]]))
}
map = fread(file.path(dir_in, "dlbclass_sample_map.tsv"))
d = d[Sample %in% map$source_id]
D = d$purity * d$ploidy + 2 * (1 - d$purity)
dlbclass = d[, .(source = "dlbclass", genome_build = "grch37",
                 ID = map$sample_id[match(Sample, map$source_id)],
                 chrom = sub("^23$", "X", sub("^24$", "Y", Chromosome)),
                 start = Start, end = End,
                 LOH_flag = fifelse(is.na(IS_SCNA), NA_integer_,
                                    fifelse(IS_SCNA == 1 & NA. %in% 0, 1L, 0L)),
                 CN = round(pmax((2^Segment_Mean * D - 2 * (1 - purity)) / purity, 0), 3),
                 log.ratio = Segment_Mean, source_id = Sample, purity, ploidy)]

# ---- schmitz_array ----
s = fread(file.path(dir_in, "raw", "schmitz_cnv_fix.seg.gz"))
smap = fread(file.path(dir_in, "schmitz_array_sample_map.tsv"))
s = s[Subject_Name %in% smap$source_id]
schmitz = s[, .(source = "schmitz_array", genome_build = "grch37",
                ID = smap$sample_id[match(Subject_Name, smap$source_id)],
                chrom = as.character(Chromosome), start = hg19.Start, end = hg19.End,
                LOH_flag = NA_integer_, CN = round(Estimated_Raw_Copy, 3),
                log.ratio = round(log2(pmax(Estimated_Raw_Copy, 0.05) / 2), 4),
                source_id = Subject_Name, purity = NA_real_, ploidy = NA_real_)]

grch37 = rbind(dlbclass, schmitz)
setorder(grch37, source, ID, chrom, start)

# ---- hg38: lift each end separately ----
chain = import.chain("inst/extdata/hg19ToHg38.over.chain")
lift_points = function(chrom, pos){
  gr = GRanges(paste0("chr", chrom), IRanges(pos, width = 1), i = seq_along(pos))
  hits = unlist(liftOver(gr, chain))
  hits = hits[!duplicated(hits$i)]  # first mapping if a point maps more than once
  out = data.table(chrom = rep(NA_character_, length(pos)), pos = NA_real_)
  out[hits$i, `:=`(chrom = as.character(seqnames(hits)), pos = start(hits))]
  out
}
# An end that doesn't map (typically a large segment starting or ending at a
# chromosome tip, outside the chain) is stepped inward until a position maps,
# as long as it stays inside the segment, so whole-arm segments aren't lost
lift_end = function(chrom, pos, other, inward){
  out = lift_points(chrom, pos)
  for(step in c(1e4, 5e4, 1e5, 5e5, 1e6, 2e6, 5e6)){
    todo = which(is.na(out$chrom))
    if(length(todo) == 0) break
    moved = pos[todo] + inward * step
    inside = if(inward > 0) moved < other[todo] else moved > other[todo]
    todo = todo[inside]; moved = moved[inside]
    if(length(todo) == 0) next
    hit = lift_points(chrom[todo], moved)
    out$chrom[todo] = hit$chrom
    out$pos[todo] = hit$pos
  }
  out
}
a = lift_end(grch37$chrom, grch37$start, grch37$end, inward = 1)
b = lift_end(grch37$chrom, grch37$end, grch37$start, inward = -1)
ok = !is.na(a$chrom) & !is.na(b$chrom) & a$chrom == b$chrom & b$pos > a$pos
ok[is.na(ok)] = FALSE
# plain vectors, assigned outside data.table's [ ] (where names can resolve to columns)
keep = which(ok)
hg38 = copy(grch37[keep])
hg38$genome_build = "hg38"
hg38$chrom = a$chrom[keep]
hg38$start = a$pos[keep]
hg38$end = b$pos[keep]
dropped = grch37[!ok, .N, by = source]
message("hg38 liftover dropped ", sum(!ok), " of ", nrow(grch37), " segments (",
        paste(dropped$source, dropped$N, collapse = ", "), ")")

# Segments of a sample must not overlap (they would be double counted).
# Liftover can make neighbouring segments overlap (rearrangements between
# builds, ends stepped inward); as in GAMBLR.utils::resolve_segment_overlaps
# (method "trim"), a segment starting before the furthest end so far on its
# chromosome is trimmed to start there, and dropped if nothing is left.
resolve_overlaps = function(x){
  x = x[order(source, ID, chrom, start, end)]
  # data.table::shift explicitly: GenomicRanges, loaded above, masks shift()
  prev_reach = x[, data.table::shift(cummax(end)), by = .(source, ID, chrom)]$V1
  trim = !is.na(prev_reach) & x$start < prev_reach
  x$start[trim] = prev_reach[trim]
  empty = x$end <= x$start
  message(x$genome_build[1], ": trimmed ", sum(trim & !empty), " overlapping segment(s), dropped ",
          sum(empty), " left empty")
  x[!empty]
}
grch37 = resolve_overlaps(grch37)
hg38 = resolve_overlaps(hg38)

out = rbind(grch37, hg38)[, ..columns]
fwrite(out, out_file, sep = "\t", na = "NA")
message("wrote ", out_file, ": ", nrow(out), " rows; samples: ",
        paste(out[genome_build == "grch37", .(n = uniqueN(ID)), by = source][, paste(source, n)], collapse = ", "))
