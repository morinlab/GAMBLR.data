# build_mutations_db.R
#
# ONE-TIME adapter: build gambl_mutations.db from an EXISTING sample_data.rda,
# without re-running the (GSC-dependent) assembly. Use this only to convert a
# pre-existing bundled object.
#
# The ONGOING build writes the DB directly from source at the end of
# data-raw/assemble_bundled_data.R via the same write_mutations_db() helper, so
# there is no sample_data.rda round-trip in the normal pipeline.
#
# Usage:
#   Rscript data-raw/build_mutations_db.R [SAMPLE_DATA_RDA] [OUT_DB]

source("data-raw/write_mutations_db.R")

args <- commandArgs(trailingOnly = TRUE)
rda <- if (length(args) >= 1) args[[1]] else "data/sample_data.rda"
out <- if (length(args) >= 2) args[[2]] else "gambl_mutations.db"

message("Loading ", rda, " ...")
e <- new.env()
load(rda, envir = e)
sample_data <- get("sample_data", envir = e)

write_mutations_db(sample_data, out_db = out,
                   source_desc = paste("one-time adapter from", rda))
