# make_sample_metadata.R
#
# Extract the lightweight sample metadata (sample_data$meta) into its own bundled
# object, `sample_metadata`, so that loading metadata does NOT drag the multi-GB
# MAF frames of `sample_data` into memory. The heavy mutation frames now live in
# gambl_mutations.db (see data-raw/build_mutations_db.R).

e <- new.env()
load(file.path("data", "sample_data.rda"), envir = e)
sample_metadata <- get("sample_data", envir = e)$meta

message(sprintf("sample_metadata: %d rows x %d cols (%.0f KB)",
                nrow(sample_metadata), ncol(sample_metadata),
                as.numeric(object.size(sample_metadata)) / 1024))

save(sample_metadata, file = file.path("data", "sample_metadata.rda"), compress = "xz")
