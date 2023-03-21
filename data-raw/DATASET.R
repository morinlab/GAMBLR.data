#!/usr/bin/env Rscript

library(usethis)
library(dplyr)

# Helper function to input the data objects
generate_data <- function(
    Row
    ){

    # What is the name of the object?
    data_object <- paste0(
            Row[["dataset"]],
            "_",
            Row[["genome_build"]],
            "_v",
            Row[["version"]]
        )

    # Where is it located?
    data_file <- paste0(
            getwd(),
            "/inst/extdata/",
            Row[["dataset"]],
            "/",
            Row[["version"]],
            "/",
            Row[["genome_build"]],
            ".tsv"
        )

    # Use the name of the object to store the associated data
    assign(
        data_object,
        read.delim(
            data_file
        )
    )
    # Generate .rda object with the data in it
    do.call(
        "use_data",
        list(as.name(data_object),
        overwrite = TRUE)
    )
}

# All possible names of datasets
datasets <- c(
    "somatic_hypermutation_locations"
)

# All possible genome builds
genome_builds <- c(
    "GRCh37",
    "GRCh38"
)

# All possible versions
# This will need to be modified with the release of new data
versions <- c(
    "0.0",
    "0.1",
    "_latest"
)

# All possible combinations of data/genome build/version
all <- expand.grid(
    dataset = datasets,
    genome_build = genome_builds,
    version = versions,
    stringsAsFactors = FALSE)

all <- all[apply(all, 1, function(x) {length(unique(x)) == 3}),]

# Generate the rdas for each object
apply(
    all,
    1,
    generate_data
)
