#load libraries
library(testthat)

test_that("Check for consistency in rows and columns", {
  expect_equal(nrow(get_manta_sv()), 1045) #all default parameters
  expect_equal(nrow(get_manta_sv(these_samples_metadata = get_gambl_metadata(seq_type_filter = "genome"))), 1045) #provide metadata for genome samples (this is default behavior)
  expect_equal(nrow(get_manta_sv(projection = "hg38")), 1083) #all default parameters, except projection
  expect_equal(nrow(get_manta_sv(these_samples_metadata = get_gambl_metadata() %>% head(100))), 4) #return manta calls for the first 100 samples in the bundled metadata (4 samples should have manta calls)
  expect_equal(nrow(unique(get_manta_sv()["tumour_sample_id"])), 520) #how many samples do we get manta calls for with default parameters
  expect_equal(nrow(unique(get_manta_sv(projection = "hg38")["tumour_sample_id"])), 532) #how many samples from hg38 do we have manta calls for
  expect_equal(nrow(get_manta_sv(region = "8:128723128-128774067")), 245) #expected manta calls for a given region
  expect_equal(ncol(get_manta_sv()), 16) #expected number of columns returned with default parameter
  expect_equal(ncol(get_manta_sv(these_sample_ids = "DOHH-2")), 16) #expected number of columns returned with common alternate parameters
  expect_equal(ncol(get_manta_sv(projection = "hg38")), 16) #expected number of columns returned with common alternate parameters
  expect_equal(ncol(get_manta_sv(region = "8:128723128-128774067")), 16) #expected number of columns returned with common alternate parameters
  expect_equal(ncol(get_manta_sv(verbose = TRUE)), 16) #expected number of columns returned with common alternate parameters
})


test_that("Comapre returns that should be identical", {
  expect_identical(get_manta_sv(region = "8:128723128-128774067"), #comparing the return based on if the user is using the `region` vs `chromosome`, `qstart` and `qend` parameters
                   get_manta_sv(chromosome = 8, qstart = 128723128, qend = 128774067))
})


test_that("Non-sense examples, expected to fail", {
  expect_error(get_manta_sv(write_to_file = TRUE)) #parameter not should not be available for the GAMBLR.data version of this function
  expect_error(get_manta_sv(from_flatfile = TRUE)) #parameter not should not be available for the GAMBLR.data version of this function
  expect_error(get_manta_sv(from_cache = TRUE)) #parameter not should not be available for the GAMBLR.data version of this function
  expect_error(get_manta_sv(this_is_not_a_parameter = TRUE)) #this parameter does not exist
  expect_error(get_manta_sv(these_sample_ids = "nonexistent_sample")) #this sample ID is made up (shocker!)
  expect_error(get_manta_sv(these_sample_ids = "02-11616_tumorA")) #this sample ID is from a non-supported seq type (capture)
})


test_that("Check classes and types of the returned object", { #check the class of the returned object + the column types
  expect_s3_class(get_manta_sv(), class = "data.frame")
  expect_type(get_manta_sv()[,"CHROM_A"], type = "character")
  expect_type(get_manta_sv()[,"START_A"], type = "double")
  expect_type(get_manta_sv()[,"END_A"], type = "double")
  expect_type(get_manta_sv()[,"CHROM_B"], type = "character")
  expect_type(get_manta_sv()[,"START_B"], type = "double")
  expect_type(get_manta_sv()[,"END_B"], type = "double")
  expect_type(get_manta_sv()[,"manta_name"], type = "character")
  expect_type(get_manta_sv()[,"SCORE"], type = "double")
  expect_type(get_manta_sv()[,"STRAND_A"], type = "character")
  expect_type(get_manta_sv()[,"STRAND_B"], type = "character")
  expect_type(get_manta_sv()[,"tumour_sample_id"], type = "character")
  expect_type(get_manta_sv()[,"normal_sample_id"], type = "character")
  expect_type(get_manta_sv()[,"VAF_tumour"], type = "double")
  expect_type(get_manta_sv()[,"DP"], type = "double")
  expect_type(get_manta_sv()[,"pair_status"], type = "character")
  expect_type(get_manta_sv()[,"FILTER"], type = "character")
})


test_that("Check if filtering options do what's advertised", {
  expect_gte(min(get_manta_sv(min_vaf = 0.5)[,"VAF_tumour"]), 0.5) #do the filtering options work as advertised
  expect_gte(min(get_manta_sv(min_score = 50)[,"SCORE"]), 50) #do the filtering options work as advertised
  expect_true(all(get_manta_sv(pass = TRUE)[,"FILTER"] == "PASS")) #do the filtering options work as advertised
  expect_true(all(get_manta_sv(pairing_status = "matched")[,"pair_status"] == "matched")) #do the filtering options work as advertised
  expect_true(all(get_manta_sv(pairing_status = "unmatched")[,"pair_status"] == "unmatched")) #do the filtering options work as advertised
})


test_that("Check for column existance", { #do we get all the columns we expect
  expect_true(all(c("CHROM_A", "START_A", "END_A", 
                "CHROM_B", "START_B", "END_B", 
                "manta_name", "SCORE", "STRAND_A", 
                "STRAND_B", "tumour_sample_id", 
                "normal_sample_id", "VAF_tumour", 
                "DP", "pair_status", "FILTER") %in% names(get_manta_sv())))  
})


test_that("Chr prefixes only show up where expected", {
  expect_true(all(!grepl("chr", get_manta_sv()[,"CHROM_A"]))) #if grch37 is requested, the chromosomes should not be prefixed
  expect_true(all(!grepl("chr", get_manta_sv()[,"CHROM_B"]))) #if grch37 is requested, the chromosomes should not be prefixed
  expect_true(all(grepl("chr", get_manta_sv(projection = "hg38")[,"CHROM_A"]))) #if hg38 is requested, the chromosomes should be prefixed
  expect_true(all(grepl("chr", get_manta_sv(projection = "hg38")[,"CHROM_B"]))) #if hg38 is requested, the chromosomes should be prefixed
})


test_that("Check if the function correctly subsets to the specified region", { 
  expect_true(all(grepl("8", get_manta_sv(these_sample_ids = "DOHH-2", region = "8:128723128-128774067")[,"CHROM_A"]))) #the return with this region should only have chromosome 8 in the CHROM_A
  expect_gte(min(get_manta_sv(these_sample_ids = "DOHH-2", region = "8:128723128-128774067")[,"START_A"]), 128723128) #check if START_A has the expected coordinates
  expect_lte(min(get_manta_sv(these_sample_ids = "DOHH-2", region = "8:128723128-128774067")[,"END_A"]), 128774067) #check if END_A has the expected coordinates
  expect_true(all(grepl("14", get_manta_sv(these_sample_ids = "DOHH-2", region = "8:128723128-128774067")[,"CHROM_B"]))) #the return with this region should only have chromosome 14 in the CHROM_B
  expect_gte(min(get_manta_sv(these_sample_ids = "DOHH-2", region = "8:128723128-128774067")[,"START_B"]), 106114282) #check if START_B has the expected coordinates
  expect_lte(min(get_manta_sv(these_sample_ids = "DOHH-2", region = "8:128723128-128774067")[,"END_B"]), 106114286) #check if END_B has the expected coordinates
})


test_that("Check the verboseness of the function", {
  expect_message(get_coding_ssm(verbose = TRUE)) #does the verbose option actually output a verbose output
})


test_that("See if no variants are returned when seq type is set to capture (yet, no SVs for capture samples in the bundled data)", {
  expect_equal(nrow(get_manta_sv(this_seq_type = "capture")), 0)
})
