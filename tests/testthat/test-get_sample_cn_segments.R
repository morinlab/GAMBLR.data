#load packages
library(testthat)

test_that("Check for rows and column consistencies", {
  expect_equal(nrow(get_sample_cn_segments()), 94244)
  expect_equal(nrow(get_sample_cn_segments(these_samples_metadata = get_gambl_metadata())), 94244)
  expect_equal(ncol(get_sample_cn_segments()), 7)
  expect_equal(nrow(get_sample_cn_segments(these_samples_metadata = get_gambl_metadata() %>% dplyr::filter(sample_id == "DOHH-2"))), 151)
  expect_equal(nrow(get_sample_cn_segments(projection = "hg38")), 33867)
  expect_equal(nrow(get_sample_cn_segments(projection = "hg38", these_samples_metadata = get_gambl_metadata())), 33867)
  expect_equal(ncol(get_sample_cn_segments(projection = "hg38")), 7)
  expect_equal(nrow(get_sample_cn_segments(projection = "hg38", these_samples_metadata = get_gambl_metadata() %>% dplyr::filter(sample_id == "DOHH-2"))), 146)
})


test_that("Check the verboseness of the function", {
  expect_message(get_sample_cn_segments(verbose = TRUE)) #does the verbose option actually output a verbose output
})


test_that("Request capture samples, should fail since this is currently not supported in the bundled data", {
  expect_error(get_sample_cn_segments(seq_type = "capture"))
})


test_that("Is the streamlined option only returning the two expected columns", {
  expect_equal(ncol(get_sample_cn_segments(streamlined = TRUE)), 2)
  expect_true(all(c("ID", "CN") %in% names(get_sample_cn_segments(streamlined = TRUE))))
  expect_equal(ncol(get_sample_cn_segments(streamlined = TRUE, projection = "hg38")), 2)
  expect_true(all(c("ID", "CN") %in% names(get_sample_cn_segments(streamlined = TRUE, projection = "hg38"))))  
})


test_that("Non-sense examples, expected to fail", {
  expect_error(get_sample_cn_segments(from_flatfile = TRUE)) #parameter not should not be available for the GAMBLR.data version of this function
  expect_error(get_sample_cn_segments(augmented = TRUE)) #parameter not should not be available for the GAMBLR.data version of this function
  expect_error(get_sample_cn_segments(this_is_not_a_parameter = TRUE)) #this parameter does not exist
})


test_that("Chr prefixes only show up where expected", {
  expect_true(all(grepl("chr", get_sample_cn_segments(with_chr_prefix = TRUE)[,"chrom"]))) #using grch37 (default), should be easy
  expect_true(all(grepl("chr", get_sample_cn_segments(projection = "hg38", with_chr_prefix = TRUE)[,"chrom"]))) #what if we request projection that is already chr prefixed, will we get double chr prefixes?
  expect_true(all(!grepl("chr", get_sample_cn_segments(with_chr_prefix = FALSE)[,"chrom"]))) #request a projection that is not prefixed, and keep with_chr_prefixes = FALSE
  expect_true(all(!grepl("chr", get_sample_cn_segments(projection = "hg38", with_chr_prefix = FALSE)[,"chrom"]))) #remove chr prefixes from hg38 return
})


test_that("Request non-existing sample IDs", {
  expect_error(get_sample_cn_segments(these_sample_ids = "nonexistent_sample")) #this sample ID is made up (shocker!)
  expect_error(get_sample_cn_segments(these_sample_ids = "02-11616_tumorA")) #this sample ID is from a non-supported seq type (capture)
})


test_that("Provide both sample IDs and metadata", {
  expect_identical(get_sample_cn_segments(these_sample_ids = "DOHH-2"),
                   get_sample_cn_segments(these_sample_ids = "DOHH-2", these_samples_metadata = get_gambl_metadata() %>% 
                                            dplyr::filter(sample_id == "DOHH-2")))
  
  expect_identical(get_sample_cn_segments(these_sample_ids = "DOHH-2"),
                   get_sample_cn_segments(these_sample_ids = "DOHH-2", these_samples_metadata = get_gambl_metadata()))
  
  expect_identical(get_sample_cn_segments(these_sample_ids = "DOHH-2"),
                   get_sample_cn_segments(these_samples_metadata = get_gambl_metadata() %>% 
                                            dplyr::filter(sample_id == "DOHH-2")))
})

test_that("See if no variants are returned when seq type is set to capture (yet, no segments for capture samples in the bundled data)", {
  expect_equal(nrow(get_sample_cn_segments(this_seq_type = "capture")), 0)
})
