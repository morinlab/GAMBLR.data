#load packages
library(testthat)

test_that("Check dimensions of all collated results", {
  expect_equal(nrow(collate_results()), 2376) #get all the collated results
  expect_equal(ncol(collate_results()), 17) #get all the columns back
  expect_equal(ncol(collate_results(join_with_full_metadata = TRUE)), 43) #horizontally expand the metadata with collated results
  expect_identical(nrow(collate_results(join_with_full_metadata = TRUE)), nrow(get_gambl_metadata(seq_type_filter = c("genome", "capture")))) #retain the full metadata and only add in QC metrics for samples that have such information
})


test_that("Check if the sample_table parameter behaves as expected", {
  expect_true(all(grepl("BNHL_10T", collate_results(sample_table = "BNHL_10T")[,"sample_id"]))) #return collate results for one specific sample ID (as string)
  expect_equal(nrow(unique(collate_results(sample_table = c("DLBCL-c_D_1141-Tumor", "DLBCL-RICOVER_478-Tumor")))), 2) #return collate results for one specific sample ID (as a vector of characters)
  expect_equal(nrow(unique(collate_results(sample_table = data.frame(sample_id = c("DLBCL-c_D_1141-Tumor", "DLBCL-RICOVER_478-Tumor"))))), 2) #test if we can supply a data frame with samples IDs to the sample_table parameter
})


test_that("Check these_samples_metadata parameter", {
  expect_true(all(grepl("BNHL_10T", collate_results(these_samples_metadata = get_gambl_metadata() %>% dplyr::filter(sample_id == "BNHL_10T"))[,"sample_id"]))) #use a subset of metadata to retreive collated resutls for
  expect_equal(collate_results(these_samples_metadata = get_gambl_metadata() %>% dplyr::filter(sample_id == "BNHL_10T")), #sample_table and these_samples_metadata should behave the same
               collate_results(sample_table = "BNHL_10T"))
})


test_that("Check join_with_full_metadata parameter", {
  expect_equal(nrow(collate_results(these_samples_metadata = get_gambl_metadata() %>% dplyr::filter(sample_id == "BNHL_10T"), join_with_full_metadata = TRUE)), 1) #get one sample back (the one specified in these_samples_metadata)
  expect_equal(ncol(collate_results(these_samples_metadata = get_gambl_metadata() %>% dplyr::filter(sample_id == "BNHL_10T"), join_with_full_metadata = TRUE)), 43) # get extended metadata columns added for the sample ID in the provided metadata
  expect_equal(nrow(collate_results(sample_table = "BNHL_10T", join_with_full_metadata = TRUE)), 4785) #non-sense parameter combinations (if the user wants to join with metadata, they need top provide a metadata table, or else all samples in the unfiltered metadata table will be used)
  expect_equal(nrow(collate_results(join_with_full_metadata = TRUE)), 4785) #get all metadata back but with added collated results
})


test_that("Check seq_type_filer parameter", {
  expect_equal(nrow(collate_results(join_with_full_metadata = TRUE, seq_type_filter = "genome")), 1685) #link with meta data for all genome samples
  expect_true(all(grepl("genome", collate_results(join_with_full_metadata = TRUE, seq_type_filter = "genome")[,"seq_type"]))) #get genome samples
  expect_equal(nrow(collate_results(join_with_full_metadata = TRUE, seq_type_filter = "capture")), 3100) #link with meta data for all capture samples
  expect_true(all(grepl("capture", collate_results(join_with_full_metadata = TRUE, seq_type_filter = "capture")[,"seq_type"]))) #get capture samples
  expect_equal(nrow(collate_results(join_with_full_metadata = TRUE, seq_type_filter = "capture")), 3100) #link with meta data for all genome/capture samples
  all_collated = collate_results(join_with_full_metadata = TRUE, seq_type_filter = c("genome", "capture")) #get genome and capture
  expect_true(all(c("genome", "capture") %in% unique(all_collated$seq_type))) #genome and capture as the only two seq types
})


test_that("Non-sense parameter combinations", {
  expect_false(all(grepl("capture", collate_results(these_samples_metadata = get_gambl_metadata(seq_type_filter = "genome"), seq_type_filter = "capture", join_with_full_metadata = TRUE)[, "seq_type"]))) #if provided these_samples_metadata, the seq type in the metadata should take precedence over seq_type_filter
})


test_that("Non-existing parameters", {
  expect_error(collate_results(from_cache = TRUE)) #try a parameter from the GAMBLR.results version of this function
  expect_error(collate_results(this_is_not_a_parameter = TRUE)) #try a non-existing parameter
})
