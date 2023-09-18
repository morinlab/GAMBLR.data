#load package
library(testthat)

test_that("Rows and column consistencies", {
  expect_equal(nrow(id_ease(these_sample_ids = "DOHH-2")), 1) #number of rows should be 1, 1 sample ID is requested
  expect_equal(nrow(id_ease()), 4785) #expected number of rows for (all) genome and capture samples
  expect_equal(nrow(id_ease(this_seq_type = "genome")), 1685) #expected number of rows for (all) genome samples
  expect_equal(nrow(id_ease(this_seq_type = "capture")), 3100) #expected number of rows for (all) capture samples
})


test_that("Matching returns", {
  expect_identical(id_ease(this_seq_type = "genome"), 
                   id_ease(these_samples_metadata = get_gambl_metadata(seq_type_filter = "genome"))) #compare return from id_ease with seq type specified to id_ease using these_samples_metadata parameter with the same seq type specified
  
  expect_identical(id_ease(this_seq_type = "capture"), 
                   id_ease(these_samples_metadata = get_gambl_metadata(seq_type_filter = "capture"))) #compare return from id_ease with seq type specified to id_ease using these_samples_metadata parameter with the same seq type specified
  
  expect_identical(id_ease(this_seq_type = c("capture", "genome")), 
                   id_ease(these_samples_metadata = get_gambl_metadata(seq_type_filter = c("capture", "genome")))) #compare return from id_ease with seq type specified to id_ease using these_samples_metadata parameter with the same seq type specified
  
  expect_identical(id_ease(these_sample_ids = "DOHH-2"), 
                   id_ease(these_samples_metadata = get_gambl_metadata() %>% 
                             dplyr::filter(sample_id == "DOHH-2"))) #test if these_sample_ids are returning identical output to these_samples_metadata (same sample ID)

  expect_identical(id_ease(these_sample_ids = c("DOHH-2", "HT", "Toledo", "MD903")), 
                   id_ease(these_samples_metadata = get_gambl_metadata() %>% 
                             dplyr::filter(sample_id %in% c("DOHH-2", "HT", "Toledo", "MD903")))) #test with multiple sample IDs
})


test_that("Handle non-sense", {
  expect_equal(nrow(id_ease(these_samples_metadata = get_gambl_metadata(), 
                            verbose = FALSE,
                            these_sample_ids = c("r2d2","c3P0","Luke","Reddy_3812T"))), 1)
})


test_that("Check the verboseness of the function", {
  expect_message(id_ease(verbose = TRUE)) #does the verbose option actually output a verbose output
})


test_that("Tests meant to fail", {
  expect_error(id_ease(this_is_not_a_parameter = TRUE)) #call a non-existing parameter
  expect_error(id_ease(these_sample_ids = data.frame(sample_id = c("DOHH-2", "HT", "MD903")))) #try to give the function a data frame with sample IDs
})
