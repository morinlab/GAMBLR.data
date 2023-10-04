#load packages
library(testthat)

test_that("Check for rows and column consistencies", {
  expect_equal(nrow(get_ssm_by_patients(these_samples_metadata = get_gambl_metadata())), 208708)
  expect_equal(nrow(get_ssm_by_patients(these_patient_ids = "DOHH-2")), 22089)
  expect_equal(ncol(get_ssm_by_patients(these_patient_ids = "DOHH-2")), 45)
  expect_equal(nrow(get_ssm_by_patients(these_samples_metadata = get_gambl_metadata() %>% dplyr::filter(sample_id == "DOHH-2"))), 22089)
  expect_equal(nrow(get_ssm_by_patients(these_samples_metadata = get_gambl_metadata(), projection = "hg38")), 200109)
  expect_equal(nrow(get_ssm_by_patients(these_patient_ids = "DOHH-2", projection = "hg38")), 21930)
  expect_equal(ncol(get_ssm_by_patients(these_patient_ids = "DOHH-2", projection = "hg38")), 45)
  expect_equal(nrow(get_ssm_by_patients(these_samples_metadata = get_gambl_metadata() %>% dplyr::filter(sample_id == "DOHH-2"), projection = "hg38")), 21930)
})


test_that("Can we specify what MAF columns we want back?", {
  expect_true(all(c("Variant_Type", "Chromosome") %in% names(get_ssm_by_patients(these_patient_ids = "DOHH-2", basic_columns = FALSE, maf_cols = c("Variant_Type", "Chromosome")))))
  expect_equal(ncol(get_ssm_by_patients(these_patient_ids = "DOHH-2", basic_columns = FALSE, maf_cols = c("Variant_Type", "Chromosome"))), 2)
  expect_true(all(c("Variant_Type", "Chromosome") %in% names(get_ssm_by_patients(these_patient_ids = "DOHH-2", projection = "hg38", basic_columns = FALSE, maf_cols = c("Variant_Type", "Chromosome")))))
  expect_equal(ncol(get_ssm_by_patients(these_patient_ids = "DOHH-2", projection = "hg38", basic_columns = FALSE, maf_cols = c("Variant_Type", "Chromosome"))), 2)
})


test_that("Non-sense examples, expected to fail", {
  expect_error(get_ssm_by_patients(these_patient_ids = "DOHH-2", from_flatfile = TRUE)) #parameter not should not be available for the GAMBLR.data version of this function
  expect_error(get_ssm_by_patients(these_patient_ids = "DOHH-2", augmented = TRUE)) #parameter not should not be available for the GAMBLR.data version of this function
  expect_error(get_ssm_by_patients(these_patient_ids = "DOHH-2", this_is_not_a_parameter = TRUE)) #this parameter does not exist
})


test_that("Provide both sample IDs and metadata", {
  expect_identical(get_ssm_by_patients(these_patient_ids = "DOHH-2"),
                   get_ssm_by_patients(these_patient_ids = "DOHH-2", these_samples_metadata = get_gambl_metadata() %>% 
                                            dplyr::filter(sample_id == "DOHH-2")))
  
  expect_identical(get_ssm_by_patients(these_patient_ids = "DOHH-2"),
                   get_ssm_by_patients(these_patient_ids = "DOHH-2", these_samples_metadata = get_gambl_metadata()))
  
  expect_identical(get_ssm_by_patients(these_patient_ids = "DOHH-2"),
                   get_ssm_by_patients(these_samples_metadata = get_gambl_metadata() %>% 
                                            dplyr::filter(sample_id == "DOHH-2")))
})

test_that("Do the min vaf filters work as advertised?", {
  expect_gte(min(get_ssm_by_patients(these_patient_ids = "DOHH-2", min_read_support = 100)[,"t_alt_count"]), 100)
})


test_that("Request capture samples, should fail since this is currently not supported in the bundled data", {
  expect_error(get_ssm_by_patients(these_patient_ids = "DOHH-2", seq_type = "capture"))
})


test_that("Check the verboseness of the function", {
  expect_message(get_ssm_by_patients(these_patient_ids = "DOHH-2", verbose = TRUE)) #does the verbose option actually output a verbose output
})
