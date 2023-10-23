#load packages
library(testthat)

test_that("Check for rows and column consistencies", {
  expect_equal(nrow(get_ssm_by_samples(this_seq_type = "genome")), 188294)
  expect_equal(nrow(get_ssm_by_samples(this_seq_type = "capture")), 20414)
  expect_equal(nrow(get_ssm_by_samples(these_samples_metadata = get_gambl_metadata(seq_type_filter = "genome"))), 188294)
  expect_equal(nrow(get_ssm_by_samples(these_sample_ids = "DOHH-2")), 22089)
  expect_equal(ncol(get_ssm_by_samples(these_sample_ids = "DOHH-2")), 45)
  expect_equal(nrow(get_ssm_by_samples(these_samples_metadata = get_gambl_metadata() %>% dplyr::filter(sample_id == "DOHH-2"))), 22089)
  expect_equal(nrow(get_ssm_by_samples(these_samples_metadata = get_gambl_metadata(), projection = "hg38")), 200109)
  expect_equal(nrow(get_ssm_by_samples(these_sample_ids = "DOHH-2", projection = "hg38")), 21930)
  expect_equal(ncol(get_ssm_by_samples(these_sample_ids = "DOHH-2", projection = "hg38")), 45)
  expect_equal(nrow(get_ssm_by_samples(these_samples_metadata = get_gambl_metadata() %>% dplyr::filter(sample_id == "DOHH-2"), projection = "hg38")), 21930)
})


test_that("Can we specify what MAF columns we want back?", {
  expect_true(all(c("Variant_Type", "Chromosome") %in% names(get_ssm_by_samples(these_sample_ids = "DOHH-2", basic_columns = FALSE, maf_cols = c("Variant_Type", "Chromosome")))))
  expect_equal(ncol(get_ssm_by_samples(these_sample_ids = "DOHH-2", basic_columns = FALSE, maf_cols = c("Variant_Type", "Chromosome"))), 2)
  expect_true(all(c("Variant_Type", "Chromosome") %in% names(get_ssm_by_samples(these_sample_ids = "DOHH-2", projection = "hg38", basic_columns = FALSE, maf_cols = c("Variant_Type", "Chromosome")))))
  expect_equal(ncol(get_ssm_by_samples(these_sample_ids = "DOHH-2", projection = "hg38", basic_columns = FALSE, maf_cols = c("Variant_Type", "Chromosome"))), 2)
})


test_that("Non-sense examples, expected to fail", {
  expect_error(get_ssm_by_samples(these_sample_ids = "DOHH-2", from_flatfile = TRUE)) #parameter not should not be available for the GAMBLR.data version of this function
  expect_error(get_ssm_by_samples(these_sample_ids = "DOHH-2", augmented = TRUE)) #parameter not should not be available for the GAMBLR.data version of this function
  expect_error(get_ssm_by_samples(these_sample_ids = "DOHH-2", this_is_not_a_parameter = TRUE)) #this parameter does not exist
})


test_that("Provide both sample IDs and metadata", {
  expect_identical(get_ssm_by_samples(these_sample_ids = "DOHH-2"),
                   get_ssm_by_samples(these_sample_ids = "DOHH-2", these_samples_metadata = get_gambl_metadata() %>% 
                                         dplyr::filter(sample_id == "DOHH-2")))
  
  expect_identical(get_ssm_by_samples(these_sample_ids = "DOHH-2"),
                   get_ssm_by_samples(these_sample_ids = "DOHH-2", these_samples_metadata = get_gambl_metadata()))
  
  expect_identical(get_ssm_by_samples(these_sample_ids = "DOHH-2"),
                   get_ssm_by_samples(these_samples_metadata = get_gambl_metadata() %>% 
                                         dplyr::filter(sample_id == "DOHH-2")))
})

test_that("Do the min vaf filters work as advertised?", {
  expect_gte(min(get_ssm_by_samples(these_sample_ids = "DOHH-2", min_read_support = 100)[,"t_alt_count"]), 100)
})


test_that("Check the verboseness of the function", {
  expect_message(get_ssm_by_sample(these_sample_ids = "DOHH-2", verbose = TRUE)) #does the verbose option actually output a verbose output
})


test_that("Check if `(this_seq_type = capture)` is working as inteded, i.e only capture samples are returned", {
  #get samples for each seq_type
  cap_samples = unique(get_gambl_metadata(seq_type_filter = "capture") %>% 
                         pull(sample_id))
  
  gen_samples = unique(get_gambl_metadata(seq_type_filter = "genome") %>% 
                         pull(sample_id))
  
  #request capture samples for tested get function
  cap_ssm = unique(get_ssm_by_samples(this_seq_type = "capture") %>% 
                     pull(Tumor_Sample_Barcode))
  
  #are all the requested samples in fact capture samples?
  expect_true(all(cap_ssm %in% cap_samples))
  
  #are the same sample set found in the genome sample pool?
  expect_false(all(cap_ssm %in% gen_samples))
})
