#load packages
library(testthat)

test_that("Returned rows an columns consitency check", {
  expect_equal(nrow(get_ssm_by_region(region = "chr8:128,723,128-128,774,067")), 71) #grch37
  expect_equal(nrow(get_ssm_by_region(region = "chr5:1000-1000000", projection = "hg38")), 79) #hg38
  expect_equal(ncol(get_ssm_by_region(region = "chr8:128,723,128-128,774,067")), 45) #grch37
  expect_equal(ncol(get_ssm_by_region(region = "chr5:1000-1000000", projection = "hg38")), 45) #hg38
})


test_that("Compare the return for region parameter vs chromosome, qstart and qend + different region formats", {
  expect_identical(get_ssm_by_region(chromosome = "chr8", qstart = 128723128, qend = 128774067),
                   get_ssm_by_region(region = "chr8:128,723,128-128,774,067"),
                   get_ssm_by_region(region = "chr8:128723128-128774067"),
                   get_ssm_by_region(region = "8:128723128-128774067"))
  
  expect_identical(get_ssm_by_region(chromosome = "chr5", qstart = 1000, qend = 1000000, projection = "hg38"),
                   get_ssm_by_region(region = "chr5:100,0-100,000,0", projection = "hg38"),
                   get_ssm_by_region(region = "chr5:1000-1000000", projection = "hg38"),
                   get_ssm_by_region(region = "5:1000-1000000", projection = "hg38"))
})


test_that("Is the streamlined option only returning the two expected columns", {
  expect_equal(ncol(get_ssm_by_region(region = "8:128723128-128774067", streamlined = TRUE)), 2)
  expect_true(all(c("Start_Position", "Tumor_Sample_Barcode") %in% names(get_ssm_by_region(region = "8:128723128-128774067", streamlined = TRUE))))
  expect_equal(ncol(get_ssm_by_region(region = "chr5:1000-1000000", streamlined = TRUE, projection = "hg38")), 2)
  expect_true(all(c("Start_Position", "Tumor_Sample_Barcode") %in% names(get_ssm_by_region(region = "chr5:1000-1000000", streamlined = TRUE, projection = "hg38"))))  
})


test_that("Do the min vaf filters work as advertised?", {
  expect_gte(min(get_ssm_by_region(region = "8:128723128-128774067", min_read_support = 5)[,"t_alt_count"]), 5)
})


test_that("Check the verboseness of the function", {
  expect_message(get_ssm_by_region(region = "8:128723128-128774067", verbose = TRUE)) #does the verbose option actually output a verbose output
})


test_that("Non-sense examples, expected to fail", {
  expect_error(get_ssm_by_region(region = "8:128723128-128774067", from_flatfile = TRUE)) #parameter not should not be available for the GAMBLR.data version of this function
  expect_error(get_ssm_by_region(region = "8:128723128-128774067", augmented = TRUE)) #parameter not should not be available for the GAMBLR.data version of this function
  expect_error(get_ssm_by_region(region = "8:128723128-128774067", this_is_not_a_parameter = TRUE)) #this parameter does not exist
})


test_that("Try to give the function more than one region", {
  expect_error(get_ssm_by_region(region = c("8:128,723,128-128,774,067", "chr8:127736231-127742951")))
})


test_that("Check if `(this_seq_type = capture)` is working as inteded, i.e only capture samples are returned", {
  #get samples for each seq_type
  cap_samples = unique(get_gambl_metadata(seq_type_filter = "capture") %>% 
                         pull(sample_id))
  
  gen_samples = unique(get_gambl_metadata(seq_type_filter = "genome") %>% 
                         pull(sample_id))
  
  #request capture samples for tested get function
  cap_ssm = unique(get_ssm_by_region(region = "chr8:128,723,128-128,774,067", this_seq_type = "capture") %>% 
                     pull(Tumor_Sample_Barcode))
  
  #are all the requested samples in fact capture samples?
  expect_true(all(cap_ssm %in% cap_samples))
  
  #are the same sample set found in the genome sample pool?
  expect_false(all(cap_ssm %in% gen_samples))
})


test_that("Does these_sample_ids and these_sample_metadata work as inteded", {
  #get ssm for DOHH2 (using these_sample_ids)
  dohh2_ssm = unique(get_ssm_by_region(region = "chr8:128,723,128-128,774,067", these_sample_ids = "DOHH-2") %>% 
                       pull(Tumor_Sample_Barcode))
  
  #is DOHH-2 the only sample returned
  expect_true(dohh2_ssm %in% "DOHH-2")
  
  #get ssm for dlbcl cell lines (using these_samples_metadata)
  dlbcl_ssm = get_ssm_by_region(region = "chr8:128,723,128-128,774,067", these_samples_metadata = get_gambl_metadata() %>% dplyr::filter(cohort == "DLBCL_cell_lines")) %>% pull(Tumor_Sample_Barcode)
  
  #are ssm returned for the expected sample IDs?
  expect_true(all(dlbcl_ssm %in% c("DOHH-2", "OCI-Ly10", "OCI-Ly3", "SU-DHL-10", "SU-DHL-4")))
})
