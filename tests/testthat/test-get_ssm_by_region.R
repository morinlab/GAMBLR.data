#load packages
library(testthat)

test_that("Returned rows an columns consitency check", {
  expect_equal(nrow(get_ssm_by_region(region = "chr8:128,723,128-128,774,067")), 71) #grch37
  expect_equal(nrow(get_ssm_by_region(region = "chr5:1000-1000000", projection = "hg38")), 79) #hg38
  expect_equal(ncol(get_ssm_by_region(region = "chr8:128,723,128-128,774,067")), 46) #grch37
  expect_equal(ncol(get_ssm_by_region(region = "chr5:1000-1000000", projection = "hg38")), 46) #hg38
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
