#load packages
library(testthat)

test_that("Check the type of the return", {
  expect_true(is.list(assign_cn_to_ssm(this_sample_id = "DOHH-2"))) #is the return in fact a list?
})


test_that("Check data types of the elements in the returned list", {
  expect_true(is.data.frame(assign_cn_to_ssm(this_sample_id = "DOHH-2")[['maf']])) #is the maf object a data frame?
  expect_true(is.data.frame(assign_cn_to_ssm(this_sample_id = "DOHH-2")[['seg']])) #is the seg object a data frame?
})


test_that("Row and column consitencies", {
  expect_equal(nrow(assign_cn_to_ssm(this_sample_id = "DOHH-2")[['maf']]), 22089) #number of SSM returned for DOHH-2 (grch37)
  expect_equal(nrow(assign_cn_to_ssm(this_sample_id = "DOHH-2", projection = "hg38")[['maf']]), 21930) #number of SSM returned for DOHH-2 (hg38)
  expect_equal(nrow(assign_cn_to_ssm(this_sample_id = "DOHH-2", coding_only = TRUE)[['maf']]), 174) #number of SSM returned for DOHH-2 (grch37, coding only)
  expect_equal(nrow(assign_cn_to_ssm(this_sample_id = "DOHH-2", projection = "hg38", coding_only = TRUE)[['maf']]), 176) #number of SSM returned for DOHH-2 (hg38, coding only)
  expect_equal(nrow(assign_cn_to_ssm(this_sample_id = "DOHH-2", include_silent = TRUE)[['maf']]), 22089) #number of SSM returned for DOHH-2 (grch37, include silent)
  expect_equal(nrow(assign_cn_to_ssm(this_sample_id = "DOHH-2", projection = "hg38", include_silent = TRUE)[['maf']]), 21930) #number of SSM returned for DOHH-2 (hg38, include silent)
  expect_equal(nrow(assign_cn_to_ssm(this_sample_id = "DOHH-2")[['seg']]), 151) #number of segs returned for DOHH-2 (grch37)
  expect_equal(nrow(assign_cn_to_ssm(this_sample_id = "DOHH-2", projection = "hg38")[['seg']]), 146) #number of segs returned for DOHH-2 (hg38)
  expect_equal(ncol(assign_cn_to_ssm(this_sample_id = "DOHH-2")[['maf']]), 48) #check number of columns in the returned MAF (grch37)
  expect_equal(ncol(assign_cn_to_ssm(this_sample_id = "DOHH-2")[['seg']]), 7) #check number of columns in the returned SEG (grch37)
})


test_that("assume_diploid parameter", {
  expect_true(all(grepl(2, assign_cn_to_ssm(this_sample_id = "DOHH-2", assume_diploid = TRUE)[['maf']][,"CN"]))) #Check if in fact all returned CN states are diploid (grch37)
  expect_true(all(grepl(2, assign_cn_to_ssm(this_sample_id = "DOHH-2", projection = "hg38", assume_diploid = TRUE)[['maf']][,"CN"]))) #Check if in fact all returned CN states are diploid (hg38)
})


test_that("include_silent parameter", {
  expect_true(all(grepl("Silent", assign_cn_to_ssm(this_sample_id = "DOHH-2", include_silent = TRUE)[['maf']][,"Variant_Classification"]))) #Check if we get silent mutations returned if parameter set to TRUE
  expect_true(all(!grepl("Silent", assign_cn_to_ssm(this_sample_id = "DOHH-2", include_silent = FALSE, coding_only = TRUE)[['maf']][,"Variant_Classification"]))) #Check if silent mutations are included if we set this parameter to FALSE
})


test_that("gene parameter", {
  expect_error(assign_cn_to_ssm(genes = "this_is_not_a_gene")) #Give the function a non-existing gene
  expect_true(all(grepl("MYC", assign_cn_to_ssm(this_sample_id = "DOHH-2", genes = "MYC")[['maf']][,"Hugo_Symbol"]))) #Only return variants within MYC (grch37)
  expect_equal(nrow(assign_cn_to_ssm(genes = "MYC", this_sample_id = "DOHH-2")[['maf']]), 2) #check if the number of rows matches the expected value
  expect_true(all(grepl("MYC", assign_cn_to_ssm(this_sample_id = "DOHH-2", genes = "MYC", projection = "hg38")[['maf']][,"Hugo_Symbol"]))) #Only return variants within MYC (hg38)
  expect_true(all(assign_cn_to_ssm(this_sample_id = "DOHH-2", genes = c("MYC", "BCL2"))[['maf']][, Hugo_Symbol] %in% c("MYC", "BCL2"))) #Check if the function can return variants for multiple genes
})


test_that("Tests meant to fail", {
  expect_error(assign_cn_to_ssm(this_sample_id = "DOHH-2", this_is_not_a_paraemter = TRUE)) # This is not a parameter
  expect_error(assign_cn_to_ssm(this_sample_id = c("DOHH-2", "HT", "Toledo", "MD903"))) #Provide more than one sample ID
  expect_error(assign_cn_to_ssm(projection = "HG38", this_sample_id = "DOHH-2")) #provide a non-valid projection
  expect_error(assign_cn_to_ssm(this_seq_type = "capture", this_sample_id = "DOHH-2")) #Try to get capture calls back
})


test_that("Only coding mutations are kept", {
  expect_true(all(assign_cn_to_ssm(this_sample_id = "DOHH-2", coding_only = TRUE)[['maf']][, Variant_Classification] %in% coding_class)) #Check if all returned variants has the variant classification dictated by the coding_class object.
})

