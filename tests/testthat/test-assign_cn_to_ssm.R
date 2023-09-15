#load packages
library(testthat)

test_that("Is the function returning a list (as avertised)", {
  expect_true(is.list(assign_cn_to_ssm(this_sample_id = "DOHH-2")))
})

test_that("Are the elements of the list of the expected data types", {
  cn_ssm = assign_cn_to_ssm(this_sample_id = "DOHH-2")
  expect_true(is.data.frame(cn_ssm$maf))
  expect_true(is.data.frame(cn_ssm$seg))
})

test_that("Are the returned rows and columns consistent", {
  expect_equal(nrow(assign_cn_to_ssm(this_sample_id = "DOHH-2")[[1]]), 22089)
  expect_equal(nrow(assign_cn_to_ssm(this_sample_id = "DOHH-2", projection = "hg38")[[1]]), 21930)
  expect_equal(nrow(assign_cn_to_ssm(this_sample_id = "DOHH-2", coding_only = TRUE)[[1]]), 174)
  expect_equal(nrow(assign_cn_to_ssm(this_sample_id = "DOHH-2", projection = "hg38", coding_only = TRUE)[[1]]), 176)
  expect_equal(nrow(assign_cn_to_ssm(this_sample_id = "DOHH-2", include_silent = TRUE)[[1]]), 22089)
  expect_equal(nrow(assign_cn_to_ssm(this_sample_id = "DOHH-2", projection = "hg38", include_silent = TRUE)[[1]]), 21930)
  expect_equal(nrow(assign_cn_to_ssm(this_sample_id = "DOHH-2")[[2]]), 151)
  expect_equal(nrow(assign_cn_to_ssm(this_sample_id = "DOHH-2", projection = "hg38")[[2]]), 146)
  expect_equal(nrow(assign_cn_to_ssm(this_sample_id = "DOHH-2", coding_only = TRUE)[[2]]), 151)
  expect_equal(nrow(assign_cn_to_ssm(this_sample_id = "DOHH-2", projection = "hg38", coding_only = TRUE)[[2]]), 146)
  expect_equal(nrow(assign_cn_to_ssm(this_sample_id = "DOHH-2", include_silent = TRUE)[[2]]), 151)
  expect_equal(nrow(assign_cn_to_ssm(this_sample_id = "DOHH-2", projection = "hg38", include_silent = TRUE)[[2]]), 146)
  expect_equal(ncol(assign_cn_to_ssm(this_sample_id = "DOHH-2")[[1]]), 48)
  expect_equal(ncol(assign_cn_to_ssm(this_sample_id = "DOHH-2")[[2]]), 7)
})

#more tests to be added...
