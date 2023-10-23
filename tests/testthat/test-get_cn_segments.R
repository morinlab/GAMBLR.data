#load packages
library(testthat)

test_that("Returned rows an columns consitency check", {
  expect_equal(nrow(get_cn_segments(region = "chr8:128,723,128-128,774,067")), 446) #grch37
  expect_equal(nrow(get_cn_segments(region = "chr8:128,723,128-128,774,067", projection = "hg38")), 282) #hg38
  expect_equal(ncol(get_cn_segments(region = "chr8:128,723,128-128,774,067")), 7) #grch37
  expect_equal(ncol(get_cn_segments(region = "chr8:128,723,128-128,774,067", projection = "hg38")), 7) #hg38
})


test_that("Compare the return for region parameter vs chromosome, qstart and qend + different region formats", {
  expect_identical(get_cn_segments(chromosome = "chr8", qstart = 128723128, qend = 128774067),
                   get_cn_segments(region = "chr8:128,723,128-128,774,067"),
                   get_cn_segments(region = "chr8:128723128-128774067"),
                   get_cn_segments(region = "8:128723128-128774067"))
  
  expect_identical(get_cn_segments(chromosome = "chr8", qstart = 128723128, qend = 128774067, projection = "hg38"),
                   get_cn_segments(region = "chr8:128,723,128-128,774,067", projection = "hg38"),
                   get_cn_segments(region = "chr8:128723128-128774067", projection = "hg38"),
                   get_cn_segments(region = "8:128723128-128774067", projection = "hg38"))
})


test_that("Chr prefixes only show up where expected", {
  expect_true(all(grepl("chr", get_cn_segments(region = "chr8:128,723,128-128,774,067", with_chr_prefix = TRUE)[,"chrom"]))) #using grch37 (default), should be easy
  expect_true(all(grepl("chr", get_cn_segments(region = "chr8:128,723,128-128,774,067", projection = "hg38", with_chr_prefix = TRUE)[,"chrom"]))) #what if we request projection that is already chr prefixed, will we get double chr prefixes?
  expect_true(all(!grepl("chr", get_cn_segments(region = "chr8:128,723,128-128,774,067", with_chr_prefix = FALSE)[,"chrom"]))) #request a projection that is not prefixed, and keep with_chr_prefixes = FALSE
  expect_true(all(!grepl("chr", get_cn_segments(region = "chr8:128,723,128-128,774,067", projection = "hg38", with_chr_prefix = FALSE)[,"chrom"]))) #remove chr prefixes from hg38 return
})


test_that("Try to give the function multiple regions", {
  expect_error(get_cn_segments(region = c("chr8:128,723,128-128,774,067", "chr1:1000-1000000"))) #give the function multiple regions with the region parameter
  expect_error(get_cn_segments(chromosome = c("chr1", "chr8"), qstart = c("10000", "10000"), qend = c("1000000", "1000000"))) #give the function multiple regions with chromosome, qstart and qend parameters.
  expect_error(get_cn_segments(region = c("chr8:128,723,128-128,774,067", "chr1:1000-1000000"), projection = "hg38")) #give the function multiple regions with the region parameter
  expect_error(get_cn_segments(chromosome = c("chr1", "chr8"), qstart = c("10000", "10000"), qend = c("1000000", "1000000"), projection = "hg38")) #give the function multiple regions with chromosome, qstart and qend parameters.
})


test_that("Try non-existing parameters", {
  expect_error(get_cn_segments(region = "8:128723128-128774067", from_flatfile = TRUE)) #this parameter does not exists in this version of the function
  expect_error(get_cn_segments(region = "8:128723128-128774067", from_flatfile = TRUE, projection = "hg38"))
})


test_that("Is the streamlined option only returning the two expected columns", {
  expect_equal(ncol(get_cn_segments(region = "8:128723128-128774067", streamlined = TRUE)), 2)
  expect_true(all(c("ID", "CN") %in% names(get_cn_segments(region = "8:128723128-128774067", streamlined = TRUE))))
  expect_equal(ncol(get_cn_segments(region = "8:128723128-128774067", streamlined = TRUE, projection = "hg38")), 2)
  expect_true(all(c("ID", "CN") %in% names(get_cn_segments(region = "8:128723128-128774067", streamlined = TRUE, projection = "hg38"))))  
})

test_that("See if no variants are returned when seq type is set to capture (yet, no segments for capture samples in the bundled data)", {
  expect_equal(nrow(get_cn_segments(region = "8:128723128-128774067", this_seq_type = "capture")), 0)
})


test_that("Does these_sample_ids and these_sample_metadata work as inteded", {
  #get CN for DOHH2 (using these_sample_ids)
  dohh2_cn = unique(get_cn_segments(region = "8:128723128-128774067", these_sample_ids = "DOHH-2") %>% 
                       pull(ID))
  
  #is DOHH-2 the only sample returned
  expect_true(dohh2_cn %in% "DOHH-2")
  
  #get CN for dlbcl cell lines (using these_samples_metadata)
  dlbcl_cn = get_cn_segments(region = "8:128723128-128774067", these_samples_metadata = get_gambl_metadata() %>% 
                                                  dplyr::filter(cohort == "DLBCL_cell_lines")) %>% 
                                  pull(ID)
  
  #are CNs returned for the expected sample IDs?
  expect_true(all(dlbcl_cn %in% c("DOHH-2", "OCI-Ly10", "OCI-Ly3", "SU-DHL-10", "SU-DHL-4")))
})
