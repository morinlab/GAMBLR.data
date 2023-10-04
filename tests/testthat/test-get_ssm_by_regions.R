#load packages
library(testthat)

test_that("Returned rows an columns consitency check", {
  expect_equal(nrow(get_ssm_by_regions(regions_bed = dplyr::mutate(GAMBLR.data::grch37_ashm_regions, name = paste(gene, region, sep = "_")))), 9323) #grch37
  expect_equal(nrow(get_ssm_by_regions(regions_bed = dplyr::mutate(GAMBLR.data::hg38_ashm_regions, name = paste(gene, region, sep = "_")), projection = "hg38")), 1515) #hg38
  expect_equal(nrow(get_ssm_by_regions(regions_list = c("chr1:6661482-6662702", "chr2:232572640-232574297", "chr17:75443766-75451177"))), 41) #grch37
  expect_equal(nrow(get_ssm_by_regions(projection = "hg38", regions_list = c("chr1:100-10000000", "chr2:50000-5000000", "chr17:60000-60000000"))), 5050) #hg38
  expect_equal(ncol(get_ssm_by_regions(regions_list = c("chr1:6661482-6662702", "chr2:232572640-232574297", "chr17:75443766-75451177"))), 3) #grch37
  expect_equal(ncol(get_ssm_by_regions(projection = "hg38", regions_list = c("chr1:100-10000000", "chr2:50000-5000000", "chr17:60000-60000000"))), 3) #hg38
})


test_that("Ensure that regions_bed and regions_vector return the same SSM (if the same regions are specified)", {
  expect_identical(get_ssm_by_regions(regions_list = c("chr1:6661482-6662702", 
                                                         "chr2:232572640-232574297", 
                                                         "chr17:75443766-75451177")),
                   get_ssm_by_regions(regions_bed = data.frame(chr = c("chr1", 
                                                                       "chr2", 
                                                                       "chr17"), 
                                                               start = c("6661482", 
                                                                         "232572640", 
                                                                         "75443766"), 
                                                               end = c("6662702", 
                                                                       "232574297", 
                                                                       "75451177"))))
})


test_that("Is the streamlined option only returning the two expected columns", {
  expect_equal(ncol(get_ssm_by_regions(regions_list = c("chr1:6661482-6662702", "chr2:232572640-232574297", "chr17:75443766-75451177"), streamlined = TRUE)), 3)
  expect_true(all(c("start", "sample_id", "region_name") %in% names(get_ssm_by_regions(regions_list = c("chr1:6661482-6662702", "chr2:232572640-232574297", "chr17:75443766-75451177"), streamlined = TRUE))))
  expect_equal(ncol(get_ssm_by_regions(regions_list = c("chr1:6661482-6662702", "chr2:232572640-232574297", "chr17:75443766-75451177"), streamlined = TRUE, projection = "hg38")), 3)
  expect_true(all(c("start", "sample_id", "region_name") %in% names(get_ssm_by_regions(regions_list = c("chr1:6661482-6662702", "chr2:232572640-232574297", "chr17:75443766-75451177"), streamlined = TRUE, projection = "hg38"))))  
})


test_that("See if the use_name_column parameter work as advertised", {
  regions_df = data.frame(chr = c("chr1", "chr2", "chr17"), 
                          start = c("6661482", "232572640", "75443766"), 
                          end = c("6662702", "232574297", "75451177"),
                          name = c("one_region", "another_region", "my_last_region"))
  
  region_ssm = get_ssm_by_regions(regions_bed = regions_df, use_name_column = TRUE)
  
  expect_true(all(c("one_region", "another_region", "my_last_region") %in% region_ssm$region_name))
})


test_that("Check the verboseness of the function", {
  regions_df = data.frame(chr = c("chr1", "chr2", "chr17"), 
                          start = c("6661482", "232572640", "75443766"), 
                          end = c("6662702", "232574297", "75451177"),
                          name = c("one_region", "another_region", "my_last_region"))
  
  expect_message(get_ssm_by_regions(regions_bed = regions_df, verbose = TRUE)) #does the verbose option actually output a verbose output
})


test_that("Do the min vaf filters work as advertised?", {
  regions_df = data.frame(chr = c("chr1", "chr2", "chr17"), 
                          start = c("6661482", "232572640", "75443766"), 
                          end = c("6662702", "232574297", "75451177"),
                          name = c("one_region", "another_region", "my_last_region"))
  expect_gte(min(get_ssm_by_regions(regions_bed = regions_df, min_read_support = 5, streamlined = FALSE)[,"t_alt_count"]), 5)
})


test_that("Is the streamlined option only returning the expected columns", {
  regions_df = data.frame(chr = c("chr1", "chr2", "chr17"), 
                          start = c("6661482", "232572640", "75443766"), 
                          end = c("6662702", "232574297", "75451177"),
                          name = c("one_region", "another_region", "my_last_region"))
  
  expect_equal(ncol(get_ssm_by_regions(regions_bed = regions_df, streamlined = FALSE)), 45)
  expect_equal(ncol(get_ssm_by_regions(regions_bed = regions_df)), 3)
  expect_equal(ncol(get_ssm_by_regions(projection = "hg38", regions_bed = regions_df, streamlined = FALSE)), 45)
  expect_equal(ncol(get_ssm_by_regions(projection = "hg38", regions_bed = regions_df)), 3)
})


test_that("Non-sense examples, expected to fail", {
  regions_df = data.frame(chr = c("chr1", "chr2", "chr17"), 
                          start = c("6661482", "232572640", "75443766"), 
                          end = c("6662702", "232574297", "75451177"),
                          name = c("one_region", "another_region", "my_last_region"))
  
  expect_error(get_ssm_by_regions(regions_bed = regions_df, from_flatfile = TRUE)) #parameter not should not be available for the GAMBLR.data version of this function
  expect_error(get_ssm_by_regions(regions_bed = regions_df, this_is_not_a_parameter = TRUE)) #this parameter does not exist
})


test_that("Check if `(this_seq_type = capture)` is working as inteded, i.e only capture samples are returned", {
  #get samples for each seq_type
  cap_samples = unique(get_gambl_metadata(seq_type_filter = "capture") %>% 
                         pull(sample_id))
  
  gen_samples = unique(get_gambl_metadata(seq_type_filter = "genome") %>% 
                         pull(sample_id))
  
  regions_df = data.frame(chr = c("1", "6"), 
                          start = c("2488019", "47125723"), 
                          end = c("2488105", "187447758"),
                          name = c("one_region", "another_region"))
  
  #request capture samples for tested get function
  cap_ssm = unique(get_ssm_by_regions(regions_bed = regions_df, this_seq_type = "capture") %>% 
                     pull(sample_id))
  
  #are all the requested samples in fact capture samples?
  expect_true(all(cap_ssm %in% cap_samples))
  
  #are the same sample set found in the genome sample pool?
  expect_false(all(cap_ssm %in% gen_samples))
})


test_that("Does these_sample_ids and these_sample_metadata work as inteded", {
  #define a region
  regions_df = data.frame(chr = c("1", "6"), 
                          start = c("2488019", "47125723"), 
                          end = c("2488105", "187447758"),
                          name = c("one_region", "another_region"))
  
  #get ssm for DOHH2 (using these_sample_ids)
  dohh2_ssm = unique(get_ssm_by_regions(regions_bed = regions_df, these_sample_ids = "DOHH-2") %>% 
                       pull(sample_id))
  
  #is DOHH-2 the only sample returned
  expect_true(dohh2_ssm %in% "DOHH-2")
  
  #get ssm for dlbcl cell lines (using these_samples_metadata)
  dlbcl_ssm = get_ssm_by_regions(regions_bed = regions_df, these_samples_metadata = get_gambl_metadata() %>% dplyr::filter(cohort == "DLBCL_cell_lines")) %>% pull(sample_id)
  
  #are ssm returned for the expected sample IDs?
  expect_true(all(dlbcl_ssm %in% c("DOHH-2", "OCI-Ly10", "OCI-Ly3", "SU-DHL-10", "SU-DHL-4")))
})
