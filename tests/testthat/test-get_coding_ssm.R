#load packages
library(testthat)

test_that("Check for rows and column consistencies", {
  expect_equal(nrow(get_coding_ssm()), 46559)
  expect_equal(nrow(get_coding_ssm(these_samples_metadata = get_gambl_metadata())), 60392)
  expect_equal(ncol(get_coding_ssm()), 46)
  expect_equal(nrow(get_coding_ssm(these_samples_metadata = get_gambl_metadata() %>% dplyr::filter(sample_id == "DOHH-2"))), 234)
  expect_equal(nrow(get_coding_ssm(projection = "hg38")), 59143)
  expect_equal(nrow(get_coding_ssm(projection = "hg38", these_samples_metadata = get_gambl_metadata())), 59143)
  expect_equal(ncol(get_coding_ssm(projection = "hg38")), 46)
  expect_equal(nrow(get_coding_ssm(projection = "hg38", these_samples_metadata = get_gambl_metadata() %>% dplyr::filter(sample_id == "DOHH-2"))), 234)
})


test_that("Does inlcude_silent work as intended", {
  expect_true(all(grepl("Silent", get_coding_ssm(include_silent = TRUE)[,"Variant_Classification"]))) #include `Variant_Classification == Silent`
  expect_true(all(!grepl("Silent", get_coding_ssm(include_silent = FALSE)[,"Variant_Classification"]))) #exclude `Variant_Classification == Silent`
})


test_that("Check the verboseness of the function", {
  expect_message(get_coding_ssm(verbose = TRUE)) #does the verbose option actually output a verbose output
})


test_that("Can we specify what MAF columns we want back?", {
  expect_true(all(c("Variant_Type", "Chromosome") %in% names(get_coding_ssm(basic_columns = FALSE, maf_cols = c("Variant_Type", "Chromosome")))))
  expect_equal(ncol(get_coding_ssm(basic_columns = FALSE, maf_cols = c("Variant_Type", "Chromosome"))), 2)
  expect_true(all(c("Variant_Type", "Chromosome") %in% names(get_coding_ssm(projection = "hg38", basic_columns = FALSE, maf_cols = c("Variant_Type", "Chromosome")))))
  expect_equal(ncol(get_coding_ssm(projection = "hg38", basic_columns = FALSE, maf_cols = c("Variant_Type", "Chromosome"))), 2)
})


test_that("Non-sense examples, expected to fail", {
  expect_error(get_coding_ssm(engine = "fread_maf")) #parameter not should not be available for the GAMBLR.data version of this function
  expect_error(get_coding_ssm(groups = c("gambl", "icgc_dart"))) #parameter not should not be available for the GAMBLR.data version of this function
  expect_error(get_coding_ssm(from_flatfile = TRUE)) #parameter not should not be available for the GAMBLR.data version of this function
  expect_error(get_coding_ssm(augmented = TRUE)) #parameter not should not be available for the GAMBLR.data version of this function
  expect_error(get_coding_ssm(this_is_not_a_parameter = TRUE)) #this parameter does not exist
})


test_that("Do the min vaf filters work as advertised?", {
  expect_gte(min(get_coding_ssm(min_read_support = 100)[,"t_alt_count"]), 100)
})


test_that("Do the built in metadata subset options work as advertised?", {
  #exclude DLBCL_cell_line samples (cohort)
  dlbcl_cell_lines_samples = get_gambl_metadata() %>% dplyr::filter(cohort == "DLBCL_cell_lines") %>% pull(sample_id)
  expect_false(all(get_coding_ssm(exclude_cohort = "DLBCL_cell_lines")[,"Tumor_Sample_Barcode"] %in% dlbcl_cell_lines_samples))
  
  #only return SSM for DLBCL_cell_line samples (cohort)
  expect_equal(get_coding_ssm(limit_cohort = "DLBCL_cell_lines"),
               get_coding_ssm(these_samples_metadata = get_gambl_metadata() %>% 
                                dplyr::filter(cohort == "DLBCL_cell_lines")))
  
  #exclude DLBCL samples (pathology)
  dlbcl_samples = get_gambl_metadata() %>% dplyr::filter(pathology == "DLBCL") %>% pull(sample_id)
  dlbcl_ssm = get_coding_ssm(limit_pathology = "DLBCL") %>% pull(Tumor_Sample_Barcode)
  expect_true(all(dlbcl_ssm %in% dlbcl_samples))
  
  #use limit_samples to restrict the return to specific sample IDs
  expect_equal(get_coding_ssm(limit_samples = c("DOH22", "OCI-Ly10")), 
               get_coding_ssm(these_samples_metadata = get_gambl_metadata() %>% 
                                dplyr::filter(sample_id %in% c("DOH22", "OCI-Ly10")))) 
  
  #use limit_samples to restrict the return to specific sample IDs
  expect_equal(get_coding_ssm(limit_samples = "DOH22", projection = "hg38"), 
               get_coding_ssm(projection = "hg38",
                              these_samples_metadata = get_gambl_metadata() %>% 
                                dplyr::filter(sample_id == "DOH22"))) 
})


test_that("Force unmatched samples", {
  expect_false(isTRUE(all.equal(get_coding_ssm(force_unmatched_samples = "DOHH-2", 
                                               these_samples_metadata = get_gambl_metadata() %>% 
                                                 dplyr::filter(cohort == "DLBCL_cell_lines")),
                                
               get_coding_ssm(these_samples_metadata = get_gambl_metadata() %>% 
                                dplyr::filter(cohort == "DLBCL_cell_lines")))))
})


test_that("Check if `(this_seq_type = capture)` is working as inteded, i.e only capture samples are returned", {
  #get samples for each seq_type
  cap_samples = unique(get_gambl_metadata(seq_type_filter = "capture") %>% 
                         pull(sample_id))
  
  gen_samples = unique(get_gambl_metadata(seq_type_filter = "genome") %>% 
                         pull(sample_id))
  
  #request capture samples for tested get function
  cap_ssm = unique(get_coding_ssm(this_seq_type = "capture") %>% 
                     pull(Tumor_Sample_Barcode))
  
  #are all the requested samples in fact capture samples?
  expect_true(all(cap_ssm %in% cap_samples))
  
  #are the same sample set found in the genome sample pool?
  expect_false(all(cap_ssm %in% gen_samples))
})
