


test_that("consistent row number", {
  expect_equal(nrow(get_gambl_metadata(seq_type_filter=c("genome","capture"))), 4785)
  expect_equal(nrow(get_gambl_metadata(seq_type_filter=c("capture"))), 3100)
  expect_equal(nrow(get_gambl_metadata(seq_type_filter=c("genome"))), 4785-3100)
})


test_that("all bundled samples in metadata", {
  expect_false(any(!unique(c(GAMBLR.data::sample_data$grch37$maf$Tumor_Sample_Barcode,
                       GAMBLR.data::sample_data$grch37$seg$ID,
                       GAMBLR.data::sample_data$grch37$bedpe$tumour_sample_id)) %in% GAMBLR.data::sample_data$meta$sample_id)
  )
  expect_false(any(!unique(c(GAMBLR.data::sample_data$hg38$maf$Tumor_Sample_Barcode,
                             GAMBLR.data::sample_data$hg38$seg$ID,
                             GAMBLR.data::sample_data$hg38$bedpe$tumour_sample_id)) %in% GAMBLR.data::sample_data$meta$sample_id)
  )
})
