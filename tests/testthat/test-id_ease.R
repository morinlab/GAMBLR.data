test_that("sane subsetting", {
  expect_length(id_ease(these_samples_metadata=GAMBLR.data::sample_data$meta,
                       these_sample_ids = slice_sample(GAMBLR.data::sample_data$meta,n=10) %>% pull(sample_id))$these_samples, 10)
  some_sample_ids = slice_sample(GAMBLR.data::sample_data$meta,n=10) %>% pull(sample_id)
  #ensure we get back the same sample_id (and no extra), as intended
  expect_equal(id_ease(these_samples_metadata=GAMBLR.data::sample_data$meta,
                       these_sample_ids = some_sample_ids)$these_samples,some_sample_ids)
  some_sample_ids = slice_sample(GAMBLR.data::sample_data$meta,n=15) %>% pull(sample_id)
  some_sample_ids = c(some_sample_ids,"bonus_feature")
  #this test currently fails. Why is this a warning instead of an error? 
  expect_equal(nrow(id_ease(these_samples_metadata=GAMBLR.data::sample_data$meta,
                       these_sample_ids = some_sample_ids)$this_metadata),length(some_sample_ids))
})



test_that("handles nonsense", {
  #this test currently fails. Why is this a warning instead of an error? 
  expect_length(id_ease(these_samples_metadata=GAMBLR.data::sample_data$meta,
                        these_sample_ids = c("r2d2","c3P0","Luke","Reddy_3812T"))$these_samples, 1)
 
})
