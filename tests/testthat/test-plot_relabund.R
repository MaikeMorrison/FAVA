
relab_matrix = xue_microbiome_sample
group = "subject"
time = "timepoint"
K = NULL
w = NULL
arrange = FALSE

w = time_weights(times = relab_matrix$timepoint, group = relab_matrix$subject)


# Test the weighting function
test_that("relab_sample_weighter works", {

  # Many groups, one matrix ----------------------------------------------------
  # no, K, yes group and time
  expect_no_error(relab_sample_weighter(relab = relab_matrix, group = "subject", time = "timepoint"))
  # K and time, no group
  expect_error(relab_sample_weighter(relab = relab_matrix, K = 524, time = "timepoint"))
  # K and group and time
  expect_no_error(relab_sample_weighter(relab = relab_matrix, K = 524, group = "subject", time = "timepoint"))

  # Specify w instead of time ---------------------------------
  expect_no_error(relab_sample_weighter(relab = relab_matrix, group = "subject", w = w, K = 524))

  # One group ------------------------------------------------------------------
  # K and group and time
  expect_no_error(relab_sample_weighter(relab = dplyr::filter(relab_matrix, subject == "XBA"),
                                                              K = 524, group = "subject", time = "timepoint"))
  # group and time, no K
  expect_no_error(relab_sample_weighter(relab = dplyr::filter(relab_matrix, subject == "XBA"),
                                        group = "subject", time = "timepoint"))
  # K, no group
  expect_no_error(relab_sample_weighter(relab = dplyr::filter(relab_matrix, subject == "XBA"),
                                        K = 524, time = "timepoint"))

  # Specify w instead of time ---------------------------------
  expect_no_error(relab_sample_weighter(relab = dplyr::filter(relab_matrix, subject == "XBA"),
                                        group = "subject", w = w, K = 524))

  expect_no_error(relab_sample_weighter(relab = dplyr::filter(relab_matrix, subject == "XBA"),
                                        w = w, K = 524))
})

# Test arrange function
test_that("arrange_categories works", {
  expect_no_error(arrange_categories(relab_matrix = relab_matrix, arrange = TRUE, group = "subject", time = "timepoint"))
  expect_no_error(arrange_categories(relab_matrix = relab_matrix, arrange = TRUE, K = 524))

  expect_true(all(expect_no_error(arrange_categories(relab_matrix = relab_matrix, arrange = TRUE, K = 524)) ==
          expect_no_error(arrange_categories(relab_matrix = relab_matrix, arrange = "both", K = 524))))

  expect_no_error(arrange_categories(relab_matrix = relab_matrix, arrange = "both", group = "subject", time = "timepoint"))
  expect_no_error(arrange_categories(relab_matrix = relab_matrix, arrange = "both", K = 524))

  expect_no_error(arrange_categories(relab_matrix = relab_matrix, arrange = "horizontal", group = "subject", time = "timepoint"))
  expect_no_error(arrange_categories(relab_matrix = relab_matrix, arrange = "horizontal", K = 524))

  expect_no_error(arrange_categories(relab_matrix = relab_matrix, arrange = "vertical", group = "subject", time = "timepoint"))
  expect_no_error(arrange_categories(relab_matrix = relab_matrix, arrange = "vertical", K = 524))
})



