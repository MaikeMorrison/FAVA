library(testthat)
library(FAVA)
Q = xue_microbiome_sample



test_that("window_list works", {
  expect_no_error(window_list(window_size = 5, length = 20, window_step = 4))
  expect_error(window_list(window_size = NA, length = 20, window_step = 4))
  expect_error(window_list(window_size = 5, length = NA, window_step = 4))
  expect_error(window_list(window_size = 5, length = 20, window_step = NA))

  expect_equal(sum(unlist(window_list(window_size = 5, length = 20, window_step = 4))), 180)
  expect_equal(sum(unlist((window_list(window_size = 5, length = 20)))), 840)
})

test_that("window_fava works - unnormalized, ungrouped", {
  expect_no_error(window_fava(relab_matrix = Q, normalized = FALSE,
                              K = 524, window_size = 20, window_step = 10))
  expect_true(is.list(window_fava(relab_matrix = Q, normalized = FALSE,
                                  K = 524, window_size = 20, window_step = 10)))
  expect_equal(window_fava(relab_matrix = Q, normalized = FALSE,
                           K = 524, window_size = 20, window_step = 10)$window_data[1,1],
               fava(Q[1:20,], K = 524))
})

test_that("window_fava works - normalized, ungrouped", {
  expect_no_error(window_fava(relab_matrix = Q, normalized = TRUE,
                              K = 524, window_size = 20, window_step = 10))
  expect_true(is.list(window_fava(relab_matrix = Q, normalized = TRUE,
                                  K = 524, window_size = 20, window_step = 10)))
  expect_equal(window_fava(relab_matrix = Q, normalized = TRUE,
                           K = 524, window_size = 20, window_step = 10)$window_data[1,1],
               fava(Q[1:20,], K = 524, normalized = TRUE))
})

test_that("window_fava works - unnormalized, grouped", {
  expect_no_error(window_fava(relab_matrix = Q, normalized = FALSE, group = "subject",
                              K = 524, window_size = 8, window_step = 2))
  expect_true(is.list(window_fava(relab_matrix = Q, normalized = FALSE, group = "subject",
                                  K = 524, window_size = 8, window_step = 2)))
  expect_equal(window_fava(relab_matrix = Q, normalized = FALSE, group = "subject",
                           K = 524, window_size = 8, window_step = 2)$window_data[1,2],
               fava(Q[1:8,], K = 524))
})

test_that("window_fava works - normalized, grouped", {
  expect_no_error(window_fava(relab_matrix = Q, normalized = TRUE, group = "subject",
                              K = 524, window_size = 8, window_step = 2))
  expect_true(is.list(window_fava(relab_matrix = Q, normalized = TRUE, group = "subject",
                                  K = 524, window_size = 8, window_step = 2)))
  expect_equal(window_fava(relab_matrix = Q, normalized = TRUE, group = "subject",
                           K = 524, window_size = 8, window_step = 2)$window_data[1,2],
               fava(Q[1:8,], K = 524, normalized = TRUE))
})

test_that("window_fava works - unnormalized, grouped, time-weighted", {
  expect_no_error(window_fava(relab_matrix = Q, group = "subject", time = "timepoint",
                              K = 524, window_size = 8, window_step = 2))
  expect_true(is.list(window_fava(relab_matrix = Q, group = "subject", time = "timepoint",
                                  K = 524, window_size = 8, window_step = 2)))
  expect_equal(as.numeric(window_fava(relab_matrix = Q, group = "subject", time = "timepoint",
                                      K = 524, window_size = 8, window_step = 2)$window_data[1,2]),
               fava(Q[1:8,], K = 524, time = "timepoint"))
})

Q_neg = Q
Q_neg$timepoint = Q_neg$timepoint - 29
test_that("window_fava works - unnormalized, grouped, time-weighted, with NEGATIVE TIMES", {
  expect_no_error(window_fava(relab_matrix = Q_neg, group = "subject", time = "timepoint",
                              K = 524, window_size = 8, window_step = 2))
  expect_true(is.list(window_fava(relab_matrix = Q_neg, group = "subject", time = "timepoint",
                                  K = 524, window_size = 8, window_step = 2)))
  expect_equal(as.numeric(window_fava(relab_matrix = Q_neg, group = "subject", time = "timepoint",
                                      K = 524, window_size = 8, window_step = 2)$window_data[1,2]),
               fava(Q_neg[1:8,], K = 524, time = "timepoint"))
})

test_that("window_fava works - unnormalized, grouped, time- and similarity-weighted, with NEGATIVE TIMES", {
  expect_no_error(window_fava(relab_matrix = Q_neg, group = "subject", time = "timepoint", S = xue_species_similarity,
                              K = 524, window_size = 6, window_step = 2))
  expect_true(is.list(window_fava(relab_matrix = Q_neg, group = "subject", time = "timepoint", S = xue_species_similarity,
                                  K = 524, window_size = 8, window_step = 2)))
  expect_equal(as.numeric(window_fava(relab_matrix = Q_neg, group = "subject", time = "timepoint", S = xue_species_similarity,
                                      K = 524, window_size = 8, window_step = 2)$window_data[1,2]),
               fava(Q_neg[1:8,], K = 524, S = xue_species_similarity, time = "timepoint"))
})


set.seed(1)
ransamp = sample(1:77,77)
test_that("window_fava works - unnormalized, grouped, time-weighted and out of order", {
  expect_no_error(window_fava(relab_matrix = Q[ransamp,], group = "subject", time = "timepoint",
                              K = 524, window_size = 8, window_step = 2))
  expect_true(is.list(window_fava(relab_matrix = Q[ransamp,], group = "subject", time = "timepoint",
                                  K = 524, window_size = 8, window_step = 2)))
  expect_equal(as.numeric(window_fava(relab_matrix = Q[ransamp,], group = "subject", time = "timepoint",
                                      K = 524, window_size = 8, window_step = 2)$window_data[1,2]),
               fava(Q[1:8,], K = 524, time = "timepoint"))
})


set.seed(1)
Q_mult = dplyr::mutate(Q, random_variable = c(rep(c(1,2), 38), 1), .before = timepoint)
test_that("window_fava works - MULTIPLE groups", {
  expect_no_error(window_fava(relab_matrix = Q_mult, group = c("subject", "random_variable"),
                              K = 524, window_size = 4, window_step = 2))
  expect_true(is.list(window_fava(relab_matrix = Q_mult, group = c("subject", "random_variable"),
                                  K = 524, window_size = 4, window_step = 2)))
})

test_that("window_fava works - MULTIPLE groups, weighted", {
  expect_no_error(window_fava(relab_matrix = Q_mult, group = c("subject", "random_variable"),
                              K = 524, window_size = 4, window_step = 2, time = "timepoint"))

  expect_no_error(window_fava(relab_matrix = Q_mult, group = c("subject", "random_variable"),
                              K = 524, window_size = 4, window_step = 2, S = xue_species_similarity))

  expect_no_error(window_fava(relab_matrix = Q_mult, group = c("subject", "random_variable"),
                              K = 524, window_size = 4, window_step = 2, time = "timepoint", S = xue_species_similarity))
})

test_that("window_fava works - MULTIPLE groups, weighted, and out of order", {
  expect_no_error(window_fava(relab_matrix = Q_mult[ransamp,], group = c("subject", "random_variable"),
                              K = 524, window_size = 4, window_step = 2, time = "timepoint"))

  expect_no_error(window_fava(relab_matrix = Q_mult[ransamp,], group = c("subject", "random_variable"),
                              K = 524, window_size = 4, window_step = 2, S = xue_species_similarity))

  expect_no_error(window_fava(relab_matrix = Q_mult[ransamp,], group = c("subject", "random_variable"),
                              K = 524, window_size = 4, window_step = 2, time = "timepoint", S = xue_species_similarity))
})


w = time_weights(times = Q$timepoint, group = Q$subject)
test_that("window_fava works - unnormalized, grouped, generic weighting vector", {
  expect_no_error(window_fava(relab_matrix = Q, group = "subject", w = w,
                              K = 524, window_size = 8, window_step = 2))
  expect_true(is.list(window_fava(relab_matrix = Q, group = "subject", w = w,
                                  K = 524, window_size = 8, window_step = 2)))
  expect_equal(window_fava(relab_matrix = Q, group = "subject", w = w,
                           K = 524, window_size = 8, window_step = 2)$window_data[1,2],
               fava(Q[1:8,], K = 524, w = w[1:8]/sum(w[1:8])))

})



# NOTE:
# We do not expect window_fava to give the same results when time is specified as when
# the analogous weighting vector is specified because, when time is specified, the
# first and last timepoints of each window are down weighted because they do not inform
# the composition before or after the window. However, when we compute the weighting vector
# of the entire interval at once, these samples are often in the middle of the interval,
# and are therefore more informative and have a greater weighting.

# XBA = dplyr::filter(Q, subject == "XBA")
#
# w = time_weights(times = XBA$timepoint)
#
# window_indices = window_list(window_size = 8, length = nrow(XBA), window_step = 2)
#
# window1 = window_indices[[1]]
#
# time_weights_1 = time_weights(unlist(XBA[window1, "timepoint"]))
#
# w[window1]/sum(w[window1])
#
# test_that("window_fava works the same when provided with time vs. equivalent weights"){
# }

window_out_grouped = window_fava(relab_matrix = Q, group = "subject", time = "timepoint",
                                 K = 524, window_size = 8, window_step = 2)

window_out_ungrouped = window_fava(relab_matrix = dplyr::filter(Q, subject == "XBA"),
                                   time = "timepoint",
                                   K = 524, window_size = 8, window_step = 2)

test_that("window_plot works - grouped", {
  expect_no_error(window_plot(window_out_grouped$window_data))
  expect_no_error(window_plot(window_out_ungrouped$window_data))
})
