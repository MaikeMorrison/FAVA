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

# test_that("window_fst works - unnormalized, ungrouped", {
#   expect_no_error(window_fava(Q = Q, normalized = FALSE,
#                                        K = 524, window_size = 20, window_step = 10))
# })

# test_that("window_fst works - normalized, ungrouped", {
#   testthat::expect_no_error(window_fava(Q = Q, normalized = TRUE,
#                             K = 524, window_size = 20, window_step = 10))
# })
#
# test_that("window_fst works - normalized, grouped", {
#   testthat::expect_no_error(window_fava(Q = Q, normalized = TRUE,
#                                        K = 524, window_size = 6, window_step = 5,
#                                        group = "subject"))
# })
#
# test_that("window_fst works - unnormalized, grouped", {
#   testthat::expect_no_error(window_fava(Q = Q, normalized = FALSE,
#                                        K = 524, window_size = 6, window_step = 5,
#                                        group = "subject"))
# })
