test_that("unweighted gini-simpson works", {
  expect_equal(gini_simpson(c(1/2, 1/2)), 0.5)
  expect_equal(gini_simpson(1), 0)
  expect_error(gini_simpson(NA))
  expect_error(gini_simpson(0))
  expect_error(gini_simpson(c(2,3,4)))
})


