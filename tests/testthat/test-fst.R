test_that("unweighted heterozygosity works", {
  expect_equal(het(c(1/2, 1/2)), 0.5)
  expect_equal(het(1), 0)
  expect_error(het(NA))
  expect_error(het(0))
  expect_error(het(c(2,3,4)))
})
