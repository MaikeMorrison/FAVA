test_that("bootstrapping works for a matrix with multiple groups", {
  expect_no_error(bootstrap_fava(matrices = xue_microbiome_sample, n_replicates = 3, group = "subject",
                                 time = "timepoint", S = xue_species_similarity))
})

test_that("bootstrapping works for a matrix with multiple groups, unweighted", {
  expect_no_error(bootstrap_fava(matrices = xue_microbiome_sample, n_replicates = 3, group = "subject", K = 524))
})

test_that("bootstrapping works for one matrix", {
  expect_no_error(bootstrap_fava(matrices = dplyr::filter(xue_microbiome_sample, subject == "XBA"),
                                 n_replicates = 3, time = "timepoint", S = xue_species_similarity, K = 524))
})

test_that("bootstrapping works for one matrix, unweighted", {
  expect_no_error(bootstrap_fava(matrices = dplyr::filter(xue_microbiome_sample, subject == "XBA"),
                                 n_replicates = 3, K = 524))
})
