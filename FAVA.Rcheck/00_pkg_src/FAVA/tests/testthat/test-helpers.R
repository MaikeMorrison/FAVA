
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
  expect_no_condition(relab_sample_weighter(relab = relab_matrix, group = "subject", time = "timepoint"))
  # K and time, no group
  expect_error(relab_sample_weighter(relab = relab_matrix, K = 524, time = "timepoint"))
  # K and group and time
  expect_no_error(relab_sample_weighter(relab = relab_matrix, K = 524, group = "subject", time = "timepoint"))

  # Specify w instead of time ---------------------------------
  expect_no_error(relab_sample_weighter(relab = relab_matrix, group = "subject", w = w, K = 524))

  # One group ------------------------------------------------------------------
  # K and group and time
  expect_no_error(relab_sample_weighter(relab = dplyr::filter(relab_matrix, subject == "XBA"),
                                                              K = 524, group = "subject",
                                            time = "timepoint"))
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
  expect_no_condition(arrange_categories(relab_matrix = relab_matrix, arrange = TRUE, group = "subject", time = "timepoint"))
  expect_no_condition(arrange_categories(relab_matrix = relab_matrix, arrange = TRUE, K = 524))

  expect_true(all(expect_no_condition(arrange_categories(relab_matrix = relab_matrix, arrange = TRUE, K = 524)) ==
          expect_no_condition(arrange_categories(relab_matrix = relab_matrix, arrange = "both", K = 524))))

  expect_no_condition(arrange_categories(relab_matrix = relab_matrix, arrange = "both", group = "subject", time = "timepoint"))
  expect_no_condition(arrange_categories(relab_matrix = relab_matrix, arrange = "both", K = 524))

  expect_no_condition(arrange_categories(relab_matrix = relab_matrix, arrange = "horizontal", group = "subject", time = "timepoint"))
  expect_no_condition(arrange_categories(relab_matrix = relab_matrix, arrange = "horizontal", K = 524))

  expect_no_condition(arrange_categories(relab_matrix = relab_matrix, arrange = "vertical", group = "subject", time = "timepoint"))
  expect_no_condition(arrange_categories(relab_matrix = relab_matrix, arrange = "vertical", K = 524))
})



# Test the similarity matrix checker
test_that("S_checker works", {

  # Correct matrix goes through
  expect_equal(S_checker(S = xue_species_similarity, K = 524, relab_matrix = xue_microbiome_sample), as.matrix(xue_species_similarity))
  expect_equal(S_checker(S = xue_species_similarity, K = 524), as.matrix(xue_species_similarity))

  # Shuffled matrix is resorted
  scramble = sample(1:524)
  expect_warning(S_checker(S = xue_species_similarity[scramble, scramble], K = 524, relab_matrix = xue_microbiome_sample))
  expect_warning(expect_equal(S_checker(S = xue_species_similarity[scramble, scramble], K = 524, relab_matrix = xue_microbiome_sample), as.matrix(xue_species_similarity)))

  # asymmetry is warned of
  asymmetric = xue_species_similarity
  asymmetric[1,2] = 1
  expect_warning(S_checker(S = asymmetric, K = 524, relab_matrix = xue_microbiome_sample))

  # non-1 diagonal elements break
  non1 = xue_species_similarity
  non1[20, 20] = 0.5
  expect_error(S_checker(S = non1, K = 524, relab_matrix = xue_microbiome_sample))

  # any negative elements break
  neg = xue_species_similarity
  neg[20, 24] = -0.2
  neg[24, 20] = -0.2
  expect_error(S_checker(S = neg, K = 524, relab_matrix = xue_microbiome_sample))

  # any greater than 1 elements break
  big = xue_species_similarity
  big[20, 24] = 2
  big[24, 20] = 2
  expect_error(S_checker(S = big, K = 524, relab_matrix = xue_microbiome_sample))


})


