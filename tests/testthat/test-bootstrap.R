# relab_matrix = xue_microbiome_sample %>% dplyr::filter(subject != "XMA") #
# relab_matrix = test_groups
# K = 524
# group =  c("subject", "Abx")#"subject"
# S = xue_species_similarity
# n_replicates = 10
# normalized = FALSE
# alternative = "two.sided"
# # save_replicates = FALSE
# w = NULL
# time = "timepoint"


library(dplyr)


test_groups = xue_microbiome_sample %>%
  mutate(Abx = ifelse(timepoint < 29, "Before", ifelse(timepoint > 34, "After", "During")),
         .before = 1) %>% filter(Abx != "During")

test_that("Bootstrapping yields correct values for observed FAVA difference",{
  # group by subject (3 groups)
  favals = fava(relab_matrix = xue_microbiome_sample, S = xue_species_similarity,
                time = "timepoint", group = "subject")
  bs_out = bootstrap_fava(relab_matrix = xue_microbiome_sample, S = xue_species_similarity,
                          time = "timepoint", group = "subject",
                          n_replicates = 3)

  expect_true(all(round(bs_out$observed_difference$Difference, 5) ==
                    round(as.numeric(dist(favals$FAVA)) * c(-1,1,1), 5)))

  # group by Abx (2 groups)
  favals = fava(relab_matrix = test_groups,
                group = "Abx", K = 524)

  bs_out = bootstrap_fava(relab_matrix = test_groups,
                          group = "Abx", K = 524,
                          n_replicates = 3)

  expect_true(all(round(bs_out$observed_difference$Difference, 5) ==
                    round(dist(favals$FAVA), 5)))

})


xue_microbiome_sample_low_dim = test_groups[,1:53]
xue_microbiome_sample_low_dim[,4:53] =  test_groups[,4:53]/rowSums(test_groups[,4:53])

test_that("bootstrapping catches mispecified data",{
  # THREE GROUPS
  # S and data have different dimensions
  expect_warning(bootstrap_fava(relab_matrix = xue_microbiome_sample_low_dim, n_replicates = 3, group = "subject",
                                K = 50, S = xue_species_similarity))

  # time and normalization provided
  expect_error(bootstrap_fava(relab_matrix = test_groups, n_replicates = 3, group = "subject",
                              K = 50, S = xue_species_similarity, time = "timepoint", normalized = TRUE))

  # time and w provided
  expect_error(bootstrap_fava(relab_matrix = test_groups, n_replicates = 3, group = "subject",
                              K = 50, S = xue_species_similarity, time = "timepoint",
                              w = time_weights(xue_microbiome_sample$timepoint, group = xue_microbiome_sample$subject)))

  # normalization and w provided
  expect_error(bootstrap_fava(relab_matrix = test_groups, n_replicates = 3, group = "subject",
                              K = 50, S = xue_species_similarity, normalized = TRUE,
                              w = time_weights(xue_microbiome_sample$timepoint, group = xue_microbiome_sample$subject)))

  # TWO GROUPS
  # S and data have different dimensions
  expect_warning(bootstrap_fava(relab_matrix = xue_microbiome_sample_low_dim,
                                group = "Abx", n_replicates = 3,
                                K = 50, S = xue_species_similarity))

  # time and normalization provided
  expect_error(bootstrap_fava(relab_matrix = test_groups,
                              group = "Abx", n_replicates = 3,
                              K = 524, S = xue_species_similarity, time = "timepoint", normalized = TRUE))

  # time and w provided
  expect_error(bootstrap_fava(relab_matrix = test_groups,
                              group = "Abx", n_replicates = 3,
                              K = 524, S = xue_species_similarity, time = "timepoint",
                              w = time_weights(xue_microbiome_sample$timepoint, group = xue_microbiome_sample$subject)))

  # normalization and w provided
  expect_error(bootstrap_fava(relab_matrix = test_groups,
                              group = "Abx", n_replicates = 3,
                              K = 524, S = xue_species_similarity, normalized = TRUE,
                              w = time_weights(xue_microbiome_sample$timepoint, group = xue_microbiome_sample$subject)))

  # ONE GROUP
  expect_error(bootstrap_fava(relab_matrix = test_groups %>% dplyr::filter(Abx == "Before"),
                              group = "Abx", n_replicates = 3,
                              K = 524))

})

test_that("bootstrapping works for a matrix with one grouping var, multiple groups", {
  # NO WEIGHTS
  expect_no_error(bootstrap_fava(relab_matrix = xue_microbiome_sample, n_replicates = 3, group = "subject",
                                 K = 524))
  # NO K SPECIFIED WHEN EVERY COLUMN USED
  expect_identical(bootstrap_fava(relab_matrix = xue_microbiome_sample, n_replicates = 3, group = "subject", time = "timepoint", seed = 1)$P_values,
                   bootstrap_fava(relab_matrix = xue_microbiome_sample, n_replicates = 3, group = "subject", time = "timepoint", K = 524, seed = 1)$P_values)
  expect_identical(bootstrap_fava(relab_matrix = xue_microbiome_sample, n_replicates = 3, group = "subject", time = "timepoint",
                                  seed = 1, S = xue_species_similarity)$P_values,
                   bootstrap_fava(relab_matrix = xue_microbiome_sample, n_replicates = 3, group = "subject", time = "timepoint", K = 524,
                                  seed = 1, S = xue_species_similarity)$P_values)
  # SIMILARITY
  expect_no_error(bootstrap_fava(relab_matrix = xue_microbiome_sample, n_replicates = 3, group = "subject",
                                 K = 524, S = xue_species_similarity))
  # TIME
  expect_no_error(bootstrap_fava(relab_matrix = xue_microbiome_sample, n_replicates = 3, group = "subject",
                                 K = 524, time = "timepoint"))
  # W
  expect_no_error(bootstrap_fava(relab_matrix = xue_microbiome_sample, n_replicates = 3, group = "subject",
                                 K = 524, w = time_weights(times = xue_microbiome_sample$timepoint, group = xue_microbiome_sample$subject)))
  # CONFIRM TIME AND W IDENTICAL w same seed
  expect_identical(bootstrap_fava(relab_matrix = xue_microbiome_sample, n_replicates = 4, group = "subject",
                                  K = 524, time = "timepoint", seed = 1)$P_values,
                   bootstrap_fava(relab_matrix = xue_microbiome_sample, n_replicates = 4, group = "subject",
                                  K = 524, w = time_weights(times = xue_microbiome_sample$timepoint, group = xue_microbiome_sample$subject), seed = 1)$P_values)
  # SIMILARITY AND TIME
  expect_no_error(bootstrap_fava(relab_matrix = xue_microbiome_sample, n_replicates = 3, group = "subject",
                                 K = 524, time = "timepoint", S = xue_species_similarity))
  # SIMILARITY AND W
  expect_no_error(bootstrap_fava(relab_matrix = xue_microbiome_sample, n_replicates = 3, group = "subject",
                                 K = 524, w = time_weights(times = xue_microbiome_sample$timepoint, group = xue_microbiome_sample$subject),
                                 S = xue_species_similarity))
  # TIME AND W
  expect_error(bootstrap_fava(relab_matrix = xue_microbiome_sample, n_replicates = 3, group = "subject",
                              K = 524, w = time_weights(times = xue_microbiome_sample$timepoint, group = xue_microbiome_sample$subject),
                              time = "timepoint"))
  # NORMALIZED
  expect_no_error(bootstrap_fava(relab_matrix = xue_microbiome_sample, n_replicates = 3, group = "subject",
                                 K = 524,
                                 normalized = TRUE))
  # NORMALIZED AND SIMILARITY
  expect_error(bootstrap_fava(relab_matrix = xue_microbiome_sample, n_replicates = 3, group = "subject",
                              K = 524, S = xue_species_similarity,
                              normalized = TRUE))
  # NORMALIZED AND TIME
  expect_error(bootstrap_fava(relab_matrix = xue_microbiome_sample, n_replicates = 3, group = "subject",
                              K = 524,time = "timepoint",
                              normalized = TRUE))
  # NORMALIZED AND W
  expect_error(bootstrap_fava(relab_matrix = xue_microbiome_sample, n_replicates = 3, group = "subject",
                              K = 524, w = time_weights(times = xue_microbiome_sample$timepoint, group = xue_microbiome_sample$subject),
                              normalized = TRUE))

})

relab_2_groups = xue_microbiome_sample %>% dplyr::filter(subject %in% c("XDA", "XMA"))

test_that("bootstrapping works for a matrix with one grouping var, two groups", {
  # NO WEIGHTS
  expect_no_error(bootstrap_fava(relab_matrix = relab_2_groups, n_replicates = 3, group = "subject", K = 524))
  # NO K SPECIFIED WHEN EVERY COLUMN USED
  expect_identical(bootstrap_fava(relab_matrix = relab_2_groups, n_replicates = 3, group = "subject", time = "timepoint", seed = 1)$P_values,
                   bootstrap_fava(relab_matrix = relab_2_groups, n_replicates = 3, group = "subject", time = "timepoint", K = 524, seed = 1)$P_values)
  expect_identical(bootstrap_fava(relab_matrix = relab_2_groups, n_replicates = 3, group = "subject", time = "timepoint",
                                  seed = 1, S = xue_species_similarity)$P_values,
                   bootstrap_fava(relab_matrix = relab_2_groups, n_replicates = 3, group = "subject", time = "timepoint", K = 524,
                                  seed = 1, S = xue_species_similarity)$P_values)
  # SIMILARITY
  expect_no_error(bootstrap_fava(relab_matrix = relab_2_groups, n_replicates = 3, group = "subject",
                                 K = 524, S = xue_species_similarity))
  # TIME
  expect_no_error(bootstrap_fava(relab_matrix = relab_2_groups, n_replicates = 3, group = "subject",
                                 K = 524, time = "timepoint"))
  # W
  expect_no_error(bootstrap_fava(relab_matrix = relab_2_groups, n_replicates = 3, group = "subject",
                                 K = 524, w = time_weights(times = relab_2_groups$timepoint, group = relab_2_groups$subject)))
  # CONFIRM TIME AND W IDENTICAL w same seed
  expect_identical(bootstrap_fava(relab_matrix = relab_2_groups, n_replicates = 4, group = "subject",
                                  K = 524, time = "timepoint", seed = 1)$P_values,
                   bootstrap_fava(relab_matrix = relab_2_groups, n_replicates = 4, group = "subject",
                                  K = 524, w = time_weights(times = relab_2_groups$timepoint, group = relab_2_groups$subject), seed = 1)$P_values)
  # SIMILARITY AND TIME
  expect_no_error(bootstrap_fava(relab_matrix = relab_2_groups, n_replicates = 3, group = "subject",
                                 K = 524, time = "timepoint", S = xue_species_similarity))
  # SIMILARITY AND W
  expect_no_error(bootstrap_fava(relab_matrix = relab_2_groups, n_replicates = 3, group = "subject",
                                 K = 524, w = time_weights(times = relab_2_groups$timepoint, group = relab_2_groups$subject),
                                 S = xue_species_similarity))
  # TIME AND W
  expect_error(bootstrap_fava(relab_matrix = relab_2_groups, n_replicates = 3, group = "subject",
                              K = 524, w = time_weights(times = relab_2_groups$timepoint, group = relab_2_groups$subject),
                              time = "timepoint"))
  # NORMALIZED
  expect_no_error(bootstrap_fava(relab_matrix = relab_2_groups, n_replicates = 3, group = "subject",
                                 K = 524,
                                 normalized = TRUE))
  # NORMALIZED AND SIMILARITY
  expect_error(bootstrap_fava(relab_matrix = relab_2_groups, n_replicates = 3, group = "subject",
                              K = 524, S = xue_species_similarity,
                              normalized = TRUE))
  # NORMALIZED AND TIME
  expect_error(bootstrap_fava(relab_matrix = relab_2_groups, n_replicates = 3, group = "subject",
                              K = 524,time = "timepoint",
                              normalized = TRUE))
  # NORMALIZED AND W
  expect_error(bootstrap_fava(relab_matrix = relab_2_groups, n_replicates = 3, group = "subject",
                              K = 524, w = time_weights(times = relab_2_groups$timepoint, group = relab_2_groups$subject),
                              normalized = TRUE))

})



test_that("bootstrapping yields expected error for one matrix", {
  expect_error(bootstrap_fava(relab_matrix = dplyr::filter(xue_microbiome_sample, subject == "XBA"),
                              n_replicates = 3, K = 524, S = xue_species_similarity, group = "subject"))
})



test_that("bootstrapping works for a matrix with two grouping vars, multiple groups", {
  # NO WEIGHTS
  expect_no_error(bootstrap_fava(relab_matrix = test_groups, n_replicates = 3, group = c("subject", "Abx"),
                                 K = 524))
  # NO K SPECIFIED WHEN EVERY COLUMN USED
  expect_identical(bootstrap_fava(relab_matrix = test_groups, n_replicates = 3, group = c("subject", "Abx"), time = "timepoint", seed = 1)$P_values,
                   bootstrap_fava(relab_matrix = test_groups, n_replicates = 3, group = c("subject", "Abx"), time = "timepoint", K = 524, seed = 1)$P_values)
  expect_identical(bootstrap_fava(relab_matrix = test_groups, n_replicates = 3, group = c("subject", "Abx"), time = "timepoint",
                                  seed = 1, S = xue_species_similarity)$P_values,
                   bootstrap_fava(relab_matrix = test_groups, n_replicates = 3, group = c("subject", "Abx"), time = "timepoint", K = 524,
                                  seed = 1, S = xue_species_similarity)$P_values)
  # SIMILARITY
  expect_no_error(bootstrap_fava(relab_matrix = test_groups, n_replicates = 3, group = c("subject", "Abx"),
                                 K = 524, S = xue_species_similarity))
  # TIME
  expect_no_error(bootstrap_fava(relab_matrix = test_groups, n_replicates = 3, group = c("subject", "Abx"),
                                 K = 524, time = "timepoint"))
  # W
  expect_no_error(bootstrap_fava(relab_matrix = test_groups, n_replicates = 3, group = c("subject", "Abx"),
                                 K = 524, w = time_weights(times = test_groups$timepoint, group = test_groups$subject)))
  # CONFIRM TIME AND W IDENTICAL w same seed
  expect_identical(bootstrap_fava(relab_matrix = test_groups, n_replicates = 4, group = c("subject", "Abx"),
                                  K = 524, time = "timepoint", seed = 1)$P_values,
                   bootstrap_fava(relab_matrix = test_groups, n_replicates = 4, group = c("subject", "Abx"),
                                  K = 524, w = time_weights(times = test_groups$timepoint,
                                                            group = paste0(test_groups$Abx, test_groups$subject)),
                                  seed = 1)$P_values)
  # SIMILARITY AND TIME
  expect_no_error(bootstrap_fava(relab_matrix = test_groups, n_replicates = 3, group = c("subject", "Abx"),
                                 K = 524, time = "timepoint", S = xue_species_similarity))
  # SIMILARITY AND W
  expect_no_error(bootstrap_fava(relab_matrix = test_groups, n_replicates = 3, group = c("subject", "Abx"),
                                 K = 524, w = time_weights(times = test_groups$timepoint, group = test_groups$subject),
                                 S = xue_species_similarity))
  # TIME AND W
  expect_error(bootstrap_fava(relab_matrix = test_groups, n_replicates = 3, group = c("subject", "Abx"),
                              K = 524, w = time_weights(times = test_groups$timepoint, group = test_groups$subject),
                              time = "timepoint"))
  # NORMALIZED
  expect_no_error(bootstrap_fava(relab_matrix = test_groups, n_replicates = 3, group = c("subject", "Abx"),
                                 K = 524,
                                 normalized = TRUE))
  # NORMALIZED AND SIMILARITY
  expect_error(bootstrap_fava(relab_matrix = test_groups, n_replicates = 3, group = c("subject", "Abx"),
                              K = 524, S = xue_species_similarity,
                              normalized = TRUE))
  # NORMALIZED AND TIME
  expect_error(bootstrap_fava(relab_matrix = test_groups, n_replicates = 3, group = c("subject", "Abx"),
                              K = 524,time = "timepoint",
                              normalized = TRUE))
  # NORMALIZED AND W
  expect_error(bootstrap_fava(relab_matrix = test_groups, n_replicates = 3, group = c("subject", "Abx"),
                              K = 524, w = time_weights(times = test_groups$timepoint, group = test_groups$subject),
                              normalized = TRUE))

})

test_abx_factor = test_groups %>% dplyr::mutate(Abx = as.factor(Abx),
                                                subject = factor(subject, ordered = TRUE))

test_that("bootstrapping works if groups are factors", {
  # no weights
  expect_no_error(bootstrap_fava(relab_matrix = test_abx_factor, n_replicates = 3,
                                 group = "subject",
                                 K = 524))
  expect_no_error(bootstrap_fava(relab_matrix = test_abx_factor, n_replicates = 3,
                                 group = c("subject", "Abx"),
                                 K = 524))
  expect_no_error(bootstrap_fava(relab_matrix = test_abx_factor, n_replicates = 3,
                                 group = "Abx",
                                 K = 524))

  # with weights
  expect_no_error(bootstrap_fava(relab_matrix = test_abx_factor, n_replicates = 3,
                                 group = "subject",time = "timepoint",
                                 K = 524))
  # normalized
  expect_no_error(bootstrap_fava(relab_matrix = test_abx_factor, n_replicates = 3,
                                 group = "subject", normalized = TRUE,
                                 K = 524))
})

test_abx_numeric = test_abx_factor %>%
  dplyr::mutate(.before=Abx, test_group = as.numeric(subject),
                test_group_2 = as.numeric(Abx))

test_that("bootstrapping works if groups are numeric", {
  # no weights
  expect_no_error(bootstrap_fava(relab_matrix = test_abx_numeric, n_replicates = 3,
                                 group = "test_group",
                                 K = 524))
  expect_no_error(bootstrap_fava(relab_matrix = test_abx_numeric, n_replicates = 3,
                                 group = c("test_group_2", "test_group"),
                                 K = 524))
  expect_no_error(bootstrap_fava(relab_matrix = test_abx_numeric, n_replicates = 3,
                                 group = c("Abx", "test_group"),
                                 K = 524))

  expect_identical(expect_no_error(bootstrap_fava(relab_matrix = test_abx_numeric, n_replicates = 3,
                                                  group = c("test_group_2", "test_group"),
                                                  K = 524, seed = 1))$P_values[,2:3],
                   expect_no_error(bootstrap_fava(relab_matrix = test_abx_numeric, n_replicates = 3,
                                                  group = c("Abx", "subject"),
                                                  K = 524, seed = 1))$P_values[,2:3])
  # with weights
  expect_no_error(bootstrap_fava(relab_matrix = test_abx_numeric, n_replicates = 3,
                                 group = "test_group",time = "timepoint",
                                 K = 524))
  # normalized
  expect_no_error(bootstrap_fava(relab_matrix = test_abx_numeric, n_replicates = 3,
                                 group = "test_group", normalized = TRUE,
                                 K = 524))
})


test_abx_factor_space = test_abx_factor %>%
  mutate(subject_space = paste0(subject, " x"),
         subject_minus = paste0("-", subject),
         .before = timepoint)
test_that("bootstrapping tolerates group names with spaces", {
  # no weights
  expect_no_error(bootstrap_fava(relab_matrix = test_abx_factor_space, n_replicates = 3,
                                 group = "subject_space",
                                 K = 524))
  test_space = bootstrap_fava(relab_matrix = test_abx_factor_space, n_replicates = 3,
                              group = "subject_space",
                              K = 524)
  expect_true(all(test_space$observed_difference$Comparison %in% test_space$bootstrap_difference$Comparison))


  expect_no_error(bootstrap_fava(relab_matrix = test_abx_factor_space, n_replicates = 3,
                                 group = "subject_minus",
                                 K = 524))
  test_minus = bootstrap_fava(relab_matrix = test_abx_factor_space, n_replicates = 3,
                              group = "subject_minus",
                              K = 524)
  expect_true(all(test_minus$observed_difference$Comparison %in% test_minus$bootstrap_difference$Comparison))
})
