relab_matrix = xue_microbiome_sample %>% dplyr::filter(subject != "XMA") #
K = 524
group =  "subject" #c("subject", "Abx")#"subject"
S = xue_species_similarity
n_replicates = 3
normalized = FALSE
alternative = "two.sided"
# save_replicates = FALSE
w = NULL
time = "timepoint"

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


library(dplyr)
test_groups = xue_microbiome_sample %>%
  mutate(Abx = ifelse(timepoint < 29, "Before", ifelse(timepoint > 34, "After", "During")),
         .before = 1)


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
                                  K = 524, w = time_weights(times = test_groups$timepoint, group = paste0(test_groups$Abx, test_groups$subject)), seed = 1)$P_values)
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

