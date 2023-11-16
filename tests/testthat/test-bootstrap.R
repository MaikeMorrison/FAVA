relab_matrix = xue_microbiome_sample #%>% filter(subject != "XMA")
K = 524
group = "subject"
S = xue_species_similarity
n_replicates = 100
normalized = FALSE
alternative = "two.sided"
save_replicates = FALSE

test_that("bootstrapping works for a matrix with groups", {
  expect_no_error(bootstrap_fava(relab_matrix = xue_microbiome_sample, n_replicates = 3, group = "subject",
                                 K = 524, S = xue_species_similarity))
})

test_that("bootstrapping works for a matrix with groups, unweighted", {
  expect_no_error(bootstrap_fava(relab_matrix = xue_microbiome_sample, n_replicates = 3, group = "subject", K = 524))
})

test_that("bootstrapping yields expected error for one matrix", {
  expect_error(bootstrap_fava(relab_matrix = dplyr::filter(xue_microbiome_sample, subject == "XBA"),
                                 n_replicates = 3, K = 524, S = xue_species_similarity, group = "subject"))
})


library(dplyr)
test_groups = xue_microbiome_sample %>%
  mutate(Abx = ifelse(timepoint < 29, "Before", ifelse(timepoint > 34, "After", "During")),
         .before = 1)

test_that("bootstrapping works for a matrix with multiple groups, unweighted", {
  expect_no_error(bootstrap_fava(relab_matrix = test_groups, n_replicates = 3,
                                 group = c("subject", "Abx"), K = 524))
})
