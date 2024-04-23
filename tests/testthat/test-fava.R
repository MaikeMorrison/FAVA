# GINI SIMPSON ----------------------------------------------------------------

test_that("unweighted gini-simpson works", {
  expect_equal(gini_simpson(c(1/2, 1/2)), 0.5)
  expect_equal(gini_simpson(1), 0)
  expect_error(gini_simpson(NA))
  expect_error(gini_simpson(0))
  expect_error(gini_simpson(c(2,3,4)))
})


S = matrix(c(1, 1, 0,
             1, 1, 0,
             0, 0, 1),
           byrow = TRUE, nrow = 3)

test_that("weighted gini-simpson works", {
  expect_equal(gini_simpson(c(1/4, 1/4, 1/2), S = S), 0.5)
  expect_equal(gini_simpson(c(1/2, 1/2, 0), S = S), 0)
  expect_error(gini_simpson(c(NA, NA, NA), S = S))
  expect_error(gini_simpson(c(0,0,0), S = S))
  expect_error(gini_simpson(c(2,3,4), S = S))
})


# SIMPLE FAVA -----------------------------------------------------------------
A = diag(3)
B = matrix(rep(1/3, 9), ncol = 3)
C = matrix(c(0, 0, 1,
             0, 1, 0,
             1/2, 1/2, 0),
           ncol = 3, byrow = TRUE)
D = matrix(c(1, 0, 0,
             0, 1, 0,
             1/2, 1/2, 0),
           ncol = 3, byrow = TRUE)

w12 = c(1/2, 1/2, 0)
w13 = c(1/2, 0, 1/2)

test_that("fava works - unweighted", {
  expect_equal(fava(A), 1)
  expect_equal(fava(B), 0)
})

test_that("fava works - similarity", {
  expect_equal(fava(A, S = S), 1)
  expect_equal(fava(B, S = S), 0)
  expect_equal(fava(C, S = S), 1)
  expect_true(is.na(fava(D, S = S)))
})

test_that("fava works - w", {
  expect_equal(fava(A, w = w12), 1)
  expect_equal(fava(B,  w = w12), 0)
  expect_equal(fava(C, w = w12), 1)
})

test_that("fava works -  both", {
  expect_equal(fava(A, S = S, w = w13), 1)
  expect_equal(fava(B, S = S, w = w13), 0)
  expect_equal(fava(C, S = S, w = w13), 1)
  expect_true(is.na(fava(D, S = S, w = w13)))
})

test_that("fava works - normalized", {
  expect_equal(fava(A, normalized = TRUE), 1)
  expect_equal(fava(B, normalized = TRUE), 0)
  expect_equal(fava(C, normalized = TRUE), 1)
  expect_equal(fava(D, normalized = TRUE), fava_norm(D))
})


# GROUPED FAVA ----------------------------------------------------------------
library(dplyr)
gAB = rbind(A, B) %>%
  data.frame() %>%
  mutate(g = c(rep("A", 3), rep("B", 3)), .before = 1)

gABC = rbind(A, B, C) %>%
  data.frame() %>%
  mutate(g = c(rep("A", 3), rep("B", 3), rep("C", 3)), .before = 1)

gABCD = rbind(A, B, C, D) %>%
  data.frame() %>%
  mutate(g = c(rep("A", 3), rep("B", 3), rep("C", 3), rep("D", 3)), .before = 1)

test_that("grouped fava works - unweighted", {
  expect_equal(fava(gAB, group = "g"), data.frame(g = c("A", "B"),
                                                  FAVA = c(1, 0)))
  expect_no_error(fava(gABCD, group = "g"))
})

test_that("grouped fava works - similarity", {
  expect_equal(fava(gABCD, group = "g", S = S),
               data.frame(g = c("A", "B", "C", "D"),
                          FAVA = c(1, 0, 1, NaN)))
})

test_that("grouped fava works - w", {
  expect_equal(fava(gABC, group = "g", w = rep(w12, 3)),
               data.frame(g = c("A", "B", "C"),
                          FAVA = c(1, 0, 1)))
})

test_that("grouped fava works - both", {
  expect_equal(fava(gABCD, group = "g", S = S, w = rep(w13, 4)),
               data.frame(g = c("A", "B", "C", "D"),
                          FAVA = c(1, 0, 1, NaN)))
})

test_that("grouped fava works - normalized", {
  expect_equal(fava(gABCD, group = "g", normalized = TRUE),
               data.frame(g = c("A", "B", "C", "D"),
                          FAVA = c(1, 0, 1, fava_norm(D))))
})

test_that("grouped fava works - real data", {
  expect_equal(fava(xue_microbiome_sample, group = "subject", K = 524),
               data.frame(subject = c("XBA", "XDA", "XMA"),
                          FAVA = c(fava(filter(xue_microbiome_sample, subject == "XBA"), K = 524),
                                   fava(filter(xue_microbiome_sample, subject == "XDA"), K = 524),
                                   fava(filter(xue_microbiome_sample, subject == "XMA"), K = 524))))
})


# time series data works -----------------------------------------------------

tABC = gABC %>%
  mutate(.before = 1, timepoint = c(1,2,3,
                               1,2,3,
                               1,2,3))

wtime = time_weights(c(1,2,3))

test_that("time series fava works - ungrouped", {
  expect_equal(fava(data.frame(C) %>%
                      mutate(timepoint = c(1,2,3),
                             .before = 1),
                    time = "timepoint"),
               fava(C, w = wtime))

  # And in the wrong order:
  expect_equal(fava(data.frame(C) %>%
                      mutate(timepoint = c(3,2,1),
                             .before = 1),
                    time = "timepoint"),
               fava(C, w = wtime))

  expect_equal(fava(data.frame(C) %>%
                      mutate(timepoint = c(1,2,3),
                             .before = 1),
                    time = "timepoint", S = S),
               fava(C, w = wtime, S = S))

  expect_equal(fava(data.frame(D) %>%
                      mutate(timepoint = c(1,2,3),
                             .before = 1),
                    time = "timepoint"),
               fava(D, w = wtime))

  expect_equal(fava(data.frame(D) %>%
                      mutate(timepoint = c(1,2,3),
                             .before = 1),
                    time = "timepoint", S = S),
               fava(D, w = wtime, S = S))
})

test_that("time series fava works - grouped", {
  expect_equal(fava(tABC, group = "g", time = "timepoint"),
               data.frame(g = c("A", "B", "C"),
                          FAVA = c(1, 0, fava(C, w = wtime))))
})

test_that("time series fava works - grouped, S", {
  expect_equal(fava(tABC, group = "g", time = "timepoint", S = S),
               data.frame(g = c("A", "B", "C"),
                          FAVA = c(1, 0, 1)))
})

# MULTIPLE GROUPS ---------------------------------------------------------------
library(dplyr)
test_groups = xue_microbiome_sample %>%
  mutate(Abx = ifelse(timepoint < 29, "Before", ifelse(timepoint > 34, "After", "During")),
         .before = 1)

test_that("fava works with multiple groups", {
  expect_no_error(fava(test_groups, group = c("subject", "Abx"), K = 524))
  expect_equal(fava(test_groups, group = c("subject", "Abx"), K = 524)[[1,4]],
               fava(filter(test_groups, Abx == "Before", subject == "XBA"), K = 524))
  expect_equal(fava(test_groups, group = c("subject", "Abx"), K = 524)[[2,4]],
               fava(filter(test_groups, Abx == "During", subject == "XBA"), K = 524))
  expect_equal(fava(test_groups, group = c("subject", "Abx"), K = 524)[[3,4]],
               fava(filter(test_groups, Abx == "After", subject == "XBA"), K = 524))
})

test_groups_2 = test_groups
test_groups_2$Actinomyces_sp_58647 = test_groups_2$Actinomyces_sp_58647 + 0.2

test_that("fava works with multiple groups when renormalizing", {
  expect_warning(fava(test_groups_2, group = c("subject", "Abx"), K = 524))

  expect_warning(test_before_a <- fava(test_groups_2, group = c("subject", "Abx"), K = 524)[[1,4]])
  expect_warning(test_before_b <- fava(filter(test_groups_2, Abx == "Before", subject == "XBA"), K = 524))
  expect_equal(test_before_a, test_before_b)

  expect_warning(test_during_a <- fava(test_groups_2, group = c("subject", "Abx"), K = 524)[[2,4]])
  expect_warning(test_during_b <- fava(filter(test_groups_2, Abx == "During", subject == "XBA"), K = 524))
  expect_equal(test_during_a, test_during_b)

  expect_warning(test_after_a <- fava(test_groups_2, group = c("subject", "Abx"), K = 524)[[3,4]])
  expect_warning(test_after_b <- fava(filter(test_groups_2, Abx == "After", subject == "XBA"), K = 524))
  expect_equal(test_after_a, test_after_b)

})




test_that("fava works with multiple groups with weightings", {
  # JUST S
  expect_no_error(fava(test_groups, group = c("subject", "Abx"), K = 524, S = xue_species_similarity))
  expect_true(all(fava(test_groups, group = c("subject", "Abx"), K = 524, S = xue_species_similarity)$FAVA > 0))

  # JUST time
  expect_no_error(fava(test_groups, group = c("subject", "Abx"), K = 524, time = "timepoint"))
  expect_true(all(fava(test_groups %>% arrange(subject), group = c("subject", "Abx"), K = 524, time = "timepoint")$FAVA > 0))


  # BOTH S and time
  expect_no_error(fava(test_groups, group = c("subject", "Abx"), K = 524, time = "timepoint",
                       S = xue_species_similarity))
  expect_true(all(fava(test_groups %>% arrange(timepoint), group = c("subject", "Abx"), K = 524,
                       S = xue_species_similarity, time = "timepoint")$FAVA > 0))


  ### SWITCH THE ORDER OF THE GROUPING VARIABLES

  # JUST S
  expect_no_error(fava(test_groups, group = c("Abx", "subject"), K = 524, S = xue_species_similarity))
  expect_true(all(fava(test_groups, group = c("Abx", "subject"), K = 524, S = xue_species_similarity)$FAVA > 0))

  # JUST time
  expect_no_error(fava(test_groups, group = c("Abx", "subject"), K = 524, time = "timepoint"))
  expect_true(all(fava(test_groups %>% arrange(subject), group = c("Abx", "subject"), K = 524, time = "timepoint")$FAVA > 0))


  # BOTH S and time
  expect_no_error(fava(test_groups, group = c("Abx", "subject"), K = 524, time = "timepoint",
                       S = xue_species_similarity))
  expect_true(all(fava(test_groups %>% arrange(timepoint), group = c("Abx", "subject"), K = 524,
                       S = xue_species_similarity, time = "timepoint")$FAVA > 0))


})
