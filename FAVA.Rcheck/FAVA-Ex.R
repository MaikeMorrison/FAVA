pkgname <- "FAVA"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
base::assign(".ExTimings", "FAVA-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('FAVA')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("bootstrap_fava")
### * bootstrap_fava

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: bootstrap_fava
### Title: Statistically compare FAVA values between pairs of relative
###   abundance matrices.
### Aliases: bootstrap_fava

### ** Examples

# Statistically compare values of FAVA between
# subjects in the xue_microbiome_sample data:

 boot_out = bootstrap_fava(relab_matrix = xue_microbiome_sample,
               n_replicates = 100,
               seed = 1,
               group = "subject",
               K = 524,
               S = xue_species_similarity)

# Table of P-values comparing values of FAVA between group 1 and group 2:
 boot_out$P_values

 # Plots of the bootstrap distributions of differences in FAVA between each pair of matrices,
 # and how the true observed differences (red dots) compare to the distribution.
 boot_out$bootstrap_distribution_plot



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("bootstrap_fava", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fava")
### * fava

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fava
### Title: Compute the Fst of a matrix of compositional vectors
### Aliases: fava

### ** Examples

# Compute the Fst of
# the following compositional vectors:
q1 = c(1,   0,   0,   0)
q2 = c(0.5, 0.5, 0,   0)
q3 = c(1/4, 1/4, 1/4, 1/4)
q4 = c(0,   0,   1,   0)
relative_abundances = matrix(c(q1, q2, q3, q4),
                  byrow = TRUE, nrow = 4)

fava(relative_abundances)

# Incoporating weights:

# Compute fava ignoring
# rows 2 and 3
row_weights = c(0.5, 0, 0, 0.5)
fava(relative_abundances, w = row_weights)

# Compute fava assuming that
# categories 1 and 2 are identical:
similarity_matrix = diag(4)
similarity_matrix[1,2] = 1
similarity_matrix[2,1] = 1
fava(relative_abundances, S = similarity_matrix)

# Assume categories 1 and 2 are identical AND
# ignore rows 2 and 4:
row_weights = c(0.5, 0, 0.5, 0)
fava(relative_abundances, w = row_weights, S = similarity_matrix)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fava", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fava_norm")
### * fava_norm

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fava_norm
### Title: Compute the normalized Fst of a matrix of compositional vectors
### Aliases: fava_norm

### ** Examples

# Compute the weighted fava of
# the following compositional vectors:
q1 = c(1,   0,   0,   0)
q2 = c(0.5, 0.5, 0,   0)
q3 = c(1/4, 1/4, 1/4, 1/4)
q4 = c(0,   0,   1,   0)
relative_abundances = matrix(c(q1, q2, q3, q4),
                  byrow = TRUE, nrow = 4)

fava_norm(relative_abundances)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fava_norm", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("gini_simpson")
### * gini_simpson

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: gini_simpson
### Title: Compute the Gini-Simpson index of a compositional vector
### Aliases: gini_simpson

### ** Examples

# Compute unweighted Gini-Simpson index:
gini_simpson(q = c(0.4, 0.3, 0.3))

# Compute Gini-Simpson index assuming that
# categories 1 and 2 are identical:
similarity_matrix = diag(3)
similarity_matrix[1,2] = 1
similarity_matrix[2,1] = 1
gini_simpson(q = c(0.4, 0.3, 0.3), S = similarity_matrix)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("gini_simpson", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("gini_simpson_mean")
### * gini_simpson_mean

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: gini_simpson_mean
### Title: Compute the mean Gini-Simpson index of the rows in a matrix of
###   compositional vectors
### Aliases: gini_simpson_mean

### ** Examples

# To compute the mean Gini-Simpson index of
# the following compositional vectors...
q1 = c(1,   0,   0,   0)
q2 = c(0.5, 0.5, 0,   0)
q3 = c(1/4, 1/4, 1/4, 1/4)
q4 = c(0,   0,   1,   0)

# we could compute the mean manually:
mean(sapply(list(q1, q2, q3, q4), gini_simpson))

# Or we could use gini_simpson_mean:
relative_abundances = matrix(c(q1, q2, q3, q4),
                  byrow = TRUE, nrow = 4)

gini_simpson_mean(relative_abundances)

# Incoporating weights:

# Compute mean Gini-Simpson index ignoring
# rows 2 and 3
row_weights = c(0.5, 0, 0, 0.5)
gini_simpson_mean(relative_abundances, w = row_weights)

# Compute mean Gini-Simpson index assuming that
# categories 1 and 2 are identical:
similarity_matrix = diag(4)
similarity_matrix[1,2] = 1
similarity_matrix[2,1] = 1
gini_simpson_mean(relative_abundances, S = similarity_matrix)

# Assume categories 1 and 2 are identical AND
# ignore rows 2 and 4:
row_weights = c(0.5, 0, 0.5, 0)
gini_simpson_mean(relative_abundances, w = row_weights, S = similarity_matrix)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("gini_simpson_mean", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("gini_simpson_pooled")
### * gini_simpson_pooled

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: gini_simpson_pooled
### Title: Compute the pooled Gini-Simpson index of the rows in a matrix of
###   compositional vectors
### Aliases: gini_simpson_pooled

### ** Examples

# To compute the pooled Gini-Simpson index of
# the following compositional vectors...
q1 = c(1,   0,   0,   0)
q2 = c(0.5, 0.5, 0,   0)
q3 = c(1/4, 1/4, 1/4, 1/4)
q4 = c(0,   0,   1,   0)

# we could compute the mean manually:
qPooled = (q1 + q2 + q3 + q4)/4
gini_simpson(qPooled)

# Or we could use gini_simpson_pooled:
relative_abundances = matrix(c(q1, q2, q3, q4),
                  byrow = TRUE, nrow = 4)

gini_simpson_pooled(relative_abundances)

# Incoporating weights:

# Compute pooled Gini-Simpson index ignoring
# rows 2 and 3
row_weights = c(0.5, 0, 0, 0.5)
gini_simpson_pooled(relative_abundances, w = row_weights)

# Compute pooled Gini-Simpson index assuming that
# categories 1 and 2 are identical:
similarity_matrix = diag(4)
similarity_matrix[1,2] = 1
similarity_matrix[2,1] = 1
gini_simpson_pooled(relative_abundances, S = similarity_matrix)

# Assume categories 1 and 2 are identical AND
# ignore rows 2 and 4:
row_weights = c(0.5, 0, 0.5, 0)
gini_simpson_pooled(relative_abundances, w = row_weights, S = similarity_matrix)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("gini_simpson_pooled", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot_relabund")
### * plot_relabund

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot_relabund
### Title: Visualize a relative abundance matrix as a stacked bar plot.
### Aliases: plot_relabund

### ** Examples


# Make an example matrix of compositional data
# Each row is an individual. Rows sum to 1.
population_A = matrix(c(
    .5, .3, .2,
    .4, .2, .4,
    .5, .4, .1,
    .6, .1, .3,
    .2, 0, .8
  ),
  nrow = 5,
  byrow = TRUE
  )

  plot_relabund(relab_matrix = population_A,
              K = 3, # How many categories per vector?
              arrange = FALSE
              )
  plot_relabund(relab_matrix = population_A,
              K = 3, # How many categories per vector?
              arrange = "horizontal"
              )
  plot_relabund(relab_matrix = population_A,
              K = 3, # How many categories per vector?
              arrange = "vertical"
              )
   plot_relabund(relab_matrix = population_A,
              K = 3, # How many categories per vector?
              arrange = TRUE  # could also be "both"
              )


# You can modify the plot as you would any ggplot2 object
plot_relabund(relab_matrix = population_A,
              K = 3, # How many categories per vector?
              arrange = TRUE
              ) +
  # Below are example, optional modifications to the default plot
  ggplot2::ggtitle("Population A") +
  ggplot2::scale_fill_brewer("Blues") +
  ggplot2::scale_color_brewer("Blues") +
  ggplot2::xlab("Individuals")
  # Note that both scale_fill and scale_color are needed to change the color of the bars.


  # Plot a dataset which has 2 populations

  population_B = matrix(c(
    .9, 0, .1,
    .6, .4, 0,
    .7, 0, .3,
    .3, .4, .3,
    .5, .3, .2
  ),
  nrow = 5,
  byrow = TRUE
  )


  populations_AB = cbind(data.frame(c("A", "A", "A", "A", "A",
                                     "B", "B", "B", "B", "B")),
                         rbind(population_A, population_B))
  colnames(populations_AB) = c("population", "category_1", "category_2", "category_3")


 plot_relabund(relab_matrix = populations_AB, group = "population")
 plot_relabund(relab_matrix = populations_AB, group = "population", arrange = "vertical")
 plot_relabund(relab_matrix = populations_AB, group = "population", arrange = "horizontal")
 plot_relabund(relab_matrix = populations_AB, group = "population", arrange = "both")





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot_relabund", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("relab_phyloseq")
### * relab_phyloseq

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: relab_phyloseq
### Title: Generate a relative abundance matrix with sample meta data and
###   OTU abundances from a phyloseq object.
### Aliases: relab_phyloseq

### ** Examples

if (requireNamespace("phyloseq", quietly = TRUE)) {
  data(GlobalPatterns, package = "phyloseq")
  phyloseq_subset = phyloseq::subset_samples(GlobalPatterns,
                                             X.SampleID %in% c("CL3", "CC1"))
  otu_table = relab_phyloseq(phyloseq_subset)
  otu_table[, 1:10]
}



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("relab_phyloseq", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("time_weights")
### * time_weights

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: time_weights
### Title: Compute a normalized weighting vector based on a vector of
###   sampling times.
### Aliases: time_weights

### ** Examples

time_vector = c(1, 8, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
                32, 33, 34, 35, 36, 37, 38, 39, 44, 50, 57, 64)

time_weights(times = time_vector)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("time_weights", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("window_fava")
### * window_fava

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: window_fava
### Title: Compute FAVA in sliding windows.
### Aliases: window_fava

### ** Examples

A = matrix(c(.3,.7,0,.1,0,.9,.2,.5,.3,.1,.8,.1,.3,.4,.3,.6,.4,0,0,.5,.5),
           ncol = 3, byrow = TRUE)
window_out = window_fava(relab_matrix = A, window_size = 4, normalized = TRUE)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("window_fava", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("window_list")
### * window_list

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: window_list
### Title: Generate sliding windows of specified length given the maximum
###   number of samples
### Aliases: window_list

### ** Examples

window_list(window_size = 6, length = 40)
window_list(window_size = 6, length = 40, window_step = 2)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("window_list", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("window_plot")
### * window_plot

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: window_plot
### Title: Generate a plot of FAVA in sliding windows.
### Aliases: window_plot

### ** Examples

A = matrix(c(.3,.7,0,.1,0,.9,.2,.5,.3,.1,.8,.1,.3,.4,.3,.6,.4,0,0,.5,.5),
           ncol = 3, byrow = TRUE)
window_out = window_fava(relab_matrix = A, window_size = 4, normalized = TRUE)
window_out$window_data
window_out$window_plot



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("window_plot", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
