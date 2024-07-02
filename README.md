
<!-- README.md is generated from README.Rmd. Please edit that file -->

<a href='https://github.com/MaikeMorrison/FAVA'/><img src='man/figures/FAVA_logo_2.png' height="200" align="right" style="float:right; height:200px;" />

# FAVA

<!-- badges: start -->
<!-- badges: end -->

The *FAVA* R package implements the statistic FAVA, an $F_{ST}$-based
Assessment of Variability across vectors of relative Abundances, as well
as a suite of helper functions which enable the visualization and
statistical analysis of relative abundance data. The *FAVA* R package
accompanies the paper, “Quantifying compositional variability in
microbial communities with FAVA” by Morrison et al.

The *FAVA* R package includes the following core functions:

- `fava`: Quantify variability across many compositional vectors in a
  single, normalized index, called FAVA

- `bootstrap_fava`: Compare values of FAVA between pairs of abundance
  matrices

- `window_fava`: Compute FAVA in sliding windows along the rows of a
  relative abundance matrix

- `plot_relabund`: Visualize a relative abundance matrix as a stacked
  bar plot

## Installation

You can install FAVA from
[GitHub](https://github.com/MaikeMorrison/FAVA) with:

``` r

# First, install devtools if you haven't already:
# install.packages("devtools")

devtools::install_github("MaikeMorrison/FAVA")

# If you wish to access the tutorial (also accessible below) from within
# the package:

devtools::install_github("MaikeMorrison/FAVA", build_vignettes = TRUE)
```

## Tutorial

A tutorial on the usage of *FAVA* with a focus on the analysis of
compositional data representing microbiome samples is available in the
`microbiome_tutorial` vignette, which is available at [this
link](https://maikemorrison.com/files/microbiome_tutorial.html) or via
the following code after package installation.

``` r
vignette("microbiome_tutorial", package = "FAVA")
```
