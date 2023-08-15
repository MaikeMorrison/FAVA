
<!-- README.md is generated from README.Rmd. Please edit that file -->

# FAVA

<a href='https://github.com/MaikeMorrison/FAVA'/><img src='man/figures/FAVA_logo_2' height="150" align="right" style="float:right; height:150px;" />

<!-- badges: start -->
<!-- badges: end -->

The *FAVA* R package implements the statistic FAVA, an
![F\_{ST}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;F_%7BST%7D "F_{ST}")-based
Assessment of Variability across vectors of relative Abundances, as well
as a suite of helper functions which enable the visualization and
statistical analysis of relative abundance data. The *FAVA* R package
package accompanies the paper, “FAVA: a tool to quantify compositional
variability in microbial communities” by Morrison et al. (2023?).

The *FAVA* R package includes the following core functions:

- `fava`: Quantify variability across many compositional vectors in a
  single, normalized index, called FAVA

- `bootstrap_fava`: Estimate the uncertainty in FAVA by generating
  bootstrap replicates of one or more relative abundance matrices and
  computing FAVA for each replicate matrix

- `window_fava`: Compute FAVA in sliding windows along the rows of a
  relative abundance matrix

- `plot_relabund`: Visualize a relative abundance matrix as a stacked
  bar plot

## Installation

You can install the development version of FAVA from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("MaikeMorrison/FAVA")
```

## Tutorial

A tutorial on the usage of *FAVA* with a focus on the analysis of
compositional data representing microbiome samples is available in the
`microbiome_tutorial` vignette, which is available at THIS LINK or via
the following code after package installation.

``` r
vignette("microbiome_tutorial", package = "FAVA")
```
