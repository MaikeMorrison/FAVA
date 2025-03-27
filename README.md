
<!-- README.md is generated from README.Rmd. Please edit that file -->

# FAVA

<a href='https://github.com/MaikeMorrison/FAVA'/><img src='man/figures/FAVA_logo_2.png' height="200" align="right" style="float:right; height:200px;" />

The *FAVA* R package implements the statistic FAVA, an FST-based
Assessment of Variability across vectors of relative Abundances, as well
as a suite of helper functions which enable the visualization and
statistical analysis of relative abundance data. The *FAVA* R package
accompanies the paper, “Quantifying compositional variability in
microbial communities with FAVA” by [Morrison et
al. (2025)](https://doi.org/10.1073/pnas.2413211122).

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

Install FAVA with:

``` r
install.packages("FAVA")
```

You can install the development version of FAVA from
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

The package website,
[maikemorrison.github.io/FAVA/](https://maikemorrison.github.io/FAVA/),
contains documentation and examples for all package functions. It also
contains a tutorial on the usage of *FAVA* for the analysis of
microbiome data. The tutorial vignette is available at [this
link](https://maikemorrison.github.io/FAVA/articles/microbiome_tutorial.html).
