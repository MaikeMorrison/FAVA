
<!-- README.md is generated from README.Rmd. Please edit that file -->

# FAVA

<!-- badges: start -->
<!-- badges: end -->

FAVA (Fst-based Assessment of Variability over relative Abundances) is a
suite of tools to quantify and compare the variability of compositional
datasets (i.e., matrices with rows that sum to 1). Included functions
simulate random compositional data, plot compositional data using
ggplot2, and calculate Fst, Fst weighted by distances between rows
and/or pairwise similarity among columns, and Fst/FstMax (a normalized
measure of variability). They also can generate bootstrap replicates of
one or more compositional datasets along with associated statistics.
This package accompanies the paper ‘’PAPER TITLE HERE’’ by Maike
Morrison, AND OTHERS. You can access the paper in NICE JOURNAL at this
link: <https://doi.org/>….

## Installation

You can install the development version of FAVA from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("MaikeMorrison/FAVA")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(FAVA)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
