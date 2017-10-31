
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- {r, echo = FALSE} -->
<!-- knitr::opts_chunk$set( -->
<!--   collapse = TRUE, -->
<!--   comment = "#>", -->
<!--   fig.path = "README_figures/" -->
<!-- ) -->
conivol
=======

``` r
getwd()
```

    ## [1] "/media/damlunx/Data/arbeit/R/conivol"

The conivol package provides functions for the chi-bar-squared distribution, the bivariate chi-bar-squared distribution, and the conic intrinsic volumes. Its main function is an estimator for the weights (conic intrinsic volumes) of the bivariate chi-bar-squared distribution from sample data, based on the expectation maximization (EM) method.

Installation
------------

You can install conivol from github with:

<!-- {r gh-installation, eval = FALSE} -->
<!-- # install.packages("devtools") -->
<!-- devtools::install_github("damelunx/conivol") -->
Functions
---------

The following functions are exported (sorted by context):

### Chi-bar-squared distribution

-   `dchibarsq`: evaluates the density
-   `pchibarsq`: evaluates the cumulative distribution function
-   `rchibarsq`: produces samples

<!-- See [this vignette](vignettes/conic-intrinsic-volumes.html) for more information about the chi-bar-squared distribution. -->
### Bivariate chi-bar-squared distribution

-   `dbichibarsq`: evaluates the density
-   `pbichibarsq`: evaluates the cumulative distribution function
-   `rbichibarsq`: produces samples

<!-- See [this vignette](vignettes/conic-intrinsic-volumes.html) for more information about the bivariate chi-bar-squared distribution. -->
### Computing with conic intrinsic volumes

-   `comp_ivols_product`: computes the intrinsic volumes of a product cone by convolving the intrinsic volumes of its elements
-   `estimate_statdim_var`: estimates the statistical dimension and the variance of the intrinsic volumes from samples of the corresponding bivariate chi-bar-squared distribution

### Computations involving products of circular cones

-   `circ_ivol`: computes the intrinsic volumes of (a product of) circular cones
-   `circ_rbichibarsq`: produces samples from the bivariate chi-bar-squared distribution with weights given by the conic intrinsic volumes of a product of circular cones
