
<!-- README.md is generated from README.Rmd. Please edit that file -->
conivol
=======

The conivol package provides functions for the chi-bar-squared distribution, the bivariate chi-bar-squared distribution, and the conic intrinsic volumes. Its main function is an estimator for the weights (conic intrinsic volumes) of the bivariate chi-bar-squared distribution from sample data, based on the expectation maximization (EM) method.

Installation
------------

You can install conivol from github with:

``` r
# install.packages("devtools")
devtools::install_github("damelunx/conivol")
```

Functions
---------

The following functions are exported (sorted by context):

### Chi-bar-squared distribution

-   `dchibarsq`: evaluates the density
-   `pchibarsq`: evaluates the cumulative distribution function
-   `rchibarsq`: produces samples

See [this vignette](vignettes/conic-intrinsic-volumes.html) for more information about the chi-bar-squared distribution.

### Bivariate chi-bar-squared distribution

-   `dbichibarsq`: evaluates the density
-   `pbichibarsq`: evaluates the cumulative distribution function
-   `rbichibarsq`: produces samples

See [this vignette](vignettes/conic-intrinsic-volumes.html) for more information about the bivariate chi-bar-squared distribution.

### Computing with conic intrinsic volumes

-   `comp_ivols_product`: computes the intrinsic volumes of a product cone by convolving the intrinsic volumes of its elements
-   `estimate_statdim_var`: estimates the statistical dimension and the variance of the intrinsic volumes from samples of the corresponding bivariate chi-bar-squared distribution

### Computations involving products of circular cones

-   `circ_ivol`: computes the intrinsic volumes of (a product of) circular cones
-   `rbichibarsq_circ`: produces samples from the bivariate chi-bar-squared distribution with weights given by the conic intrinsic volumes of a product of circular cones

### General polyhedral cones

-   `rbichibarsq_polyh`: produces samples from the bivariate chi-bar-squared distribution with weights given by the conic intrinsic volumes of general polyhedral cones

### Estimating the weights of the bivariate chi-bar-squared distribution

-   `prepare_data`: evaluates the sample data of the bivariate chi-bar-squared data (find the corresponding chi-squared density values); this potentially time-consuming step is called during find\_ivols\_EM and can be computed outside and passed as parameter to avoid multiple calls should find\_ivols\_EM be called more than once
-   `init_v`: find an initial estimate of the weights, potentially based on first and/or second moment
-   `comp_loglike`: compute the log-likelihood of a weight vector for specific sample data
-   `find_ivols_EM`: produces EM-type iterates that may or may not converge to the maximum likelihood estimate for the weights of the bivariate chi-bar-squared distribution from sample data; as the likelihood function is quite flat around its maximum, the function supports several ways to introduce some (well-founded) bias and thus improve the estimate

Example
-------

TBD
