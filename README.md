conivol: An R package for the (bivariate) chi-bar-squared distribution and conic intrinsic volumes
================

<!-- README.md is generated from README.Rmd. Please edit that file -->
This R package provides functions for the chi-bar-squared distribution, the bivariate chi-bar-squared distribution, and the conic intrinsic volumes. It supports standard functions for the density/cdf/sampling of the (bivariate) chi-bar-squared distribution, calculations and known formulas for special classes of intrinsic volumes of cones, sampling functions for ellipsoidal cones and general polyhedral cones, as well as functions for estimating intrinsic volumes either from direct samples of the intrinsic volumes distribution (in the case of polyhedral cones) or from samples of the corresponding bivariate chi-bar-squared distribution. The package supports point estimates as well as Bayesian estimates via JAGS and Stan.

Installation
------------

You can install conivol from github with:

``` r
# install.packages("devtools")
devtools::install_github("damelunx/conivol")
```

Vignettes
---------

The following vignettes introduce the theory of intrinsic volumes, explain the idea behind the algorithm to reconstruct the intrinsic volumes from samples of the bivariate chi-bar-squared distribution, and the Bayesian approach to this reconstruction.

-   [Conic intrinsic volumes and (bivariate) chi-bar-squared distribution](../doc/conic-intrinsic-volumes.html): introduces conic intrinsic volumes and (bivariate) chi-bar-squared distributions, as well as the computations involving polyhedral cones
-   [Estimating conic intrinsic volumes via EM algorithm](../doc/estim-conic-intrinsic-volumes-with-EM.html): describes the details of the algorithm for finding the intrinsic volumes of closed convex cones from samples of the associated bivariate chi-bar-squared distribution
-   [Bayesian estimates for conic intrinsic volumes](../doc/bayesian.html): describes the Bayesian approach for reconstructing intrinsic volumes from sampling data, which can either be samples from the intrinsic volumes distribution (in the case of polyhedral cones), or from the bivariate chi-bar-squared distribution, and which can be with or without enforcing log-concavity of the intrinsic volumes

Functions
---------

In the following we list up the functions that are exported in the package (sorted by context), and include some examples to illustrate their use. See the above vignettes for more details about the underlying theory and algorithms.

### (Bivariate) Chi-bar-squared distribution:

`conivol` provides the standard support of the distributions (chi-bar-squared and bivariate chi-bar-squared) in the form of functions for the densities, cumulative distribution functions, and sampling methods:

-   `dchibarsq`, `dbichibarsq`: evaluate the densities of (bivariate) chi-bar-squared distributions,
-   `pchibarsq`, `pbichibarsq`: evaluate the dumulative distribution functions of (bivariate) chi-bar-squared distributions,
-   `rchibarsq`, `rbichibarsq`: produce samples of the (bivariate) chi-bar-squared distributions.

**Usage:**

``` r
# vector of weights
v <- rep(1,8)/8

# points to evaluate densities/cdfs (in the bavariate case two-column matrix)
x <- 0:10
xmat <- matrix(c(0:10,0:10),ncol=2)

# evaluate densities:
dchibarsq(x, v)
#>  [1]        Inf 0.13419199 0.12061209 0.10708571 0.09155903 0.07548244
#>  [7] 0.06030318 0.04691767 0.03569914 0.02665556 0.01958501
dbichibarsq(xmat, v)
#>  [1]          Inf 0.0175809440 0.0182933272 0.0123633200 0.0070024293
#>  [6] 0.0036001390 0.0017409915 0.0008070907 0.0003627570 0.0001592392
#> [11] 0.0000686107

# evaluate cdfs:
pchibarsq(x, v)
#>  [1] 0.1250000 0.3027630 0.4297465 0.5437900 0.6432281 0.7267271 0.7944988
#>  [8] 0.8479373 0.8890605 0.9200628 0.9430303
pbichibarsq(xmat, v)
#>  [1] 0.000000000 0.001292866 0.010039908 0.028749442 0.055055648
#>  [6] 0.085009193 0.115062662 0.142780036 0.166851524 0.186835988
#> [11] 0.202856633

# draw samples
rchibarsq(10,v)
#>  [1] 2.1121801 0.0000000 1.2120577 0.7420672 0.9077184 4.2519674 0.0000000
#>  [8] 1.9902648 1.1891448 1.3105069
rbichibarsq(10,v)
#>            [,1]     [,2]
#>  [1,] 1.2380520 7.450968
#>  [2,] 5.2617809 0.000000
#>  [3,] 1.6761319 7.201253
#>  [4,] 7.4492338 1.172290
#>  [5,] 6.0487199 2.855874
#>  [6,] 7.9452104 0.000000
#>  [7,] 2.5227745 3.167460
#>  [8,] 1.4987002 0.359360
#>  [9,] 0.2345424 2.180875
#> [10,] 2.8772112 6.265570
```

### Special classes of cones:

-   `prod_ivols`: computes the intrinsic volumes of a product cone by convolving the intrinsic volumes of its elements
-   `circ_ivols`: computes the intrinsic volumes of (a product of) circular cones
-   `ellips_semiax`, `ellips_rbichibarsq`: computes the semiaxes / produces samples from the bivariate chi-bar-squared distribution of an ellipsoidal cone
-   `weyl_matrix`, `weyl_ivols`: computes a matrix representation / computes the intrinsic volumes of (a product of) Weyl chambers

### General polyhedral cones:

-   `polyh_reduce_gen`, `polyh_reduce_ineq`: compute a reduced representation of a polyhedral cone given by generators / inequalities
-   `polyh_rivols_gen`, `polyh_rivols_ineq`: produce samples from the intrinsic volumes distribution of a polyhedral cone given by generators / inequalities
-   `polyh_rbichibarsq_gen`, `polyh_rbichibarsq_ineq`: produce samples from the bivariate chi-bar-squared distribution with weights given by the conic intrinsic volumes of a polyhedral cone given by generators / inequalities
-   `polyh_bayes`: generates functions for computing quantiles of marginals of the posterior distribution and for sampling from the posterior distribution, given samples of the intrinsic volumes distribution (based on analytic solution)
-   `polyh_stan`: generates inputs for Stan (data list and model string or external file) for sampling from the posterior distribution, given samples of the intrinsic volumes distribution using a model that naturally implies log-concavity (and cannot be solved analytically)

### Estimating the weights of the bivariate chi-bar-squared distribution:

-   `estim_statdim_var`: estimates the statistical dimension and the variance of the intrinsic volumes from samples of the corresponding bivariate chi-bar-squared distribution
-   `init_ivols`: find an initial estimate of the weights, potentially based on first and/or second moment
-   `loglike_ivols`: compute the log-likelihood of a weight vector for specific sample data
-   `prepare_em`: evaluates the sample data of the bivariate chi-bar-squared data (find the corresponding chi-squared density values)
-   `estim_em`: produces EM-type iterates that may or may not converge to the maximum likelihood estimate for the weights of the bivariate chi-bar-squared distribution from sample data
-   `estim_jags`, `estim_stan`: generates inputs for JAGS / Stan (data list and model string or external file) for sampling from the posterior distribution of the intrinsic volumes, given samples of the bivariate chi-bar-squared distribution

### Computing with conic intrinsic volumes

-   `comp_ivols_product`: computes the intrinsic volumes of a product cone by convolving the intrinsic volumes of its elements
-   `estimate_statdim_var`: estimates the statistical dimension and the variance of the intrinsic volumes from samples of the corresponding bivariate chi-bar-squared distribution

### Computations involving products of circular cones

-   `circ_ivol`: computes the intrinsic volumes of (a product of) circular cones
-   `circ_rbichibarsq`: produces samples from the bivariate chi-bar-squared distribution with weights given by the conic intrinsic volumes of a product of circular cones

### General polyhedral cones

-   `polyh_reduce`: computes a reduced representation of a polyhedral cone given by generators
-   `polyh_samp_ivol_gen`: produces samples from the intrinsic volumes distribution of a polyhedral cone given by generators
-   `polyh_samp_ivol_ineq`: produces samples from the intrinsic volumes distribution of a polyhedral cone given by inequalities
-   `polyh_rbichibarsq_gen`: produces samples from the bivariate chi-bar-squared distribution with weights given by the conic intrinsic volumes of a polyhedral cone given by generators
-   `polyh_rbichibarsq_ineq`: produces samples from the bivariate chi-bar-squared distribution with weights given by the conic intrinsic volumes of a polyhedral cone given by inequalities

### Estimating the weights of the bivariate chi-bar-squared distribution

-   `prepare_data`: evaluates the sample data of the bivariate chi-bar-squared data (find the corresponding chi-squared density values); this potentially time-consuming step is called during find\_ivols\_EM and can be computed outside and passed as parameter to avoid multiple calls should find\_ivols\_EM be called more than once
-   `init_v`: find an initial estimate of the weights, potentially based on first and/or second moment
-   `comp_loglike`: compute the log-likelihood of a weight vector for specific sample data
-   `find_ivols_EM`: produces EM-type iterates that may or may not converge to the maximum likelihood estimate for the weights of the bivariate chi-bar-squared distribution from sample data; as the likelihood function is quite flat around its maximum, the function supports several ways to introduce some (well-founded) bias and thus improve the estimate. See the example below and [this vignette](vignettes/estim-conic-intrinsic-volumes-with-EM.html) for a description of this parameter tuning.

Example
-------

<!-- ```{r echo=FALSE, message=FALSE} -->
<!-- library(conivol) -->
<!-- library(Rmosek) -->
<!-- library(tidyverse) -->
<!-- ``` -->
<!-- We consider the product of two circular cones, one 5-dimensional and another -->
<!-- 8-dimensional. -->
<!-- ```{r} -->
<!-- # specify the dimensions -->
<!-- D <- c(5,8) -->
<!-- # specify the angles -->
<!-- alpha <- c( 0.7*pi/2, 0.8*pi/2 ) -->
<!-- # get the exact intrinsic volumes -->
<!-- v <- circ_ivol(D, alpha, product = TRUE) -->
<!-- # plot the values and its logarithms -->
<!-- d <- sum(D) -->
<!-- ggplot(tibble(k=0:d, v=v), aes(x=k,y=v)) + -->
<!--     geom_line() + -->
<!--     theme_bw() -->
<!-- ggplot(tibble(k=0:d, `log(v)`=log(v)), aes(x=k,y=`log(v)`)) + -->
<!--     geom_line() + -->
<!--     theme_bw() -->
<!-- ``` -->
<!-- The goal is to reconstruct these values from samples from the bivariate -->
<!-- chi-bar-squared distribution with the above intrinsic volumes as weights. -->
<!-- ```{r} -->
<!-- # set sample size -->
<!-- n <- 10^5 -->
<!-- # obtain sample of the specified size -->
<!-- set.seed(1234) -->
<!-- m_samp <- circ_rbichibarsq(n,D,alpha) -->
<!-- # scatter plot of the sample -->
<!-- ggplot(as_tibble(m_samp), aes(V1,V2)) + geom_point(alpha=.02) + -->
<!--     theme_bw() + -->
<!--     theme(axis.title.x=element_blank(),axis.title.y=element_blank()) -->
<!-- ``` -->
