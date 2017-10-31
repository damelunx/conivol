conivol: An R package for the (bivariate) chi-bar-squared distribution and conic intrinsic volumes
================

<!-- README.md is generated from README.Rmd. Please edit that file -->
This R package provides functions for the chi-bar-squared distribution, the bivariate chi-bar-squared distribution, and the conic intrinsic volumes. It supports standard functions for the density/cdf/sampling of the (bivariate) chi-bar-squared distribution, calculations and known formulas for special classes of intrinsic volumes of cones, sampling functions for ellipsoidal cones and general polyhedral cones, as well as functions for estimating intrinsic volumes either from direct samples of the intrinsic volumes distribution (in the case of polyhedral cones) or from samples of the corresponding bivariate chi-bar-squared distribution. The package supports point estimates as well as Bayesian estimates via JAGS and Stan.

Installation
------------

You can install `conivol` from github with:

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

# evaluate densities
dchibarsq(x, v)
#>  [1]        Inf 0.13419199 0.12061209 0.10708571 0.09155903 0.07548244
#>  [7] 0.06030318 0.04691767 0.03569914 0.02665556 0.01958501
dbichibarsq(xmat, v)
#>  [1]          Inf 0.0175809440 0.0182933272 0.0123633200 0.0070024293
#>  [6] 0.0036001390 0.0017409915 0.0008070907 0.0003627570 0.0001592392
#> [11] 0.0000686107

# evaluate cdfs
pchibarsq(x, v)
#>  [1] 0.1250000 0.3027630 0.4297465 0.5437900 0.6432281 0.7267271 0.7944988
#>  [8] 0.8479373 0.8890605 0.9200628 0.9430303
pbichibarsq(xmat, v)
#>  [1] 0.000000000 0.001292866 0.010039908 0.028749442 0.055055648
#>  [6] 0.085009193 0.115062662 0.142780036 0.166851524 0.186835988
#> [11] 0.202856633

# draw samples
rchibarsq(10,v)
#>  [1] 1.2116036 1.0261713 1.5308342 4.8864677 5.7462174 0.6808057 9.0402589
#>  [8] 0.0000000 4.3338161 1.4764667
rbichibarsq(10,v)
#>            [,1]        [,2]
#>  [1,] 7.3228926  0.44216238
#>  [2,] 3.9767080  0.00000000
#>  [3,] 2.5451880 13.62707115
#>  [4,] 0.4220977  5.62844702
#>  [5,] 0.2141612  5.22237939
#>  [6,] 7.5640044  0.07058617
#>  [7,] 1.3955951  1.77952743
#>  [8,] 1.8853436  8.22697856
#>  [9,] 9.9682674  1.72485695
#> [10,] 8.2812407  4.17028224
```

### Special classes of cones:

Some special classes of cones admit a direct computation of intrinsic volumes or a particularly simple sampling procedure for the bivariate chi-bar-squared distribution. These situations are covered by the following functions:

-   `prod_ivols`: computes the intrinsic volumes of a product cone by convolving the intrinsic volumes of its elements,
-   `circ_ivols`: computes the intrinsic volumes of (a product of) circular cones,
-   `ellips_semiax`, `ellips_rbichibarsq`: computes the semiaxes / produces samples from the bivariate chi-bar-squared distribution of an ellipsoidal cone,
-   `weyl_matrix`, `weyl_ivols`: computes a matrix representation / computes the intrinsic volumes of (a product of) Weyl chambers.

It is worth pointing out that with the class of direct product of Weyl chambers we obtain an interesting family of cones for which we have exact formulas; the same remark applies to the family of direct products of circular cones. We will use these functions to test the reconstruction algorithms. Ellipsoidal cones do not admit (simple) direct calculations of its intrinsic volumes.

**Usage:**

``` r
# the intrinsic volumes of the 4-dimensional circular cone of radius pi/4 a.k.a. Lorentz cone
circ_ivols(4,pi/4)
#> [1] 0.09084506 0.25000000 0.31830989 0.25000000 0.09084506

# the intrinsic volumes of the direct product of two 4-dimensional Lorentz cones
prod_ivols( list(circ_ivols(4,pi/4), circ_ivols(4,pi/4)) )
#> [1] 0.008252824 0.045422528 0.120333759 0.204577472 0.242826832 0.204577472
#> [7] 0.120333759 0.045422528 0.008252824

# computing the semiaxes of the ellipsoidal cone given by the linear image of the Lorentz cone
A <- matrix(sample(1:25),5,5)
ellips_semiax(A)
#> [1] 3.8936127 1.4359247 0.8710291 0.1561323

# draw samples of the bivariate chi-bar-squared distribution of the ellipsoidal cone
ellips_rbichibarsq(10, A)
#> $semiax
#> [1] 3.8936127 1.4359247 0.8710291 0.1561323
#> 
#> $samples
#>               [,1]       [,2]
#>  [1,] 1.151604e+00  0.6251532
#>  [2,] 9.327465e+00 -4.1944500
#>  [3,] 9.194609e-01  9.2061330
#>  [4,] 6.282846e-02  1.5018279
#>  [5,] 7.820756e-01  0.5276355
#>  [6,] 2.609564e+00  3.5429103
#>  [7,] 8.051955e-01  1.4340914
#>  [8,] 6.755320e-01 11.5557560
#>  [9,] 1.316501e+00 -0.1950536
#> [10,] 2.729160e-21  6.0779454

# compute the matrix of the product of some Weyl chambers
weyl_matrix( rep(3,4), c("BC","BCp","D","Dp"), product=TRUE)
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12]
#>  [1,]   -1    1    0    0    0    0    0    0    0   0.0   0.0     0
#>  [2,]    0   -1    1    0    0    0    0    0    0   0.0   0.0     0
#>  [3,]    0    0   -1    0    0    0    0    0    0   0.0   0.0     0
#>  [4,]    0    0    0    1    0    0    0    0    0   0.0   0.0     0
#>  [5,]    0    0    0    1    1    0    0    0    0   0.0   0.0     0
#>  [6,]    0    0    0    1    1    1    0    0    0   0.0   0.0     0
#>  [7,]    0    0    0    0    0    0   -1    1    0   0.0   0.0     0
#>  [8,]    0    0    0    0    0    0   -1   -1    1   0.0   0.0     0
#>  [9,]    0    0    0    0    0    0    0    0   -1   0.0   0.0     0
#> [10,]    0    0    0    0    0    0    0    0    0   0.5  -0.5     0
#> [11,]    0    0    0    0    0    0    0    0    0   0.5   0.5     0
#> [12,]    0    0    0    0    0    0    0    0    0   0.5   0.5     1

# compute the corresponding intrinsic volumes of the above cone
weyl_ivols( rep(3,4), c("BC","BCp","D","Dp"), product=TRUE)
#>  [1] 6.781684e-05 1.245569e-03 9.691780e-03 4.227024e-02 1.151364e-01
#>  [6] 2.064842e-01 2.502080e-01 2.064842e-01 1.151364e-01 4.227024e-02
#> [11] 9.691780e-03 1.245569e-03 6.781684e-05
```

### General polyhedral cones:

-   `polyh_reduce_gen`, `polyh_reduce_ineq`: compute a reduced representation of a polyhedral cone given by generators / inequalities
-   `polyh_rivols_gen`, `polyh_rivols_ineq`: produce samples from the intrinsic volumes distribution of a polyhedral cone given by generators / inequalities
-   `polyh_rbichibarsq_gen`, `polyh_rbichibarsq_ineq`: produce samples from the bivariate chi-bar-squared distribution with weights given by the conic intrinsic volumes of a polyhedral cone given by generators / inequalities
-   `polyh_bayes`: generates functions for computing quantiles of marginals of the posterior distribution and for sampling from the posterior distribution, given samples of the intrinsic volumes distribution (based on analytic solution)
-   `polyh_stan`: generates inputs for Stan (data list and model string or external file) for sampling from the posterior distribution, given samples of the intrinsic volumes distribution using a model that naturally implies log-concavity (and cannot be solved analytically)

**Usage:**

``` r
# example
```

### Estimating the weights of the bivariate chi-bar-squared distribution:

-   `estim_statdim_var`: estimates the statistical dimension and the variance of the intrinsic volumes from samples of the corresponding bivariate chi-bar-squared distribution
-   `init_ivols`: find an initial estimate of the weights, potentially based on first and/or second moment
-   `loglike_ivols`: compute the log-likelihood of a weight vector for specific sample data
-   `prepare_em`: evaluates the sample data of the bivariate chi-bar-squared data (find the corresponding chi-squared density values)
-   `estim_em`: produces EM-type iterates that may or may not converge to the maximum likelihood estimate for the weights of the bivariate chi-bar-squared distribution from sample data
-   `estim_jags`, `estim_stan`: generates inputs for JAGS / Stan (data list and model string or external file) for sampling from the posterior distribution of the intrinsic volumes, given samples of the bivariate chi-bar-squared distribution

**Usage:**

``` r
# example
```
