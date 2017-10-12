#' conivol: A package for the (bivariate) chi-bar-squared distribution and conic intrinsic volumes
#'
#' The conivol package provides functions for the chi-bar-squared distribution,
#' the bivariate chi-bar-squared distribution, and the conic intrinsic volumes.
#' Its main function is an estimator for the weights (conic intrinsic volumes)
#' of the bivariate chi-bar-squared distribution from sample data,
#' based on the expectation maximization (EM) method.
#'
#'
#' @section Chi-bar-squared distribution:
#' \itemize{
#'   \item \code{\link[conivol]{dchibarsq}}: evaluates the density
#'   \item \code{\link[conivol]{pchibarsq}}: evaluates the cumulative distribution function
#'   \item \code{\link[conivol]{rchibarsq}}: produces samples
#' }
#'
#' @section Bivariate chi-bar-squared distribution:
#' \itemize{
#'   \item \code{\link[conivol]{dbichibarsq}}: evaluates the density
#'   \item \code{\link[conivol]{pbichibarsq}}: evaluates the cumulative distribution function
#'   \item \code{\link[conivol]{rbichibarsq}}: produces samples
#' }
#'
#' @section Computing with conic intrinsic volumes:
#' \itemize{
#'   \item \code{\link[conivol]{comp_ivols_product}}: computes the intrinsic volumes of a product cone
#'                                    by convolving the intrinsic volumes of
#'                                    its elements
#'   \item \code{\link[conivol]{estimate_statdim_var}}: estimates the statistical dimension and
#'                                      the variance of the intrinsic volumes
#'                                      from samples of the corresponding
#'                                      bivariate chi-bar-squared distribution
#' }
#'
#' @section Computations involving products of circular cones:
#' \itemize{
#'   \item \code{\link[conivol]{circ_ivols}}: computes the intrinsic volumes of (a product of)
#'                           circular cones
#'   \item \code{\link[conivol]{circ_rbichibarsq}}: produces samples from the bivariate
#'                                  chi-bar-squared distribution with weights
#'                                  given by the conic intrinsic volumes of
#'                                  a product of circular cones
#' }
#'
#' @section Computations involving products of (polars of) Weyl chambers:
#' \itemize{
#'   \item \code{\link[conivol]{weyl_matrix}}: TBD
#'   \item \code{\link[conivol]{weyl_ivols}}: TBD
#' }
#'
#' @section General polyhedral cones:
#' \itemize{
#'   \item \code{\link[conivol]{polyh_reduce}}: computes a reduced representation
#'                                   of a polyhedral cone given by generators
#'   \item \code{\link[conivol]{polyh_samp_ivols_gen}}: produces samples from the
#'                                   intrinsic volumes distribution of
#'                                   a polyhedral cone given by generators
#'   \item \code{\link[conivol]{polyh_samp_ivols_ineq}}: produces samples from the
#'                                   intrinsic volumes distribution of
#'                                   a polyhedral cone given by inequalities
#'   \item \code{\link[conivol]{polyh_ivols_Bayes}}: TBD
#'   \item \code{\link[conivol]{polyh_rbichibarsq_gen}}: produces samples from the bivariate
#'                                   chi-bar-squared distribution with weights
#'                                   given by the conic intrinsic volumes of
#'                                   a polyhedral cone given by generators
#'   \item \code{\link[conivol]{polyh_rbichibarsq_ineq}}: produces samples from the bivariate
#'                                   chi-bar-squared distribution with weights
#'                                   given by the conic intrinsic volumes of
#'                                   a polyhedral cone given by inequalities
#' }
#'
#' @section Estimating the weights of the bivariate chi-bar-squared distribution:
#' \itemize{
#'   \item \code{\link[conivol]{prepare_data}}: evaluates the sample data of the bivariate chi-bar-squared
#'                              data (find the corresponding chi-squared density values);
#'                              this potentially time-consuming step is called during
#'                              \code{find_ivols_EM} and can be computed outside and passed
#'                              as parameter to avoid multiple calls, should \code{find_ivols_EM}
#'                              be called more than once
#'   \item \code{\link[conivol]{init_v}}: find an initial estimate of the weights, potentially
#'                        based on first and/or second moment
#'   \item \code{\link[conivol]{comp_loglike}}: compute the log-likelihood of a weight vector
#'                              for specific sample data
#'   \item \code{\link[conivol]{find_ivols_EM}}: produces EM-type iterates that may or may not converge
#'                               to the maximum likelihood estimate for the weights
#'                               of the bivariate chi-bar-squared distribution
#'                               from sample data; as the likelihood function is
#'                               quite flat around its maximum, the function supports
#'                               several ways to introduce some (well-founded) bias
#'                               and thus improve the estimate
#'   \item \code{\link[conivol]{find_ivols_jags}} TBD
#' }
#'
#' See the corresponding object documentation for more information. See also the following vignettes:
#' \itemize{
#'   \item \href{../doc/conic-intrinsic-volumes.html}{conic-intrinsic-volumes}:
#'         introduces conic intrinsic volumes and (bivariate) chi-bar-squared distributions
#'   \item \href{../doc/estim-conic-intrinsic-volumes-with-EM.html}{estim-conic-intrinsic-volumes-with-EM}:
#'         describes the details of the algorithm for finding the intrinsic volumes of closed
#'         convex cones from samples of the associated bivariate chi-bar-squared distribution
#'   \item \href{../doc/goodness-of-fit.html}{goodness-of-fit}:
#'         analyzes the goodness of fit of a vector of intrinsic volumes for given sample
#'         data from a bivariate chi-bar-squared distribution
#' }
#'
#' @docType package
#' @name conivol
NULL














