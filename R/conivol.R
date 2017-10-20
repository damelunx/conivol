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
#' @section Special classes of cones:
#' \itemize{
#'   \item \code{\link[conivol]{circ_ivols}}: computes the intrinsic volumes of (a product of)
#'                           circular cones
#'   \item \code{\link[conivol]{ellips_semiax}}: computes the semiaxes of an ellipsoidal cone
#'   \item \code{\link[conivol]{ellips_rbichibarsq}}: produces samples from the bivariate
#'                           chi-bar-squared distribution of an ellipsoidal cone
#'   \item \code{\link[conivol]{weyl_matrix}}: computes a matrix representation of
#'                           (a product of) Weyl chambers
#'   \item \code{\link[conivol]{weyl_ivols}}: computes the intrinsic volumes of (a product of)
#'                           Weyl chambers
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
#'   \item \code{\link[conivol]{polyh_ivols_bayes}}: generates functions for computing quantiles of marginals
#'                                   of the posterior distribution and for sampling
#'                                   from the posterior distribution,
#'                                   given samples of the intrinsic volumes distribution
#'   \item \code{\link[conivol]{polyh_ivols_stan}}: generates inputs for Stan
#'                                   (data list and model string or external file)
#'                                   for sampling from the posterior distribution,
#'                                   given samples of the intrinsic volumes distribution
#'                                   using a model that naturally implies log-concavity
#'                                   (and cannot be solved analytically)
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
#'   \item \code{\link[conivol]{prepare_data_em}}: evaluates the sample data of the bivariate chi-bar-squared
#'                              data (find the corresponding chi-squared density values);
#'                              this potentially time-consuming step is called during
#'                              \code{find_ivols_em} and can be computed outside and passed
#'                              as parameter to avoid multiple calls, should \code{find_ivols_em}
#'                              be called more than once
#'   \item \code{\link[conivol]{init_v}}: find an initial estimate of the weights, potentially
#'                              based on first and/or second moment
#'   \item \code{\link[conivol]{comp_loglike}}: compute the log-likelihood of a weight vector
#'                              for specific sample data
#'   \item \code{\link[conivol]{find_ivols_em}}: produces EM-type iterates that may or may not converge
#'                              to the maximum likelihood estimate for the weights
#'                              of the bivariate chi-bar-squared distribution
#'                              from sample data; as the likelihood function is
#'                              quite flat around its maximum, the function supports
#'                              several ways to introduce some bias
#'                              and thus improve the estimate
#'   \item \code{\link[conivol]{find_ivols_jags}}: generates inputs for JAGS
#'                              (data list and model string or external file)
#'                              for sampling from the posterior distribution,
#'                              given samples of the bivariate chi-bar-squared distribution
#'   \item \code{\link[conivol]{find_ivols_stan}}: generates inputs for Stan
#'                              (data list and model string or external file)
#'                              for sampling from the posterior distribution,
#'                              given samples of the bivariate chi-bar-squared distribution
#' }
#'
#'
#'
#' See the corresponding object documentation for more information. See also the following vignettes:
#' \itemize{
#'   \item \href{../doc/conic-intrinsic-volumes.html}{conic-intrinsic-volumes}:
#'         introduces conic intrinsic volumes and (bivariate) chi-bar-squared distributions
#'   \item \href{../doc/estim-conic-intrinsic-volumes-with-EM.html}{estim-conic-intrinsic-volumes-with-EM}:
#'         describes the details of the algorithm for finding the intrinsic volumes of closed
#'         convex cones from samples of the associated bivariate chi-bar-squared distribution
#'   \item \href{../doc/bayesian.html}{bayesian}:
#'         describes the Bayesian approach for evaluating sampling data, either from
#'         the intrinsic volumes distribution or from the chi-bar-squared distribution,
#'         and with or without enforcing log-concavity for the intrinsic volumes
#' }
#'
#' @docType package
#' @name conivol
NULL














