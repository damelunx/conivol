#' conivol: A package for the (bivariate) chi-bar-squared distribution and conic intrinsic volumes
#'
#' The conivol package provides functions for the chi-bar-squared distribution,
#' the bivariate chi-bar-squared distribution, and the conic intrinsic volumes.
#'
#' \code{conivol} supports standard functions for the density/cdf/sampling of the (bivariate)
#' chi-bar-squared distribution, calculations and known formulas for special classes
#' of intrinsic volumes of cones, sampling functions for ellipsoidal cones and
#' general polyhedral cones, as well as functions for estimating intrinsic volumes
#' either from direct samples of the intrinsic volumes distribution
#' (in the case of polyhedral cones) or from samples of the corresponding
#' bivariate chi-bar-squared distribution. The package supports point estimates
#' as well as Bayesian estimates via JAGS and Stan.
#'
#' @section (Bivariate) Chi-bar-squared distribution:
#' \itemize{
#'   \item \code{\link[conivol]{dchibarsq}} / \code{\link[conivol]{pchibarsq}} /
#'         \code{\link[conivol]{rchibarsq}}: evaluates the density /
#'         evaluates the cumulative distribution function / produces samples
#'         of the chi-bar-squared distribution
#'   \item \code{\link[conivol]{dbichibarsq}} / \code{\link[conivol]{pbichibarsq}} /
#'         \code{\link[conivol]{rbichibarsq}}: evaluates the density /
#'         evaluates the cumulative distribution function / produces samples
#'         of the bivariate chi-bar-squared distribution
#' }
#'
#' @section Special classes of cones:
#' \itemize{
#'   \item \code{\link[conivol]{prod_ivols}}: computes the intrinsic volumes of a product cone
#'                                    by convolving the intrinsic volumes of
#'                                    its elements
#'   \item \code{\link[conivol]{circ_ivols}}: computes the intrinsic volumes of (a product of)
#'                           circular cones
#'   \item \code{\link[conivol]{ellips_semiax}} / \code{\link[conivol]{ellips_rbichibarsq}}:
#'                           computes the semiaxes / produces samples from the bivariate
#'                           chi-bar-squared distribution of an ellipsoidal cone
#'   \item \code{\link[conivol]{weyl_matrix}} / \code{\link[conivol]{weyl_ivols}}:
#'                           computes a matrix representation / computes the
#'                           intrinsic volumes of (a product of) Weyl chambers
#' }
#'
#' @section General polyhedral cones:
#' \itemize{
#'   \item \code{\link[conivol]{polyh_reduce_gen}} / \code{\link[conivol]{polyh_reduce_ineq}}:
#'                      compute a reduced representation
#'                      of a polyhedral cone given by generators / inequalities
#'   \item \code{\link[conivol]{polyh_rivols_gen}} / \code{\link[conivol]{polyh_rivols_ineq}}:
#'                      produce samples from the intrinsic volumes distribution of
#'                      a polyhedral cone given by generators / inequalities
#'   \item \code{\link[conivol]{polyh_rbichibarsq_gen}} / \code{\link[conivol]{polyh_rbichibarsq_ineq}}:
#'                      produce samples from the bivariate chi-bar-squared distribution
#'                      with weights given by the conic intrinsic volumes of
#'                      a polyhedral cone given by generators / inequalities
#'   \item \code{\link[conivol]{polyh_bayes}}: generates functions for computing quantiles of marginals
#'                      of the posterior distribution and for sampling
#'                      from the posterior distribution,
#'                      given samples of the intrinsic volumes distribution
#'                      (based on analytic solution)
#'   \item \code{\link[conivol]{polyh_stan}}: generates inputs for Stan
#'                      (data list and model string or external file)
#'                      for sampling from the posterior distribution,
#'                      given samples of the intrinsic volumes distribution
#'                      using a model that naturally implies log-concavity
#'                      (and cannot be solved analytically)
#' }
#'
#' @section Estimating the weights of the bivariate chi-bar-squared distribution:
#' \itemize{
#'   \item \code{\link[conivol]{estim_statdim_var}}: estimates the statistical dimension and
#'                      the variance of the intrinsic volumes from samples of the
#'                      corresponding bivariate chi-bar-squared distribution
#'   \item \code{\link[conivol]{init_ivols}}: find an initial estimate of the weights, potentially
#'                      based on first and/or second moment
#'   \item \code{\link[conivol]{loglike_ivols}}: compute the log-likelihood of a weight vector
#'                      for specific sample data
#'   \item \code{\link[conivol]{prepare_em}}: evaluates the sample data of the bivariate chi-bar-squared
#'                      data (find the corresponding chi-squared density values)
#'   \item \code{\link[conivol]{estim_em}}: produces EM-type iterates that may or may not converge
#'                      to the maximum likelihood estimate for the weights
#'                      of the bivariate chi-bar-squared distribution
#'                      from sample data
#'   \item \code{\link[conivol]{estim_jags}} / \code{\link[conivol]{estim_stan}}:
#'                      generates inputs for JAGS / Stan (data list and model
#'                      string or external file) for sampling from the posterior distribution
#'                      of the intrinsic volumes,
#'                      given samples of the bivariate chi-bar-squared distribution
#' }
#'
#' @section See also:
#' \describe{
#'   \item{github}{\href{https://github.com/damelunx/conivol}{github.com/damelunx/conivol}}
#'   \item{vignette}{\href{../doc/conic-intrinsic-volumes.html}{Conic intrinsic
#'         volumes and (bivariate) chi-bar-squared distribution}:
#'         introduces conic intrinsic volumes and (bivariate) chi-bar-squared distributions,
#'         as well as the computations involving polyhedral cones}
#'   \item{vignette}{\href{../doc/estim-conic-intrinsic-volumes-with-EM.html}{Estimating conic intrinsic volumes via EM algorithm}:
#'         describes the details of the algorithm for finding the intrinsic volumes of closed
#'         convex cones from samples of the associated bivariate chi-bar-squared distribution}
#'   \item{vignette}{\href{../doc/bayesian.html}{Bayesian estimates for conic intrinsic volumes}:
#'         describes the Bayesian approach for reconstructing intrinsic volumes
#'         from sampling data, which can either be samples from the intrinsic
#'         volumes distribution (in the case of polyhedral cones), or from the
#'         bivariate chi-bar-squared distribution, and which can be with or without
#'         enforcing log-concavity of the intrinsic volumes}
#' }
#'
#' @docType package
#' @name conivol
NULL














