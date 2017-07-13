#' conivol: A package for the (bivariate) chi-bar-squared distribution and conic intrinsic volumes
#'
#' The conivol package [...]
#'
#'
#' @section Chi-bar-squared distribution:
#' \itemize{
#'   \item \code{dchibarsq}: evaluates the density
#'   \item \code{pchibarsq}: evaluates the cumulative distribution function
#'   \item \code{rchibarsq}: produces samples
#' }
#'
#' @section Bivariate chi-bar-squared distribution:
#' \itemize{
#'   \item \code{dbichibarsq}: evaluates the density
#'   \item \code{pbichibarsq}: evaluates the cumulative distribution function
#'   \item \code{rbichibarsq}: produces samples
#' }
#'
#' @section Computing with conic intrinsic volumes:
#' \itemize{
#'   \item \code{comp_ivols_product}: computes the intrinsic volumes of a product cone
#'                                    by convolving the intrinsic volumes of
#'                                    its elements
#'   \item \code{estimate_statdim_var}: estimates the statistical dimension and
#'                                      the variance of the intrinsic volumes
#'                                      from samples of the corresponding
#'                                      bivariate chi-bar-squared distribution
#' }
#'
#' @section Computations involving products of circular cones:
#' \itemize{
#'   \item \code{circ_ivol}: computes the intrinsic volumes of (a product of)
#'                           circular cones
#'   \item \code{rbichibarsq_circ}: produces samples from the bivariate
#'                                  chi-bar-squared distribution with weights
#'                                  given by the conic intrinsic volumes of
#'                                  a product of circular cones
#' }
#'
#' @section Estimating the weights of the bivariate chi-bar-squared distribution:
#' \itemize{
#'   \item \code{bichibarsq_find_weights}: produces EM-type iterates that may
#'                                         or may not converge to the maximum
#'                                         likelihood estimate for the weights
#'                                         of the bivariate chi-bar-squared distribution
#'                                         from sample data
#' }
#'
#' See the corresponding help functions for more information.
#'
#' @docType package
#' @name conivol
NULL














