#' Compute intrinsic volumes of product cone from intrinsic volumes of components
#'
#' \code{prod_ivols} computes the intrinsic volumes of a product cone
#' from the intrinsic volumes of its components. That is, it returns the
#' convolution of the vectors of intrinsic volumes given in \code{V}.
#'
#' @param V list of vectors of intrinsic volumes
#'
#' @return The output of \code{prod_ivols(V)} is the vector that is
#'         obtained from convolving the components of \code{V}.
#'
#' @section See also:
#' Package: \code{\link[conivol]{conivol}}
#'
#' @examples
#' prod_ivols(list( c(0.5,0.5), c(0.1,0.4,0.5) ))
#' prod_ivols(list( c(0.5,0.5), c(0.5,0.5), c(0.5, 0.5) ))
#'
#' @export
#'
prod_ivols <- function(V) {
    lapply(V, conivol:::.test_vector)
    for (i in 2:length(V))
        V[[1]] <- convolve(V[[1]],rev(V[[i]]),type="o")
    # test which one are numerically zero, then set them equal to zero
    # (to avoid negative entries)
    I <- which(sapply(V[[1]], function(t){isTRUE(all.equal(t,0,tolerance=conivol:::.adj_tol))}))
    V[[1]][I]=0
    return(V[[1]])
}


#' Estimate the statistical dimension and variance from sample data
#'
#' \code{estim_statdim_var} estimates the statistical dimension and the
#' variance of intrinsic volumes from samples from the corresponding
#' bivariate chi-bar-squared distribution.
#'
#' @param d the dimension of the bivariate chi-bar squared distribution.
#' @param m_samp two-column matrix whose rows from iid samples from the bivariate
#'               chi-bar-squared distribution corresponding to the cone
#'
#' @return The output of \code{estim_statdim_var} is a two-element list
#'         consisting of the esimated statistical dimension \code{delta}
#'         and variance \code{var},
#'
#' @section See also:
#' \code{\link[conivol]{rbichibarsq}}, \code{\link[conivol]{ellips_rbichibarsq}},
#' \code{\link[conivol]{polyh_rbichibarsq_gen}}, \code{\link[conivol]{polyh_rbichibarsq_ineq}}
#'
#' Package: \code{\link[conivol]{conivol}}
#'
#' @examples
#' v <- circ_ivols(10,pi/8)
#' m_samp <- rbichibarsq(10^4,v)
#' estim_statdim_var(10, m_samp)
#'
#' # true values:
#' list( delta=sum((1:length(v)-1)*v),
#'       var=sum((1:length(v)-1)^2*v)-sum((1:length(v)-1)*v)^2 )
#'
#' @export
#'
estim_statdim_var <- function(d, m_samp) {
    md <- colMeans(m_samp)
    mv <- colMeans(m_samp^2)
    delta <- (md[1] + d-md[2])/2
    var <- sqrt( (1+mv[1]-(delta+1)^2) * (1+mv[2]-(d-delta+1)^2) )
    # var <- (1+mv[1]-(delta+1)^2)

    return(list(delta=delta,var=var))
}





