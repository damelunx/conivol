#' The Chi-bar-squared Distribution
#'
#' Density, distribution function, and random generation for the chi-bar-squared
#' distribution with mixing weights (conic intrinsic volumes) \code{v}.
#'
#' @name Chibarsquare
#'
#' @section See also:
#' \code{\link[conivol]{dbichibarsq}}, \code{\link[conivol]{pbichibarsq}}, \code{\link[conivol]{rbichibarsq}},
#' \code{\link[stats]{dchisq}}, \code{\link[stats]{pchisq}}, \code{\link[stats]{qchisq}}, \code{\link[stats]{rchisq}}
#'
#' Package: \code{\link[conivol]{conivol}}
#'
NULL


#' @describeIn Chibarsquare Evaluates the density of the chi-bar-squared distribution.
#'
#' @param x,q vector of quantiles.
#' @param v vector of mixing weights (conic intrinsic volumes).
#'
#' @return The output of \code{dchibarsq(x, v)} is a vector of length \code{length(x)}
#'         such that the entries are the values of the density of the
#'         chi-bar-squared distribution with weights \code{v} in the
#'         corresponding entry of \code{x}.
#'
#' @examples
#' dchibarsq(3,c(0.1,0.6,0.3))
#' dchibarsq(c(0.2,3,4.5),c(0.1,0.6,0.3))
#'
#' @export
#'
dchibarsq <- function(x, v) {
    conivol:::.test_vector(v)
    out <- vector("double",length(x))
    d <- length(v)-1
    # check where x==0 and set to Inf if v[1]!=0 or v[2]!=0
    I <- which( sapply( x, function(t) isTRUE(all.equal(t,0)) ) )
    if (length(I)>0) {
        if (v[1]!=0 || v[2]!=0)
            out[I] <- Inf
        else
            out[I] <- v[3] * 0.5
        out[-I] <- colSums( sapply( x[-I], function(t) dchisq(t,1:d) ) * v[1:d] )
    } else {
        out <- colSums( sapply( x, function(t) dchisq(t,1:d) ) * v[1:d] )
    }
    return(out)
}


#' @describeIn Chibarsquare Evaluates the cdf of the chi-bar-squared distribution.
#'
#' @return The output of \code{pchibarsq(q, v)} is a vector of length \code{length(q)}
#'         such that the entries are the values of the cdf of the
#'         chi-bar-squared distribution with weights \code{v} in the
#'         corresponding entry of \code{q}.
#'
#' @examples
#' pchibarsq(3,c(0.1,0.6,0.3))
#' pchibarsq(c(0.2,3,4.5),c(0.1,0.6,0.3))
#'
#' @export
#'
pchibarsq <- function(q, v) {
    conivol:::.test_vector(v)
    out <- vector("double",length(q))
    d <- length(v)-1

    out[which(q>=0)] <- v[1] + colSums( sapply( q, function(t) pchisq(t,1:d) ) * v[1:d] )

    return(out)
}


#' @describeIn Chibarsquare Samples from the chi-bar-squared distribution.
#'
#' @param n number of observations.
#'
#' @return The output of \code{rchibarsq(n, v)} is a vector of length \code{n}
#'         such that the entries are iid samples from the
#'         chi-bar-squared distribution with weights \code{v}.
#'
#' @examples
#' rchibarsq(10,c(0.1,0.6,0.3))
#'
#' @export
#'
rchibarsq <- function(n, v) {
    conivol:::.test_vector(v)
    return( rchisq(n, sample(length(v), n, replace=TRUE, prob=v)-1) )
}


