#' The conic intrinsic volumes of (products of) circular cones.
#'
#' \code{circ_ivols} computes the conic intrinsic volumes of circular cones,
#' whose dimensions and angles are given in the vectors \code{d} and
#' \code{alpha} (vectors must be of same lengths); if the length of the vectors
#' is one, a single vector is returned; if the length of the vectors is greater
#' than one and \code{prcoduct==FALSE}, a list of vectors is returned;
#' if the length of the vectors is greater than one and \code{product==TRUE},
#' a single vector with the intrinsic volumes of the product cone is returned.
#'
#' @param d vector of dimensions; must be same length as \code{alpha}.
#' @param alpha vector of angles; must be same length as \code{d}.
#' @param product logical; if \code{TRUE}, intrinsic volumes of product cone are returned.
#'
#' @return If \code{length(d)==1} or \code{(length(d)>1 & product==TRUE)}
#'         then a single vectors will be returned.
#'         If \code{(length(d)>1 & product==FALSE)} then a list of
#'         vectors will be returned.
#'
#' @section See also:
#' \code{\link[conivol]{circ_rbichibarsq}}
#'
#' Package: \code{\link[conivol]{conivol}}
#'
#' @examples
#' circ_ivols(5, pi/4)
#' circ_ivols(c(5,5), c(pi/4,pi/8))
#' circ_ivols(c(5,5), c(pi/4,pi/8), product = TRUE)
#'
#' @export
#'
circ_ivols <- function(d, alpha, product = FALSE) {
    if (length(d)!=length(alpha))
        stop("Inputs d and alpha must be of same length.")

    V <- list()
    for (i in 1:length(d)) {
        if (d[i]<=1 || alpha[i]<0 || alpha[i]>pi/2)
            v <- NA
        else if (alpha[i]==0)
            v <- c(0.5, 0.5, rep(0,d[i]-2))
        else if (alpha[i]==pi/2)
            v <- c(rep(0,d[i]-2), 0.5, 0.5)
        else {
            v <- rep(0,d[i]+1)
            v[1] <- exp( lgamma(d[i]/2)-lgamma((d[i]+1)/2)-lgamma(1/2)+log(d[i]-1)-log(2)+
                             log( integrate(function(x){sin(x)^(d[i]-2)},0,pi/2-alpha[i])$value ) )

            v[d[i]+1] <- exp( lgamma(d[i]/2)-lgamma((d[i]+1)/2)-lgamma(1/2)+log(d[i]-1)-log(2)+
                                  log( integrate(function(x){sin(x)^(d[i]-2)},0,alpha[i])$value ) )

            k <- 1:(d[i]-1)
            v[2:d[i]] <- exp( lgamma(d[i]/2)-lgamma((k+1)/2)-lgamma((d[i]-k+1)/2)+
                                  (k-1)*log(sin(alpha[i]))+(d[i]-k-1)*log(cos(alpha[i]))-log(2) )
        }
        V[[i]] <- v
    }

    if (length(d)==1)
         return(V[[1]])
    else if (product)
         return(conivol::comp_ivols_product(V))
    else return(V)
}



#' Sample from bivariate chi-bar-squared distribution of products of circular cones.
#'
#' \code{circ_rbichibarsq} generates an \code{n} by \code{2} matrix
#' such that the rows form iid samples from the bivariate chi-bar-squared
#' distribution of products of circular cones.
#'
#' @param n number of samples
#' @param d dimension
#' @param alpha angle
#'
#' @return The output of \code{circ_rbichibarsq(n,d,alpha)} is a
#'         \code{n} by \code{2} matrix such that the rows form iid samples from the
#'         bivariate chi-bar-squared distribution with weights given by the
#'         intrinsic volumes of a product of circular cones corresponding
#'         to the vector of dimensions \code{d} and the vector of angles \code{alpha}.
#'
#' @section See also:
#' \code{\link[conivol]{circ_ivols}}, \code{\link[conivol]{polyh_rbichibarsq}},
#' \code{\link[conivol]{rbichibarsq}}, \code{\link[conivol]{rchibarsq}},
#' \code{\link[stats]{rchisq}},
#'
#' Package: \code{\link[conivol]{conivol}}
#'
#' @examples
#' circ_rbichibarsq(20,5,pi/3)
#' circ_rbichibarsq(20,c(5,3),c(pi/3,pi/4))
#'
#' @export
#'
circ_rbichibarsq <- function(n,d,alpha) {
    if (length(d)!=length(alpha))
        stop("Inputs d and alpha must be of same length.")
    if (any(d<2) || any(alpha<0) || any(alpha>pi/2))
        stop("Dimensions d must be >=2 and alpha must be between 0 and pi/2.")

    return( conivol::rbichibarsq(n, conivol::circ_ivols(d,alpha,TRUE)) )

    # SLOW:
    # if (length(d)==1)
    #     return( rbichibarsq(n, circ_ivols(d,alpha)) )
    # else
    #     return( Reduce('+', lapply(circ_ivols(d,alpha), function(v){rbichibarsq(n,v)}) ) )

    # SLOW:
    # out <- matrix(0,n,2)
    # for (i in 1:length(d)) {
    #     normsq_xbar <- rchisq(n,d[i]-1)
    #     normsq_xd <- rchisq(n,1)
    #     normsq_x <- normsq_xbar + normsq_xd
    #     sign_x <- sample(c(-1,1), n, replace=TRUE)
    #     ssq <- sin(alpha[i])^2
    #     csq <- cos(alpha[i])^2
    #
    #     # determine the points in C and update the output correspondingly
    #     I <- which(sign_x==1  & ssq*normsq_xd>=csq*normsq_xbar)
    #     out[I,1] <- out[I,1] + normsq_x[I]
    #
    #     # determine the points in C° and update the output correspondingly
    #     J <- which(sign_x==-1 & csq*normsq_xd<=ssq*normsq_xbar)
    #     out[J,2] <- out[J,2] + normsq_x[J]
    #
    #     # take care of the points not in C or C°
    #     proj <- normsq_x[-c(I,J)] *
    #         sin(alpha[i] + asin(sign_x[-c(I,J)]*sqrt(normsq_xd[-c(I,J)]/normsq_x[-c(I,J)])))^2
    #     out[-c(I,J),1] <- out[-c(I,J),1] + proj
    #     out[-c(I,J),2] <- out[-c(I,J),2] + normsq_x[-c(I,J)]-proj
    # }
    # return(out)
}


