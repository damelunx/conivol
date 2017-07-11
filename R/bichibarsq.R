#' The Bivariate Chi-bar-squared Distribution
#'
#' Density, distribution function, and random generation for the bivariate
#' chi-bar-squared distribution with mixing weights (conic intrinsic volumes) \code{v}.
#'
#' @name Bichibarsquare
NULL


#' @describeIn Bichibarsquare Evaluates the density of the bivariate chi-bar-squared distribution.
#'
#' @param x two-entry vectors or two-column matrices.
#' @param v vector of mixing weights (conic intrinsic volumes).
#'
#' @return The output of \code{dbichibarsq(x, v)} is either a number, if \code{x}
#'         is a vector of length \code{2}, or a vector of length \code{dim(x)[1]}
#'         such that the entries are the values of the density of the
#'         bivariate chi-bar-squared distribution with weights \code{v} in the
#'         corresponding entry (row) of \code{x}.
#'
#' @examples
#' dbichibarsq(c(2,3),c(0.1,0.6,0.3))
#' dbichibarsq(matrix(c(2,3,4,5),2,2),c(0.1,0.6,0.3))
#'
dbichibarsq <- function(x, v) {
    .test_vector(v)
    d <- length(v)-1
    D <- dim(x)
    # test whether input is vector or matrix; first cover vector case
    if ( is.null(D) ) {
        if (length(x)!=2)
            stop("Vector x must have length 2; or x can be two-column matrix.")
        # test whether point is on x-axis or y-axis
        if ( isTRUE(all.equal(x[2],0)) & ( v[1]>0 || v[2]>0 ) )
            return( Inf )
        else if ( isTRUE(all.equal(x[1],0)) & ( v[d+1]>0 || v[d]>0 ) )
            return( Inf )
        else
            return( sum( v * dchisq(x[1],0:d) * dchisq(x[2],rev(0:d)) ) )
    } else {
        if (D[2]!=2)
            stop("Matrix x must have two columns.")
        out <- vector("double",D[1])

        # there are quite a few poles in this density, so we need to be careful
        # check which points are in origin, only on x-axis, only on y-axis...

        onXaxis <- sapply( x[ ,2], function(t) isTRUE(all.equal(t,0)) )
        onYaxis <- sapply( x[ ,1], function(t) isTRUE(all.equal(t,0)) )
        I <- which( onXaxis & onYaxis )
        J <- which( onXaxis & !onYaxis )
        K <- which( !onXaxis & onYaxis )

        if ( length(I)>0 ) {
            if ( v[1]>0 || v[d+1]>0 || v[2]>0 || v[d]>0 )
                out[I] <- Inf
        }

        if ( length(J)>0 & d>0 ) {
            if ( d==1 ) {
                if ( v[2]>0 )
                    out[J] <- Inf
            } else if (d==2) {
                if ( v[3]>0 || v[2]>0 )
                    out[J] <- Inf
            } else {
                if ( v[d+1]>0 || v[d]>0 )
                    out[J] <- Inf
                else
                    out[J] <- colSums( sapply( x[J,1], function(t) dchisq(t,0:d-2) )
                                       * sapply( x[J,2], function(t) dchisq(t,rev(2:d)) )
                                       * v[1:d-1] )
            }
        }

        if ( length(K)>0 & d>0 ) {
            if ( d==1 ) {
                if ( v[1]>0 )
                    out[K] <- Inf
            } else if (d==2) {
                if ( v[1]>0 || v[2]>0 )
                    out[K] <- Inf
            } else {
                if ( v[1]>0 || v[2]>0 )
                    out[K] <- Inf
                else
                    out[K] <- colSums( sapply( x[K,1], function(t) dchisq(t,2:d) )
                                       * sapply( x[K,2], function(t) dchisq(t,rev(0:d-2)) )
                                       * v[3:d+1] )
            }
        }

        L <- c(I,J,K)
        if ( length(L)>0 )
            out[L] <- colSums( sapply( x[L,1], function(t) dchisq(t,0:d) )
                               * sapply( x[L,2], function(t) dchisq(t,rev(0:d)) )
                               * v )
        else
            out <- colSums( sapply( x[ ,1], function(t) dchisq(t,0:d) )
                            * sapply( x[ ,2], function(t) dchisq(t,rev(0:d)) )
                            * v )

        return(out)
    }
}


#' @describeIn Bichibarsquare Evaluates the cdf of the bivariate chi-bar-squared distribution.
#'
#' @return The output of \code{pbichibarsq(x, v)} is either a number, if \code{x}
#'         is a vector of length \code{2}, or a vector of length \code{dim(x)[1]}
#'         such that the entries are the values of the cdf of the
#'         bivariate chi-bar-squared distribution with weights \code{v} in the
#'         corresponding entry (row) of \code{x}.
#'
#'
#' @examples
#' pbichibarsq(3,c(0.1,0.6,0.3))
#' pbichibarsq(c(0.2,3,4.5),c(0.1,0.6,0.3))
#'
pbichibarsq <- function(x, v) {
    .test_vector(v)
    d <- length(v)-1
    D <- dim(x)
    if ( is.null(D) ) {
        if (length(x)!=2)
            stop("Vector x must have length 2; or x can be two-column matrix.")

        #careful, pchisq(0,0) gives 0 instead of 1 (!)
        return( v[1]*as.double(x[1]>=0)*pchisq(x[2],d) +
                    sum( v[2:d]*pchisq(x[1],2:d)*pchisq(x[2],rev(2:d)) ) +
                    v[d+1]*pchisq(x[1],d)*as.double(x[2]>=0) )
    } else {
        if (D[2]!=2)
            stop("Matrix x must have two columns.")
        out <- vector("double",D[1])

        # find out which rows are nonnegative
        I <- which( (x[1:D[1]]>=0) & (x[D[1]+(1:D[1])]>=0) )

        out[I] <- v[1]*pchisq(x[I,2],d) +
            colSums( sapply( x[I, ], function(t) pchisq(t[1],1:d-1)*pchisq(t[2],rev(1:d-1)) ) * v[2:d] ) +
            v[d+1]*pchisq(x[I,1],d)

        return(out)
    }
}


#' @describeIn Bichibarsquare Samples from the bivariate chi-bar-squared distribution.
#'
#' @param n number of observations.
#'
#' @return The output of \code{rbichibarsq(n, v)} is either a number, if \code{n==1},
#'         or a two-column matrix with \code{n} rows such that the rows are
#'         iid samples from the bivariate chi-bar-squared distribution with weights \code{v}.
#'
#' @examples
#' rbichibarsq(10,c(0.1,0.6,0.3))
#'
rbichibarsq <- function(n, v) {
    .test_vector(v)
    d <- length(v)-1
    Z <- sample(d+1, n, replace=TRUE, prob=v)-1
    return( matrix( c(rchisq(n,Z),rchisq(n,d-Z)), n, 2 ) )
}

