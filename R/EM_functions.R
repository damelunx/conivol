#' Finding the weights of the bivariate chi-bar-squared distribution.
#'
#' \code{bichibarsq_find_weights} produces EM-type iterates from a two-column
#' matrix whose rows form iid samples from a bivariate chi-bar-squared
#' distribution, which may or may not (depending on the starting point) converge
#' to the maximum likelihood estimate of the mixing weights of the distribution.
#'
#' @param d the dimension of the bivariate chi-bar squared distribution.
#' @param N the number of iterates that shall be produced.
#' @param m_samp two-column matrix whose rows from iid samples from a bivariate
#'               chi-bar-squared distribution.
#' @param v_start the starting point for the EM iterates; if none are provided,
#'                the starting point is found in a way specified by the input \code{mode}.
#' @param mode specifies the way through which the starting point is found:
#'             \describe{
#'               \item{\code{mode==0}:}{uniform distribution}
#'               \item{\code{mode==1}:}{discretized normal distribution}
#'               \item{\code{mode==2}:}{circular cone fitting statdim}
#'               \item{\code{mode==3}:}{circular cone fitting variance}
#'               \item{\code{mode==4}:}{circular cone fitting statdim and variance (geometric mean)}
#'             }
#'             The starting point will be returned as the first column in the
#'             output matrix.
#' @param lambda nonnegative parameters which, if positive, enforce the
#'               log-concavity inequalities. Enforcing these may have negative
#'               effects on the performance, as the update step may become non-
#'               convex. \code{lambda} can be a scalar or vector of length \code{d-1}.
#' @param selfdual logical; if \code{TRUE}, the symmetry equations \code{v[k+1]==v[d-k+1]},
#'                 with \code{k=0,...,d}, are enforced. These equations hold for
#'                 the intrinsic volumes of self-dual cones.
#'
#' @return The output of \code{bichibarsq_find_weights} is a \code{(d+1)}-by-\code{(N+1)}
#'         matrix whose columns constitute EM-type iterates, which may or may not
#'         converge to the maximum likelihood estimate of the mixing weights of
#'         the bivariate chi-bar-squared distribution.
#'
#' @examples
#' m_samp <- rbichibarsq_circ(10^6,c(5,5),c(pi/3,pi/4))
#' bichibarsq_find_weights( m_samp, 10 )
#' bichibarsq_find_weights( m_samp, 10, mode=1 )
#'
bichibarsq_find_weights <- function(m_samp, d, N=20, v_start=NULL, mode=0,
                                    lambda=0, selfdual=FALSE) {

    out <- matrix(0,d+1,N+1)

    # set the starting point for EM
    if (length(v_start)==d+1)
        v <- v_start
    else {
        # estimate statistical dimension
        md <- colMeans(m_samp)
        delta <- (md[1] + d-md[2])/2
        # estimate variance
        mv <- colMeans(m_samp^2)
        var <- sqrt( (1+mv[1]-(delta+1)^2) * (1+mv[2]-(d-delta+1)^2) )

        v <- .init_v(d,mode,delta=delta,var=var)
    }
    out[ ,1] <- v

    # find the values of the chi-squared densities at the sample points
    D <- .prepare_proj_data(d, m_samp)
    I0 <- which(D$ind==0)
    I1 <- which(D$ind==0 | D$ind==1)
    I2 <- which(D$ind==0 | D$ind==2)

    # prepare log-concavity enforcing parameters
    lambda_tmp <- rep_len(lambda, d-1)
    lambda0 <- c(0,0,lambda_tmp[2:(d-2)],0,0)
    lambda1 <- c(0,0,lambda_tmp[2:(d-1)],0,0)
    lambda2 <- c(0,0,lambda_tmp[1:(d-2)],0,0)

    # prepare Mosek inputs
    const <- list( c0=rep(0,d-1), c1=rep(0,d), c2=rep(0,d) )
    mos_inp <- .create_mosek_input(d,const,0,0,selfdual)
    opts <- list()
    opts$verbose <- 0

    # loop of EM algorithm
    for (i in 1:N) {
        # compute constants for next step
        v0 <- v[2:d]/(1-v[1]-v[d+1])
        v1 <- v[2:(d+1)]/(1-v[1])
        v2 <- v[1:d]/(1-v[d+1])

        denom0 <- rowSums( sweep( D$prim[I0,1:(d-1)] * D$pol[I0,rev(1:(d-1))] , MARGIN=2,v0,"*") )
        denom1 <- rowSums( sweep( D$prim[I1,1:d] , MARGIN=2,v1,"*") )
        denom2 <- rowSums( sweep( D$pol[I2,1:d]  , MARGIN=2,v2,"*") )

        c0 <- lambda0[1:(d-1)]-2*lambda0[2:d]+lambda0[3:(d+1)] +
            colSums( sweep( D$prim[I0,1:(d-1)] * D$pol[I0,rev(1:(d-1))] ,
                                                 MARGIN=2, v0, "*") / denom0 )
        c1 <- lambda1[1:d]-2*lambda1[2:(d+1)]+lambda1[3:(d+2)] +
            colSums( sweep( D$prim[I1,1:d] ,     MARGIN=2, v1, "*") / denom1 )
        c2 <- lambda2[1:d]-2*lambda2[2:(d+1)]+lambda2[3:(d+2)] +
            colSums( sweep( D$pol[I2,rev(1:d)] , MARGIN=2, v2, "*") / denom2 )

        # update mosek input
        mos_inp <- .update_mosek_input(mos_inp, d, c0, c1, c2, v[1], v[d+1])

        # solve maximization step
        w0 <- mosek(mos_inp$mode0, opts)$sol$itr$xx
        w1 <- mosek(mos_inp$mode1, opts)$sol$itr$xx
        w2 <- mosek(mos_inp$mode2, opts)$sol$itr$xx

        # set the next iterate
        v[1]   <- w2[1]  *(1-w1[d]) / (1-w2[1]*w1[d])
        v[d+1] <- w1[d]*(1-w2[1])   / (1-w2[1]*w1[d])
        v[2:d] <- w0 * (1-v[1]-v[d+1])

        out[ ,i+1] <- v
    }
    return(out)
}



