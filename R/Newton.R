#' Finding the weights of the bivariate chi-bar-squared distribution using Newton's method.
#'
#' \code{bichibarsq_find_weights_Newton} produces Newton-type iterates from a two-column
#' matrix whose rows form iid samples from a bivariate chi-bar-squared
#' distribution, which may or may not (depending on the starting point) converge
#' to the maximum likelihood estimate of the mixing weights of the distribution.
#'
#' @param d the dimension of the bivariate chi-bar squared distribution.
#' @param m_samp two-column matrix whose rows from iid samples from a bivariate
#'               chi-bar-squared distribution.
#' @param N the number of iterates that shall be produced.
#' @param v_init the starting point for the EM iterates; if none are provided,
#'                the starting point is found in a way specified by the input \code{init_mode}.
#' @param init_mode specifies the way through which the initial estimate is found:
#'             \describe{
#'               \item{\code{init_mode==0}:}{uniform distribution}
#'               \item{\code{init_mode==1}:}{discretized normal distribution}
#'               \item{\code{init_mode==2}:}{circular cone fitting statdim}
#'               \item{\code{init_mode==3}:}{circular cone fitting variance}
#'               \item{\code{init_mode==4}:}{circular cone fitting statdim and variance (geometric mean)}
#'             }
#'             The starting point will be returned as the first row in the
#'             output matrix.
#' @param gamma positive weight to determine the step length of the iteration.
#' @param lambda nonnegative parameters which, if positive, enforce the
#'               log-concavity inequalities. Enforcing these may have negative
#'               effects on the performance, as the update step may become non-
#'               convex. \code{lambda} can be a scalar or vector of length \code{d-1}.
#' @param extrapolate specifies the way the edge cases are handled:
#'             \describe{
#'               \item{\code{extrapolate==0}:}{extrapolate \code{v_d} if no \code{x\in C} and
#'                                      extrapolate \code{v_0} if no \code{x\in C°}}
#'               \item{\code{extrapolate==1}:}{extrapolate \code{v_d}, do not extrapolate \code{v_0}}
#'               \item{\code{extrapolate==2}:}{extrapolate \code{v_0}, do not extrapolate \code{v_d}}
#'               \item{\code{extrapolate==3}:}{extrapolate \code{v_0} and \code{v_d}}
#'               \item{\code{extrapolate<0}:}{neither extrapolate \code{v_0} nor \code{v_d}}
#'             }
#' @param selfdual logical; if \code{TRUE}, the symmetry equations \code{v[k+1]==v[d-k+1]},
#'                 with \code{k=0,...,d}, are enforced. These equations hold for
#'                 the intrinsic volumes of self-dual cones.
#'
#' @param .data output of \code{bichibarsq_prepare_data(d, m_samp)}; this can be called
#'              outside and passed as input to avoid re-executing this
#'              potentially time-consuming step.
#'
#' @return The output of \code{bichibarsq_find_weights_Newton} is a \code{(N+1)}-by-\code{(d+1)}
#'         matrix whose rows constitute EM-type iterates, which may or may not
#'         converge to the maximum likelihood estimate of the mixing weights of
#'         the bivariate chi-bar-squared distribution.
#'
#' @examples
#' m_samp <- rbichibarsq_circ(10^6,c(5,5),c(pi/3,pi/4))
#' bichibarsq_find_weights_Newton( 10, m_samp )
#' bichibarsq_find_weights_Newton( 10, m_samp, init_mode=1 )
#'
bichibarsq_find_weights_Newton <- function(d, m_samp, N=20, v_init=NULL, init_mode=0,
                                       lambda=0, gamma=1, extrapolate=0, selfdual=FALSE, .data=NULL) {
    out <- matrix(0,N+1,d+1)

    # set the starting point for EM
    if (length(v_init)==d+1)
        v <- v_init
    else {
        est <- estimate_statdim_var(m_samp)
        v <- init_v(d,init_mode,delta=est$delta,var=est$var)
    }
    out[1, ] <- v

    # find the values of the chi-squared densities at the sample points
    if (is.null(.data))
        data <- bichibarsq_prepare_data(d, m_samp)

    # decide whether v0 or vd shall be extrapolated
    extrap_prim = (data$prop_prim==0 & extrapolate==0) | extrapolate==1 | extrapolate==3
    extrap_pol  = (data$prop_pol ==0 & extrapolate==0) | extrapolate==2 | extrapolate==3

    # prepare log-concavity enforcing parameters
    lambda_v <- c(0,rep_len(lambda, d-1),0)
    # prepare index sets for potential normalization
    I_even <- as.logical(rep_len(1:0,d+1))
    I_odd  <- as.logical(rep_len(0:1,d+1))

    if (!extrap_pol & !extrap_prim) {
        for (i in 1:N) {
            denom <- colSums( data$dens * v[2:d] )
            c <- rowSums( sweep( 1/data$n * data$dens * v[2:d] , MARGIN=2, denom, "/") ) +
                        lambda_v[1:(d-1)]-2*lambda_v[2:d]+lambda_v[3:(d+1)]
            H <- diag(lambda_v[1:(d-1)]-2*lambda_v[2:d]+lambda_v[3:(d+1)]) +
                    matrix( rowSums( sweep(
                        apply( 1/data$n * data$dens * v[2:d] , 2, function(x) return(x %x% x) ) ,
                        MARGIN=2, denom^2, "/") ), d-1, d-1 )

            v[1]   <- (1+gamma)*v[1]
            v[2:d] <- v[2:d] * (1 + gamma * solve(H,c) )
            v[d+1] <- (1+gamma)*v[d+1]
            v[I_even] <- 0.5 * v[I_even]/sum(v[I_even])
            v[I_odd]  <- 0.5 * v[I_odd] /sum(v[I_odd])
            out[i+1, ] <- v
        }
    } else if (extrap_pol & !extrap_prim) {
        for (i in 1:N) {
            denom <- colSums( data$dens * v[2:d] )
            c <- rowSums( sweep( 1/data$n * data$dens * v[2:d] , MARGIN=2, denom, "/") ) +
                lambda_v[1:(d-1)]-2*lambda_v[2:d]+lambda_v[3:(d+1)]
            H <- diag(lambda_v[1:(d-1)]-2*lambda_v[2:d]+lambda_v[3:(d+1)]) +
                matrix( rowSums( sweep(
                    apply( 1/data$n * data$dens * v[2:d] , 2, function(x) return(x %x% x) ) ,
                    MARGIN=2, denom^2, "/") ), d-1, d-1 )

            v[2:d] <- v[2:d] * (1 + gamma * solve(H,c) )
            v[d+1] <- (1+gamma)*v[d+1]
            v[1] <- exp( spline(x=1:d,y=log(v[2:(d+1)]),method="natural",xout=0)$y )
            v[I_even] <- 0.5 * v[I_even]/sum(v[I_even])
            v[I_odd]  <- 0.5 * v[I_odd] /sum(v[I_odd])
            out[i+1, ] <- v
        }
    } else if (!extrap_pol & extrap_prim) {
        for (i in 1:N) {
            denom <- colSums( data$dens * v[2:d] )
            c <- rowSums( sweep( 1/data$n * data$dens * v[2:d] , MARGIN=2, denom, "/") ) +
                lambda_v[1:(d-1)]-2*lambda_v[2:d]+lambda_v[3:(d+1)]
            H <- diag(lambda_v[1:(d-1)]-2*lambda_v[2:d]+lambda_v[3:(d+1)]) +
                matrix( rowSums( sweep(
                    apply( 1/data$n * data$dens * v[2:d] , 2, function(x) return(x %x% x) ) ,
                    MARGIN=2, denom^2, "/") ), d-1, d-1 )

            v[1]   <- (1+gamma)*v[1]
            v[2:d] <- v[2:d] * (1 + gamma * solve(H,c) )
            v[d+1] <- exp( spline(x=0:(d-1),y=log(v[1:d]),method="natural",xout=d)$y )
            v[I_even] <- 0.5 * v[I_even]/sum(v[I_even])
            v[I_odd]  <- 0.5 * v[I_odd] /sum(v[I_odd])
            out[i+1, ] <- v
        }
    } else if (extrap_pol & extrap_prim) {
        for (i in 1:N) {
            denom <- colSums( data$dens * v[2:d] )
            c <- rowSums( sweep( 1/data$n * data$dens * v[2:d] , MARGIN=2, denom, "/") ) +
                lambda_v[1:(d-1)]-2*lambda_v[2:d]+lambda_v[3:(d+1)]
            H <- diag(lambda_v[1:(d-1)]-2*lambda_v[2:d]+lambda_v[3:(d+1)]) +
                matrix( rowSums( sweep(
                    apply( 1/data$n * data$dens * v[2:d] , 2, function(x) return(x %x% x) ) ,
                    MARGIN=2, denom^2, "/") ), d-1, d-1 )

            v[2:d] <- v[2:d] * (1 + gamma * solve(H,c) )
            v[c(1,d+1)] <- exp( spline(x=1:(d-1),y=log(v[2:d]),method="natural",xout=c(0,d+1))$y )
            v[I_even] <- 0.5 * v[I_even]/sum(v[I_even])
            v[I_odd]  <- 0.5 * v[I_odd] /sum(v[I_odd])
            out[i+1, ] <- v
        }
    }

    return(out)
}
