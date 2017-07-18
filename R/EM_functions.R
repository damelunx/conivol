#' Evaluate the sample data for maximum likelihood estimation.
#'
#' \code{bichibarsq_prepare_data} takes a two-column matrix whose rows form
#' iid samples from a bivariate chi-bar-squared distribution and
#' prepares the data used in maximum likelihood estimation.
#'
#' @param d the dimension of the bivariate chi-bar squared distribution.
#' @param m_samp two-column matrix whose rows from iid samples from a bivariate
#'               chi-bar-squared distribution.
#'
#' @return The output of \code{bichibarsq_prepare_data} is a list of four elements:
#'         \describe{
#'           \item{\code{n}:}{the number of overall sample points
#'                            (including those in primal or polar cone)}
#'           \item{\code{prop_prim}:}{proportion of points in primal cone}
#'           \item{\code{prop_pol}:}{proportion of points}
#'           \item{\code{dens}:}{the density values of the effective sample points
#'                               (neither in primal nor polar cone); that is,
#'                               \code{dens} is a \code{(d-1)} row matrix
#'                               such that the \code{k}th row contains the products
#'                               of the density values of the chi_k^2 and chi_(d-k)^2
#'                               distributions evaluated in the effective sample points;
#'                               the row-form of the matrix is more convenient for
#'                               the computations}
#'         }
#'
#' @examples
#' m_samp <- rbichibarsq_circ(10^6,c(5,5),c(pi/3,pi/4))
#' bichibarsq_prepare_data( 10, m_samp )
#'
bichibarsq_prepare_data <- function(d, m_samp) {
    I1 <- which( sapply(m_samp[ ,1], function(t){isTRUE(all.equal(t,0,tolerance=.adj_tol))}) )
    I2 <- which( sapply(m_samp[ ,2], function(t){isTRUE(all.equal(t,0,tolerance=.adj_tol))}) )

    n <- dim(m_samp)[1]

    out <- list()
    out$n         <- n
    out$prop_prim <- max(length(I2),1)/n
    out$prop_pol  <- max(length(I1),1)/n
    out$dens <- apply( m_samp[-c(I1,I2), ], 1, function(x){dchisq(x[1],1:(d-1))*dchisq(x[2],(d-1):1)} )

    return(out)
}


#' Finding an initial estimate of the intrinsic volumes.
#'
#' \code{init_v} find an initial estimate of the intrinsic volumes via
#' moment-fitting.
#'
#' @param d the dimension of the bivariate chi-bar squared distribution.
#' @param init_mode specifies the way through which the initial estimate is found:
#'             \describe{
#'               \item{\code{init_mode==0}:}{uniform distribution}
#'               \item{\code{init_mode==1}:}{discretized normal distribution}
#'               \item{\code{init_mode==2}:}{circular cone fitting statdim}
#'               \item{\code{init_mode==3}:}{circular cone fitting variance}
#'               \item{\code{init_mode==4}:}{circular cone fitting statdim and variance (geometric mean)}
#'             }
#' @param delta an estimate of the statistical dimension of the cone.
#' @param var an estimate of the variane of the intrinsic volumes.
#'
#' @return The output of \code{init_v} is a \code{(d+1)}-column vector.
#'
#' @examples
#' m_samp <- rbichibarsq_circ(10^6,c(5,5),c(pi/3,pi/4))
#' est <- estimate_statdim_var(m_samp)
#' init_v( 10 )
#' init_v( 10, 1, delta=est$delta, var=est$var )
#' init_v( 10, 2, delta=est$delta )
#' init_v( 10, 3, var=est$var )
#' init_v( 10, 4, delta=est$delta, var=est$var )
#'
# creating the starting point for the EM algorithm
#
init_v <- function(d,init_mode=0,delta=d/2,var=d/4) {
    if (init_mode==1) {
        v <- sapply( 0:d, function(k){pnorm((k+0.5-delta)/sqrt(var)) - pnorm((k-0.5-delta)/sqrt(var))})
        return(v/sum(v))
    } else if (init_mode==2) {
        alpha <- asin(sqrt(delta/d))
        return(circ_ivol(d,alpha))
    } else if (init_mode==3) {
        alpha <- asin(sqrt(2*var/(d-2)))/2
        if ((alpha<pi/4 && delta>d/2) || (alpha>pi/4 && delta<d/2))
            alpha <- pi/2-alpha;
        return(circ_ivol(d,alpha))
    } else if (init_mode==4) {
        alpha1 <- asin(sqrt(delta/d))
        alpha2 <- asin(sqrt(2*var/(d-2)))/2
        if ((alpha2<pi/4 && delta>d/2) || (alpha2>pi/4 && delta<d/2))
            alpha2 <- pi/2-alpha2;
        v1 <- circ_ivol(d,alpha1)
        v2 <- circ_ivol(d,alpha2)
        v <- sqrt(v1*v2)
        return(v/sum(v))
    } else {
        return(1/(d+1)*rep(1,d+1))
    }
}



#' Finding the weights of the bivariate chi-bar-squared distribution using EM algorithm.
#'
#' \code{bichibarsq_find_weights_EM} produces EM-type iterates from a two-column
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
#' @param lambda nonnegative parameters which, if positive, enforce the
#'               log-concavity inequalities. Enforcing these may have negative
#'               effects on the performance. \code{lambda} can be a scalar or vector of length \code{d-1}.
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
#' @return The output of \code{bichibarsq_find_weights_EM} is a \code{(N+1)}-by-\code{(d+1)}
#'         matrix whose rows constitute EM-type iterates, which may or may not
#'         converge to the maximum likelihood estimate of the mixing weights of
#'         the bivariate chi-bar-squared distribution.
#'
#' @examples
#' m_samp <- rbichibarsq_circ(10^6,c(5,5),c(pi/3,pi/4))
#' bichibarsq_find_weights_EM( 10, m_samp )
#' bichibarsq_find_weights_EM( 10, m_samp, init_mode=1 )
#'
bichibarsq_find_weights_EM <- function(d, m_samp, N=20, v_init=NULL, init_mode=0,
                                    lambda=0, extrapolate=0, selfdual=FALSE, .data=NULL) {
    #################################
    # some general internal options:
    # the number of times the algorithm shall try to enforce the log-concavity constraints
    NO_TRIES_LAMBDA <- 3
    # decide which MOSEK response codes shall be accepted
    check_mosek_success <- function(code)
        return( identical(code,0) | identical(code,10000) | identical(code,10001) )
    # set return messages from MOSEK (verbose=0 enforces silent mode)
    opts <- list(verbose=0)
    #################################

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

    # decide whether v0 or vd should/have to be extrapolated, add Machine epsilon
    # if extrapolation is necessary but prohibited
    extrap_prim = (data$prop_prim==0 & extrapolate==0) | extrapolate==1 | extrapolate==3
    extrap_pol  = (data$prop_pol ==0 & extrapolate==0) | extrapolate==2 | extrapolate==3
    if (data$prop_prim==0 & !extrap_prim) data$prop_prim <- .Machine$double.eps
    if (data$prop_pol ==0 & !extrap_pol)  data$prop_pol  <- .Machine$double.eps

    # prepare Mosek inputs
    mos_inp <- .create_mosek_input_EM(rep(0,d+1),extrap_pol,extrap_prim,selfdual)
    # prepare log-concavity enforcing parameters
    lambda_v <- c(0,0,rep_len(lambda, d-1),0,0)
    # prepare index sets for potential normalization
    I_even <- as.logical(rep_len(1:0,d+1))
    I_odd  <- as.logical(rep_len(0:1,d+1))

    if (!extrap_pol & !extrap_prim) {
        for (i in 1:N) {
            denom <- colSums( data$dens * v[2:d] )
            const_pre <- c( data$prop_pol,
                            rowSums( sweep( 1/data$n * data$dens * v[2:d] , MARGIN=2, denom, "/") ) ,
                            data$prop_prim )
            c_lambda <- (NO_TRIES_LAMBDA:0)/NO_TRIES_LAMBDA
            i_rel <- 0
            success <- FALSE
            while(!success & i_rel<=NO_TRIES_LAMBDA) {
                i_rel <- i_rel+1
                const <- const_pre + c_lambda[i_rel] * ( lambda_v[1:(d+1)]-2*lambda_v[2:(d+2)]+lambda_v[3:(d+3)] )
                mos_inp <- .update_mosek_input_EM(mos_inp,const)
                mos_out <- mosek(mos_inp, opts)
                success <- check_mosek_success(mos_out$response$code)
            }
            v <- mos_out$sol$itr$xx
            out[i+1, ] <- v
        }
    } else if (extrap_pol & !extrap_prim) {
        for (i in 1:N) {
            denom <- colSums( data$dens * v[2:d] )
            const_pre <- c( rowSums( sweep( 1/data$n * data$dens * v[2:d] , MARGIN=2, denom, "/") ) ,
                            data$prop_prim )
            c_lambda <- (NO_TRIES_LAMBDA:0)/NO_TRIES_LAMBDA
            i_rel <- 0
            success <- FALSE
            while(!success & i_rel<=NO_TRIES_LAMBDA) {
                i_rel <- i_rel+1
                const <- const_pre + c_lambda[i_rel] * ( lambda_v[2:(d+1)]-2*lambda_v[3:(d+2)]+lambda_v[4:(d+3)] )
                mos_inp <- .update_mosek_input_EM(mos_inp, const)
                mos_out <- mosek(mos_inp, opts)
                success <- check_mosek_success(mos_out$response$code)
            }
            v[2:(d+1)] <- mos_out$sol$itr$xx
            v[1] <- exp( spline(x=1:d,y=log(v[2:(d+1)]),method="natural",xout=0)$y )
            v[I_even] <- 0.5 * v[I_even]/sum(v[I_even])
            v[I_odd]  <- 0.5 * v[I_odd] /sum(v[I_odd])
            out[i+1, ] <- v
        }
    } else if (!extrap_pol & extrap_prim) {
        for (i in 1:N) {
            denom <- colSums( data$dens * v[2:d] )
            const_pre <- c( data$prop_pol,
                            rowSums( sweep( 1/data$n * data$dens * v[2:d] , MARGIN=2, denom, "/") ) )
            c_lambda <- (NO_TRIES_LAMBDA:0)/NO_TRIES_LAMBDA
            i_rel <- 0
            success <- FALSE
            while(!success & i_rel<=NO_TRIES_LAMBDA) {
                i_rel <- i_rel+1
                const <- const_pre + c_lambda[i_rel] * ( lambda_v[1:d]-2*lambda_v[2:(d+1)]+lambda_v[3:(d+2)] )
                mos_inp <- .update_mosek_input_EM(mos_inp,const)
                mos_out <- mosek(mos_inp, opts)
                success <- check_mosek_success(mos_out$response$code)
            }
            v[1:d] <- mos_out$sol$itr$xx
            v[d+1] <- exp( spline(x=0:(d-1),y=log(v[1:d]),method="natural",xout=d)$y )
            v[I_even] <- 0.5 * v[I_even]/sum(v[I_even])
            v[I_odd]  <- 0.5 * v[I_odd] /sum(v[I_odd])
            out[i+1, ] <- v
        }
    } else if (extrap_pol & extrap_prim) {
        for (i in 1:N) {
            denom <- colSums( data$dens * v[2:d] )
            const_pre <- rowSums( sweep( 1/data$n * data$dens * v[2:d] , MARGIN=2, denom, "/") )
            c_lambda <- (NO_TRIES_LAMBDA:0)/NO_TRIES_LAMBDA
            i_rel <- 0
            success <- FALSE
            while(!success & i_rel<=NO_TRIES_LAMBDA) {
                i_rel <- i_rel+1
                const <- const_pre + c_lambda[i_rel] * ( lambda_v[2:d]-2*lambda_v[3:(d+1)]+lambda_v[4:(d+2)] )
                mos_inp <- .update_mosek_input_EM(mos_inp,const)
                mos_out <- mosek(mos_inp, opts)
                success <- check_mosek_success(mos_out$response$code)
            }
            v[2:d] <- mos_out$sol$itr$xx
            v[c(1,d+1)] <- exp( spline(x=1:(d-1),y=log(v[2:d]),method="natural",xout=c(0,d+1))$y )
            v[I_even] <- 0.5 * v[I_even]/sum(v[I_even])
            v[I_odd]  <- 0.5 * v[I_odd] /sum(v[I_odd])
            out[i+1, ] <- v
        }
    }

    return(out)
}



