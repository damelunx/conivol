#' Evaluate the sample data for maximum likelihood estimation.
#'
#' \code{prepare_data} takes a two-column matrix whose rows form
#' iid samples from a bivariate chi-bar-squared distribution and
#' prepares the data used in maximum likelihood estimation.
#'
#' @param d the dimension of the bivariate chi-bar squared distribution.
#' @param m_samp two-column matrix whose rows from iid samples from a bivariate
#'               chi-bar-squared distribution.
#'
#' @return The output of \code{prepare_data} is a list of four elements:
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
#' prepare_data( 10, m_samp )
#'
#' @export
#'
prepare_data <- function(d, m_samp) {
    I1 <- which( sapply(m_samp[ ,1], function(t){isTRUE(all.equal(t,0,tolerance=conivol:::.conivol_adj_tol))}) )
    I2 <- which( sapply(m_samp[ ,2], function(t){isTRUE(all.equal(t,0,tolerance=conivol:::.conivol_adj_tol))}) )

    n <- dim(m_samp)[1]

    out <- list()
    out$n         <- n
    out$prop_prim <- max(length(I2),1)/n
    out$prop_pol  <- max(length(I1),1)/n
    out$dens <- apply( m_samp[-c(I1,I2), ], 1, function(x){dchisq(x[1],1:(d-1))*dchisq(x[2],(d-1):1)} )

    return(out)
}


#' Evaluate the log-likelihood of the estimated intrinsic volumes.
#'
#' \code{comp_loglike} evaluates the (normalized) log-likelihood of a vector
#' with respect to given data, the output of \code{prepare_data}.
#'
#' @param v vector of mixing weights (conic intrinsic volumes).
#' @param data output of \code{prepare_data(d, m_samp)}.
#'
#' @return The output of \code{comp_loglike} is the value of the normalized
#'         log-likelihood of the mixing weights \code{v} with respect to the
#'         sample data given in \code{data}
#'
#' @examples
#' D <- c(5,5)
#' alpha <- c(pi/3,pi/4)
#' d <- sum(D)
#' N <- 10^5
#' v_exact <- circ_ivol( D, alpha, product=TRUE )
#'
#' # collect sample data
#' m_samp <- rbichibarsq_circ(N,D,alpha)
#' data <- prepare_data(d, m_samp)
#' est <- estimate_statdim_var(d, m_samp)
#' v1 <- init_v( d )
#' v2 <- init_v( d, 1, delta=est$delta, var=est$var )
#' v3 <- init_v( d, 2, delta=est$delta )
#' v4 <- init_v( d, 3, var=est$var )
#' v5 <- init_v( d, 4, delta=est$delta, var=est$var )
#'
#' # evaluate log-likelihood function
#' comp_loglike(v_exact, data)
#' comp_loglike(v1, data)
#' comp_loglike(v2, data)
#' comp_loglike(v3, data)
#' comp_loglike(v4, data)
#' comp_loglike(v5, data)
#'
#' @export
#'
comp_loglike <- function(v, data){
    conivol:::.conivol_test_vector(v)
    d <- length(v)-1
    if (dim(data$dens)[1]!=d-1)
        stop("Wrong format.")
    return(
        data$prop_pol  * log(v[1]) +
            sum( 1/data$n * log( colSums( data$dens * v[2:d] ) ) ) +
            data$prop_prim * log(v[d+1])
    )
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
#' est <- estimate_statdim_var(d, m_samp)
#' init_v( 10 )
#' init_v( 10, 1, delta=est$delta, var=est$var )
#' init_v( 10, 2, delta=est$delta )
#' init_v( 10, 3, var=est$var )
#' init_v( 10, 4, delta=est$delta, var=est$var )
#'
# creating the starting point for the EM algorithm
#
#' @export
#'
init_v <- function(d,init_mode=0,delta=d/2,var=d/4) {
    if (init_mode==1) {
        v <- sapply( 0:d, function(k){pnorm((k+0.5-delta)/sqrt(var)) - pnorm((k-0.5-delta)/sqrt(var))})
        return(v/sum(v))
    } else if (init_mode==2) {
        alpha <- asin(sqrt(delta/d))
        return(conivol::circ_ivol(d,alpha))
    } else if (init_mode==3) {
        alpha <- asin(sqrt(2*var/(d-2)))/2
        if ((alpha<pi/4 && delta>d/2) || (alpha>pi/4 && delta<d/2))
        alpha <- pi/2-alpha;
        return(conivol::circ_ivol(d,alpha))
    } else if (init_mode==4) {
        alpha1 <- asin(sqrt(delta/d))
        alpha2 <- asin(sqrt(2*var/(d-2)))/2
        if ((alpha2<pi/4 && delta>d/2) || (alpha2>pi/4 && delta<d/2))
            alpha2 <- pi/2-alpha2;
        v1 <- conivol::circ_ivol(d,alpha1)
        v2 <- conivol::circ_ivol(d,alpha2)
        v <- sqrt(v1*v2)
        return(v/sum(v))
    } else {
        return(1/(d+1)*rep(1,d+1))
    }
}


# create the input for mosek for EM step
#
.conivol_create_mosek_input_EM <- function(const,extrap_pol,extrap_prim,selfdual) {
    d <- length(const)-1

    mos_inp <- list()
    mos_inp$sense <- "max"

    # setting optimizer
    if (!extrap_pol & !extrap_prim) {
        mos_inp$c <- c(rep(0,d+1))
        opro <- matrix(list(), nrow=5, ncol=d+1)
        rownames(opro) <- c("type","j","f","g","h")
        for (i in 1:(d+1)) {
            opro[ ,i] <- list("LOG", i, const[i], 1.0, 0.0)
        }
    } else if ((extrap_pol & !extrap_prim) | (!extrap_pol & extrap_prim)) {
        mos_inp$c <- c(rep(0,d))
        opro <- matrix(list(), nrow=5, ncol=d)
        rownames(opro) <- c("type","j","f","g","h")
        for (i in 1:d) {
            opro[ ,i] <- list("LOG", i, const[i], 1.0, 0.0)
        }
    } else if (extrap_pol & extrap_prim) {
        mos_inp$c <- c(rep(0,d-1))
        opro <- matrix(list(), nrow=5, ncol=d-1)
        rownames(opro) <- c("type","j","f","g","h")
        for (i in 1:(d-1)) {
            opro[ ,i] <- list("LOG", i, const[i], 1.0, 0.0)
        }
    }
    mos_inp$scopt <- list(opro=opro)

    # variable constraints
    if (!extrap_pol & !extrap_prim) {
        blx <- rep(0, d+1)            # v_i >= 0
        bux <- rep(.5, d+1)           # v_i <= 0.5
    } else if ((extrap_pol & !extrap_prim) | (!extrap_pol & extrap_prim)) {
        blx <- rep(0, d)              # v_i >= 0
        bux <- rep(.5, d)             # v_i <= 0.5
    } else if (extrap_pol & extrap_prim) {
        blx <- rep(0, d-1)            # v_i >= 0
        bux <- rep(.5, d-1)           # v_i <= 0.5
    }
    mos_inp$bx <- rbind(blx, bux)

    # constraint matrix:
    if (!selfdual) {
        nrow <- 2
    } else {
        if (!extrap_pol & !extrap_prim)
            nrow <- 2+floor((d+1)/2)
        else
            nrow <- 2+floor((d-1)/2)
    }
    if (!extrap_pol & !extrap_prim) {
        A <- Matrix::Matrix(c( rep_len(1:0,d+1), rep_len(0:1,d+1) , rep(0,(nrow-2)*(d+1)) ), nrow=nrow, byrow=TRUE, sparse=TRUE)
        if (nrow>2)
            for (i in 1:(nrow-2))
                A[2+i,c(i,d+2-i)] <- c(1,-1)
    } else if (extrap_pol & !extrap_prim) {
        A <- Matrix::Matrix(c( rep_len(1:0,d), rep_len(0:1,d) , rep(0,(nrow-2)*d) ), nrow=nrow, byrow=TRUE, sparse=TRUE)
        if (nrow>2)
            for (i in 1:(nrow-2))
                A[2+i,c(1+i,d+2-i)] <- c(1,-1)
    } else if (!extrap_pol & extrap_prim) {
        A <- Matrix::Matrix(c( rep_len(1:0,d), rep_len(0:1,d) , rep(0,(nrow-2)*d) ), nrow=nrow, byrow=TRUE, sparse=TRUE)
        if (nrow>2)
            for (i in 1:(nrow-2))
                A[2+i,c(i,d+2-i-1)] <- c(1,-1)
    } else if (extrap_pol & extrap_prim) {
        A <- Matrix::Matrix(c( rep_len(1:0,d-1), rep_len(0:1,d-1) , rep(0,(nrow-2)*(d-1)) ), nrow=nrow, byrow=TRUE, sparse=TRUE)
        if (nrow>2)
            for (i in 1:(nrow-2))
                A[2+i,c(i,d+2-i)] <- c(1,-1)
    }
    mos_inp$A <- A

    # constraint rhs:
    buc <- c( .5, .5, rep(0,nrow-2) )
    if (!extrap_pol & !extrap_prim)
        blc <- buc
    else if (extrap_pol & !extrap_prim)
        blc <- c( 0, .5, rep(0,nrow-2) )
    else if (!extrap_pol & extrap_prim) {
        if (d%%2 == 0)
            blc <- c( 0, .5, rep(0,nrow-2) )
        else
            blc <- c( .5, 0, rep(0,nrow-2) )
    } else if (extrap_pol & extrap_prim) {
        if (d%%2 == 0)
            blc <- c( 0, .5, rep(0,nrow-2) )
        else
            blc <- c( 0, 0, rep(0,nrow-2) )
    }

    mos_inp$bc <- rbind(blc, buc);

    return( mos_inp )
}


.conivol_update_mosek_input_EM <- function(mos_inp,const) {
    mos_inp$scopt$opro[3, ] <- const
    return(mos_inp)
}

#' Finding the weights of the bivariate chi-bar-squared distribution using EM algorithm.
#'
#' \code{find_ivols_EM} produces EM-type iterates from a two-column
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
#'               \item{\code{extrapolate==0}:}{extrapolate \code{v_d} if no \code{x in C} and
#'                                      extrapolate \code{v_0} if no \code{x in C°}}
#'               \item{\code{extrapolate==1}:}{extrapolate \code{v_d}, do not extrapolate \code{v_0}}
#'               \item{\code{extrapolate==2}:}{extrapolate \code{v_0}, do not extrapolate \code{v_d}}
#'               \item{\code{extrapolate==3}:}{extrapolate \code{v_0} and \code{v_d}}
#'               \item{\code{extrapolate<0}:}{neither extrapolate \code{v_0} nor \code{v_d}}
#'             }
#' @param selfdual logical; if \code{TRUE}, the symmetry equations \code{v[k+1]==v[d-k+1]},
#'                 with \code{k=0,...,d}, are enforced. These equations hold for
#'                 the intrinsic volumes of self-dual cones.
#'
#' @param data output of \code{prepare_data(d, m_samp)}; this can be called
#'              outside and passed as input to avoid re-executing this
#'              potentially time-consuming step.
#'
#' @return The output of \code{find_ivols_EM} is a list of an \code{(N+1)}-by-\code{(d+1)}
#'         matrix whose rows constitute EM-type iterates, which may or may not
#'         converge to the maximum likelihood estimate of the mixing weights of
#'         the bivariate chi-bar-squared distribution, and the corresponding values
#'         of the log-likelihood function.
#'
#' @examples
#' m_samp <- rbichibarsq_circ(10^6,c(5,5),c(pi/3,pi/4))
#' find_ivols_EM( 10, m_samp )
#' find_ivols_EM( 10, m_samp, init_mode=1 )
#'
#' @export
#'
find_ivols_EM <- function(d, m_samp, N=20, v_init=NULL, init_mode=0,
                          lambda=0, extrapolate=0, selfdual=FALSE, data=NULL) {
    if (!requireNamespace("Rmosek", quietly = TRUE))
        stop( paste0("\n Could not find package 'Rmosek'.",
            "\n If MOSEK is not available, try using 'find_ivols_GD' and 'find_ivols_Newton' instead of 'find_ivols_EM'.",
            "\n See the help entries for more information.") )
    if (!requireNamespace("Matrix", quietly = TRUE))
        stop("\n Could not find package 'Matrix'.")
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

    # find the values of the chi-squared densities at the sample points
    if (is.null(data))
        data <- conivol::prepare_data(d, m_samp)

    # decide whether v0 or vd should/have to be extrapolated, add Machine epsilon
    # if extrapolation is necessary but prohibited
    extrap_prim = (data$prop_prim==0 & extrapolate==0) | extrapolate==1 | extrapolate==3
    extrap_pol  = (data$prop_pol ==0 & extrapolate==0) | extrapolate==2 | extrapolate==3
    if (data$prop_prim==0 & !extrap_prim) data$prop_prim <- .Machine$double.eps
    if (data$prop_pol ==0 & !extrap_pol)  data$prop_pol  <- .Machine$double.eps

    out_iterates <- matrix(0,N+1,d+1)
    out_loglike  <- vector("double", N+1)

    # set the starting point for EM
    if (length(v_init)==d+1)
        v <- v_init
    else {
        est <- conivol::estimate_statdim_var(d, m_samp)
        v <- conivol::init_v(d,init_mode,delta=est$delta,var=est$var)
    }
    out_iterates[1, ] <- v
    out_loglike[1] <- comp_loglike(v,data)

    # prepare Mosek inputs
    mos_inp <- .conivol_create_mosek_input_EM(rep(0,d+1),extrap_pol,extrap_prim,selfdual)
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
                const <- const_pre + c_lambda[i_rel] * ( 2*lambda_v[2:(d+2)]-lambda_v[1:(d+1)]-lambda_v[3:(d+3)] )
                mos_inp <- .conivol_update_mosek_input_EM(mos_inp,const)
                mos_out <- Rmosek::mosek(mos_inp, opts)
                success <- check_mosek_success(mos_out$response$code)
            }
            v <- mos_out$sol$itr$xx
            out_iterates[i+1, ] <- v
            out_loglike[i+1] <- comp_loglike(v,data)
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
                const <- const_pre + c_lambda[i_rel] * ( 2*lambda_v[3:(d+2)]-lambda_v[2:(d+1)]-lambda_v[4:(d+3)] )
                mos_inp <- .conivol_update_mosek_input_EM(mos_inp, const)
                mos_out <- Rmosek::mosek(mos_inp, opts)
                success <- check_mosek_success(mos_out$response$code)
            }
            v[2:(d+1)] <- mos_out$sol$itr$xx
            v[1] <- exp( spline(x=1:d,y=log(v[2:(d+1)]),method="natural",xout=0)$y )
            v[I_even] <- 0.5 * v[I_even]/sum(v[I_even])
            v[I_odd]  <- 0.5 * v[I_odd] /sum(v[I_odd])
            out_iterates[i+1, ] <- v
            out_loglike[i+1] <- comp_loglike(v,data)
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
                const <- const_pre + c_lambda[i_rel] * ( 2*lambda_v[2:(d+1)]-lambda_v[1:d]-lambda_v[3:(d+2)] )
                mos_inp <- .conivol_update_mosek_input_EM(mos_inp,const)
                mos_out <- Rmosek::mosek(mos_inp, opts)
                success <- check_mosek_success(mos_out$response$code)
            }
            v[1:d] <- mos_out$sol$itr$xx
            v[d+1] <- exp( spline(x=0:(d-1),y=log(v[1:d]),method="natural",xout=d)$y )
            v[I_even] <- 0.5 * v[I_even]/sum(v[I_even])
            v[I_odd]  <- 0.5 * v[I_odd] /sum(v[I_odd])
            out_iterates[i+1, ] <- v
            out_loglike[i+1] <- comp_loglike(v,data)
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
                const <- const_pre + c_lambda[i_rel] * ( 2*lambda_v[3:(d+1)]-lambda_v[2:d]-lambda_v[4:(d+2)] )
                mos_inp <- .conivol_update_mosek_input_EM(mos_inp,const)
                mos_out <- Rmosek::mosek(mos_inp, opts)
                success <- check_mosek_success(mos_out$response$code)
            }
            v[2:d] <- mos_out$sol$itr$xx
            v[c(1,d+1)] <- exp( spline(x=1:(d-1),y=log(v[2:d]),method="natural",xout=c(0,d+1))$y )
            v[I_even] <- 0.5 * v[I_even]/sum(v[I_even])
            v[I_odd]  <- 0.5 * v[I_odd] /sum(v[I_odd])
            out[i+1, ] <- v
            out_loglike[i+1] <- comp_loglike(v,data)
        }
    }

    return(list(iterates=out_iterates, loglike=out_loglike))
}


#' Finding the weights of the bivariate chi-bar-squared distribution using gradient descent.
#'
#' \code{find_ivols_GD} produces gradient descent iterates from a two-column
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
#' @param step_len positive weight to determine the step length of the iteration.
#' @param lambda nonnegative parameters which, if positive, enforce the
#'               log-concavity inequalities. Enforcing these may have negative
#'               effects on the performance, as the update step may become non-
#'               convex. \code{lambda} can be a scalar or vector of length \code{d-1}.
#' @param extrapolate specifies the way the edge cases are handled:
#'             \describe{
#'               \item{\code{extrapolate==0}:}{extrapolate \code{v_d} if no \code{x in C} and
#'                                      extrapolate \code{v_0} if no \code{x in C°}}
#'               \item{\code{extrapolate==1}:}{extrapolate \code{v_d}, do not extrapolate \code{v_0}}
#'               \item{\code{extrapolate==2}:}{extrapolate \code{v_0}, do not extrapolate \code{v_d}}
#'               \item{\code{extrapolate==3}:}{extrapolate \code{v_0} and \code{v_d}}
#'               \item{\code{extrapolate<0}:}{neither extrapolate \code{v_0} nor \code{v_d}}
#'             }
#' @param selfdual logical; if \code{TRUE}, the symmetry equations \code{v[k+1]==v[d-k+1]},
#'                 with \code{k=0,...,d}, are enforced. These equations hold for
#'                 the intrinsic volumes of self-dual cones.
#'
#' @param data output of \code{prepare_data(d, m_samp)}; this can be called
#'              outside and passed as input to avoid re-executing this
#'              potentially time-consuming step.
#'
#' @return The output of \code{find_ivols_GD} is a list of an \code{(N+1)}-by-\code{(d+1)}
#'         matrix whose rows constitute gradient descent-type iterates, which may or may not
#'         converge to the maximum likelihood estimate of the mixing weights of
#'         the bivariate chi-bar-squared distribution, and the corresponding values
#'         of the log-likelihood function.
#'
#' @examples
#' m_samp <- rbichibarsq_circ(10^6,c(5,5),c(pi/3,pi/4))
#' find_ivols_GD( 10, m_samp )
#' find_ivols_GD( 10, m_samp, init_mode=1 )
#'
#' @export
#'
find_ivols_GD <- function(d, m_samp, N=20, v_init=NULL, init_mode=0,
                          lambda=0, step_len=1, extrapolate=0, selfdual=FALSE, data=NULL) {
    # find the values of the chi-squared densities at the sample points
    if (is.null(data))
        data <- conivol::prepare_data(d, m_samp)

    out_iterates <- matrix(0,N+1,d+1)
    out_loglike  <- vector("double", N+1)

    # set the starting point for GD
    if (length(v_init)==d+1)
        v <- v_init
    else {
        est <- conivol::estimate_statdim_var(d, m_samp)
        v <- conivol::init_v(d,init_mode,delta=est$delta,var=est$var)
    }
    out_iterates[1, ] <- v
    out_loglike[1] <- comp_loglike(v,data)

    # decide whether v0 or vd shall be extrapolated
    extrap_prim = (data$prop_prim==0 & extrapolate==0) | extrapolate==1 | extrapolate==3
    extrap_pol  = (data$prop_pol ==0 & extrapolate==0) | extrapolate==2 | extrapolate==3

    # prepare log-concavity enforcing parameters
    lambda_v     <- c(0,0,rep_len(lambda, d-1),0,0)
    # prepare index sets for potential normalization
    I_even <- as.logical(rep_len(1:0,d+1))
    I_odd  <- as.logical(rep_len(0:1,d+1))

    for (i in 1:N) {
        # choose k0 and k1
        v0 <- v
        v1 <- v
        v0[c(TRUE,as.logical(rep_len(1:0,d-1)),TRUE)] <- 0
        v1[c(TRUE,as.logical(rep_len(0:1,d-1)),TRUE)] <- 0
        k0 <- which.max(v0)-1
        k1 <- which.max(v1)-1
        # start finding gradient
        grad <- 2*lambda_v[2:(d+2)]-lambda_v[1:(d+1)]-lambda_v[3:(d+3)] -
            (2*lambda_v[k0+2]-lambda_v[k0+1]-lambda_v[k0+3]) * v/v[k0+1] * rep_len(1:0,d+1) -
            (2*lambda_v[k1+2]-lambda_v[k1+1]-lambda_v[k1+3]) * v/v[k1+1] * rep_len(0:1,d+1)
        # adding to the middle part
        denom <- colSums( data$dens * v[2:d] )
        num   <- 1/data$n * v[2:d] * (
            matrix(data$dens[k0, ] %x% rep_len(0:1,d-1), dim(data$dens)) +
                matrix(data$dens[k1, ] %x% rep_len(1:0,d-1), dim(data$dens)) -
                data$dens )
        grad[2:d] <- grad[2:d] - rowSums( sweep( num , MARGIN=2, denom, "/") )
        # adding to the outer part
        grad[1] <- grad[1] + data$prop_pol - sum( 1/data$n * v[1] * data$dens[k0, ] / denom )
        if (d%%2 == 0)
            grad[d+1] <- grad[d+1] + data$prop_prim - sum( 1/data$n * v[d+1] * data$dens[k0, ] / denom )
        else
            grad[d+1] <- grad[d+1] + data$prop_prim - sum( 1/data$n * v[d+1] * data$dens[k1, ] / denom )

        # gradient descent step
        I <- c(k0+1,k1+1)
        v[-I] <- v[-I]+step_len*v[-I]*grad[-I]
        v[k0+1] <- 0.5-sum( (v*rep_len(1:0,d+1))[-I] )
        v[k1+1] <- 0.5-sum( (v*rep_len(0:1,d+1))[-I] )

        if (extrap_pol & !extrap_prim) {
            v[1] <- exp( spline(x=1:d,y=log(v[2:(d+1)]),method="natural",xout=0)$y )
            v[I_even] <- 0.5 * v[I_even]/sum(v[I_even])
        } else if (!extrap_pol & extrap_prim) {
            v[d+1] <- exp( spline(x=0:(d-1),y=log(v[1:d]),method="natural",xout=d)$y )
            v[I_even] <- 0.5 * v[I_even]/sum(v[I_even])
            v[I_odd]  <- 0.5 * v[I_odd] /sum(v[I_odd])
        } else if (extrap_pol & extrap_prim) {
            v[c(1,d+1)] <- exp( spline(x=1:(d-1),y=log(v[2:d]),method="natural",xout=c(0,d+1))$y )
            v[I_even] <- 0.5 * v[I_even]/sum(v[I_even])
            v[I_odd]  <- 0.5 * v[I_odd] /sum(v[I_odd])
        }
        if (selfdual) {
            v <- (v+rev(v))/2
            v[I_even] <- 0.5 * v[I_even]/sum(v[I_even])
            v[I_odd]  <- 0.5 * v[I_odd] /sum(v[I_odd])
        }
        out_iterates[i+1, ] <- v
        # out_loglike[i+1] <- comp_loglike(v,data)
    }
    return(list(iterates=out_iterates, loglike=out_loglike))
}




#' Finding the weights of the bivariate chi-bar-squared distribution using Newton's method.
#'
#' \code{find_ivols_Newton} produces Newton-type iterates from a two-column
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
#' @param step_len positive weight to determine the step length of the iteration.
#' @param lambda nonnegative parameters which, if positive, enforce the
#'               log-concavity inequalities. Enforcing these may have negative
#'               effects on the performance, as the update step may become non-
#'               convex. \code{lambda} can be a scalar or vector of length \code{d-1}.
#' @param extrapolate specifies the way the edge cases are handled:
#'             \describe{
#'               \item{\code{extrapolate==0}:}{extrapolate \code{v_d} if no \code{x in C} and
#'                                      extrapolate \code{v_0} if no \code{x in C°}}
#'               \item{\code{extrapolate==1}:}{extrapolate \code{v_d}, do not extrapolate \code{v_0}}
#'               \item{\code{extrapolate==2}:}{extrapolate \code{v_0}, do not extrapolate \code{v_d}}
#'               \item{\code{extrapolate==3}:}{extrapolate \code{v_0} and \code{v_d}}
#'               \item{\code{extrapolate<0}:}{neither extrapolate \code{v_0} nor \code{v_d}}
#'             }
#' @param selfdual logical; if \code{TRUE}, the symmetry equations \code{v[k+1]==v[d-k+1]},
#'                 with \code{k=0,...,d}, are enforced. These equations hold for
#'                 the intrinsic volumes of self-dual cones.
#'
#' @param data output of \code{prepare_data(d, m_samp)}; this can be called
#'              outside and passed as input to avoid re-executing this
#'              potentially time-consuming step.
#'
#' @return The output of \code{find_ivols_Newton} is a list of an \code{(N+1)}-by-\code{(d+1)}
#'         matrix whose rows constitute Newton-type iterates, which may or may not
#'         converge to the maximum likelihood estimate of the mixing weights of
#'         the bivariate chi-bar-squared distribution, and the corresponding values
#'         of the log-likelihood function.
#'
#' @examples
#' m_samp <- rbichibarsq_circ(10^6,c(5,5),c(pi/3,pi/4))
#' find_ivols_Newton( 10, m_samp )
#' find_ivols_Newton( 10, m_samp, init_mode=1 )
#'
#' @export
#'
find_ivols_Newton <- function(d, m_samp, N=20, v_init=NULL, init_mode=0,
                              lambda=0, step_len=1, extrapolate=0, selfdual=FALSE, data=NULL) {
    # find the values of the chi-squared densities at the sample points
    if (is.null(data))
        data <- conivol::prepare_data(d, m_samp)

    out_iterates <- matrix(0,N+1,d+1)
    out_loglike  <- vector("double", N+1)

    # set the starting point for Newton
    if (length(v_init)==d+1)
        v <- v_init
    else {
        est <- conivol::estimate_statdim_var(d, m_samp)
        v <- conivol::init_v(d,init_mode,delta=est$delta,var=est$var)
    }
    out_iterates[1, ] <- v

    # decide whether v0 or vd shall be extrapolated
    extrap_prim = (data$prop_prim==0 & extrapolate==0) | extrapolate==1 | extrapolate==3
    extrap_pol  = (data$prop_pol ==0 & extrapolate==0) | extrapolate==2 | extrapolate==3

    # prepare log-concavity enforcing parameters
    lambda_v <- c(0,0,rep_len(lambda, d-1),0,0)
    # prepare index sets for potential normalization
    I_even <- as.logical(rep_len(1:0,d+1))
    I_odd  <- as.logical(rep_len(0:1,d+1))

    for (i in 1:N) {
        # choose k0 and k1
        v0 <- v
        v1 <- v
        v0[c(TRUE,as.logical(rep_len(1:0,d-1)),TRUE)] <- 0
        v1[c(TRUE,as.logical(rep_len(0:1,d-1)),TRUE)] <- 0
        k0 <- which.max(v0)-1
        k1 <- which.max(v1)-1
        # start finding gradient
        grad <- 2*lambda_v[2:(d+2)]-lambda_v[1:(d+1)]-lambda_v[3:(d+3)] -
            (2*lambda_v[k0+2]-lambda_v[k0+1]-lambda_v[k0+3]) * v/v[k0+1] * rep_len(1:0,d+1) -
            (2*lambda_v[k1+2]-lambda_v[k1+1]-lambda_v[k1+3]) * v/v[k1+1] * rep_len(0:1,d+1)
        # adding to the middle part
        denom <- colSums( data$dens * v[2:d] )
        num   <- 1/data$n * v[2:d] * (
            matrix(data$dens[k0, ] %x% rep_len(0:1,d-1), dim(data$dens)) +
                matrix(data$dens[k1, ] %x% rep_len(1:0,d-1), dim(data$dens)) -
                data$dens
        )
        grad[2:d] <- grad[2:d] - rowSums( sweep( num , MARGIN=2, denom, "/") )
        # adding to the outer part
        grad[1] <- grad[1] + data$prop_pol - sum( 1/data$n * v[1] * data$dens[k0, ] / denom )
        if (d%%2 == 0)
            grad[d+1] <- grad[d+1] + data$prop_prim - sum( 1/data$n * v[d+1] * data$dens[k0, ] / denom )
        else
            grad[d+1] <- grad[d+1] + data$prop_prim - sum( 1/data$n * v[d+1] * data$dens[k1, ] / denom )

        # computing the negative Hessian
        negHess <- diag( 2*lambda_v[2:(d+2)]-lambda_v[1:(d+1)]-lambda_v[3:(d+3)] ) +
                    (2*lambda_v[k0+2]-lambda_v[k0+1]-lambda_v[k0+3]) *
                    matrix( (v/v[k0+1]*rep_len(1:0,d+1)) %x% (v/v[k0+1]*rep_len(1:0,d+1)), d+1, d+1) +
                    (2*lambda_v[k1+2]-lambda_v[k1+1]-lambda_v[k1+3]) *
                    matrix( (v/v[k1+1]*rep_len(0:1,d+1)) %x% (v/v[k1+1]*rep_len(0:1,d+1)), d+1, d+1)
        # adding corner values
        negHess[1,1]     <- negHess[1,1]     + data$prop_pol
        negHess[d+1,d+1] <- negHess[d+1,d+1] + data$prop_prim
        # adding main middle part
        negHess[2:d,2:d] <- negHess[2:d,2:d] + matrix( rowSums( sweep( 1/data$n *
            apply( v[2:d] * (matrix(data$dens[k0, ] %x% rep_len(0:1,d-1), dim(data$dens)) +
                             matrix(data$dens[k1, ] %x% rep_len(1:0,d-1), dim(data$dens)) - data$dens) , 2,
                   function(x) return(x %x% x) ) ,
            MARGIN=2, denom^2, "/") ), d-1, d-1 )
        # adding to first and last row
        negHess[1, ] <- negHess[1, ] + 1/data$n * rowSums( sweep( 1/data$n *
            v * (matrix(data$dens[k0, ] %x% rep_len(1:0,d+1), dim(data$dens)+c(2,0)) +
                 matrix(data$dens[k1, ] %x% rep_len(0:1,d+1), dim(data$dens)+c(2,0)) ) ,
            MARGIN=2, denom^2/(v[1]*data$dens[k0, ]), "/") )
        if (d%%2==0)
            negHess[d+1, ] <- negHess[d+1, ] + 1/data$n * rowSums( sweep( 1/data$n *
                v * (matrix(data$dens[k0, ] %x% rep_len(1:0,d+1), dim(data$dens)+c(2,0)) +
                     matrix(data$dens[k1, ] %x% rep_len(0:1,d+1), dim(data$dens)+c(2,0)) ) ,
                MARGIN=2, denom^2/(v[d+1]*data$dens[k0, ]), "/") )
        else
            negHess[d+1, ] <- negHess[d+1, ] + 1/data$n * rowSums( sweep( 1/data$n *
                v * (matrix(data$dens[k0, ] %x% rep_len(1:0,d+1), dim(data$dens)+c(2,0)) +
                     matrix(data$dens[k1, ] %x% rep_len(0:1,d+1), dim(data$dens)+c(2,0)) ) ,
                MARGIN=2, denom^2/(v[d+1]*data$dens[k1, ]), "/") )
        # adding to first and last column
        negHess[2:d, 1] <- negHess[2:d, 1] + 1/data$n * rowSums( sweep( 1/data$n *
            v[2:d] * (matrix(data$dens[k0, ] %x% rep_len(0:1,d-1), dim(data$dens)) +
                      matrix(data$dens[k1, ] %x% rep_len(1:0,d-1), dim(data$dens)) ) ,
            MARGIN=2, denom^2/(v[1]*data$dens[k0, ]), "/") )
        if (d%%2==0)
            negHess[2:d, d+1] <- negHess[2:d, d+1] + 1/data$n * rowSums( sweep( 1/data$n *
                v[2:d] * (matrix(data$dens[k0, ] %x% rep_len(0:1,d-1), dim(data$dens)) +
                          matrix(data$dens[k1, ] %x% rep_len(1:0,d-1), dim(data$dens)) ) ,
                MARGIN=2, denom^2/(v[d+1]*data$dens[k0, ]), "/") )
        else
            negHess[2:d, d+1] <- negHess[2:d, d+1] + 1/data$n * rowSums( sweep( 1/data$n *
                v[2:d] * (matrix(data$dens[k0, ] %x% rep_len(0:1,d-1), dim(data$dens)) +
                          matrix(data$dens[k1, ] %x% rep_len(1:0,d-1), dim(data$dens)) ) ,
                MARGIN=2, denom^2/(v[d+1]*data$dens[k1, ]), "/") )

        # compute Newton step
        I <- c(k0+1,k1+1)
        v[-I] <- v[-I] + step_len*v[-I]*solve( negHess[-I,-I], grad[-I] )
        v[k0+1] <- 0.5-sum( (v*rep_len(1:0,d+1))[-I] )
        v[k1+1] <- 0.5-sum( (v*rep_len(0:1,d+1))[-I] )

        if (extrap_pol & !extrap_prim) {
            v[1] <- exp( spline(x=1:d,y=log(v[2:(d+1)]),method="natural",xout=0)$y )
            v[I_even] <- 0.5 * v[I_even]/sum(v[I_even])
        } else if (!extrap_pol & extrap_prim) {
            v[d+1] <- exp( spline(x=0:(d-1),y=log(v[1:d]),method="natural",xout=d)$y )
            v[I_even] <- 0.5 * v[I_even]/sum(v[I_even])
            v[I_odd]  <- 0.5 * v[I_odd] /sum(v[I_odd])
        } else if (extrap_pol & extrap_prim) {
            v[c(1,d+1)] <- exp( spline(x=1:(d-1),y=log(v[2:d]),method="natural",xout=c(0,d+1))$y )
            v[I_even] <- 0.5 * v[I_even]/sum(v[I_even])
            v[I_odd]  <- 0.5 * v[I_odd] /sum(v[I_odd])
        }
        if (selfdual) {
            v <- (v+rev(v))/2
            v[I_even] <- 0.5 * v[I_even]/sum(v[I_even])
            v[I_odd]  <- 0.5 * v[I_odd] /sum(v[I_odd])
        }
        out_iterates[i+1, ] <- v
        out_loglike[i+1] <- comp_loglike(v,data)
    }
    return(list(iterates=out_iterates, loglike=out_loglike))
}
