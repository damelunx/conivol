#' Evaluate the sample data for maximum likelihood estimation
#'
#' \code{prepare_em} takes a two-column matrix whose rows form
#' iid samples from a bivariate chi-bar-squared distribution and
#' prepares the data used in maximum likelihood estimation.
#'
#' @param d the dimension of the bivariate chi-bar squared distribution.
#' @param m_samp two-column matrix whose rows from iid samples from a bivariate
#'               chi-bar-squared distribution.
#'
#' @return The output of \code{prepare_em} is a list of four elements:
#'         \itemize{
#'           \item \code{n}: the number of overall sample points
#'                            (including those in primal or polar cone)
#'           \item \code{prop_prim}: proportion of points in primal cone
#'           \item \code{prop_pol}: proportion of points
#'           \item \code{dens}: the density values of the effective sample points
#'                               (neither in primal nor polar cone); that is,
#'                               \code{dens} is a \code{(d-1)} row matrix
#'                               such that the \code{k}th row contains the products
#'                               of the density values of the chi_k^2 and chi_(d-k)^2
#'                               distributions evaluated in the effective sample points;
#'                               the row-form of the matrix is more convenient for
#'                               the computations
#'         }
#'
#' @section See also:
#' \code{\link[conivol]{rbichibarsq}}, \code{\link[conivol]{circ_rbichibarsq}},
#' \code{\link[conivol]{rbichibarsq_polyh}}, \code{\link[conivol]{loglike_ivols}},
#' \code{\link[conivol]{estim_em}}
#'
#' Package: \code{\link[conivol]{conivol}}
#'
#' @examples
#' D <- c(5,5)
#' alpha <- c(pi/3,pi/4)
#' d <- sum(D)
#' N <- 10^5
#' v_exact <- circ_ivols( D, alpha, product=TRUE )
#' m_samp <- rbichibarsq(N,v_exact)
#' prepare_em( d, m_samp )
#'
#' @export
#'
prepare_em <- function(d, m_samp) {
    I1 <- which( sapply(m_samp[ ,1], function(t){isTRUE(all.equal(t,0,tolerance=conivol:::.adj_tol))}) )
    I2 <- which( sapply(m_samp[ ,2], function(t){isTRUE(all.equal(t,0,tolerance=conivol:::.adj_tol))}) )

    n <- dim(m_samp)[1]

    out <- list()
    out$n         <- n
    out$prop_prim <- max(length(I2),1)/n
    out$prop_pol  <- max(length(I1),1)/n
    if (length(I1)+length(I2)>0)
        out$dens <- apply( m_samp[-c(I1,I2), ], 1, function(x){dchisq(x[1],1:(d-1))*dchisq(x[2],(d-1):1)} )
    else
        out$dens <- apply( m_samp, 1, function(x){dchisq(x[1],1:(d-1))*dchisq(x[2],(d-1):1)} )

    return(out)
}


#' Evaluate the log-likelihood of the estimated intrinsic volumes
#'
#' \code{loglike_ivols} evaluates the (normalized) log-likelihood of a vector
#' with respect to given data, the output of \code{prepare_em}.
#'
#' @param v vector of mixing weights (conic intrinsic volumes).
#' @param data output of \code{prepare_em(d, m_samp)}.
#' @param mode specifies whether the first and last values should be taken into account:
#'             \describe{
#'               \item{\code{mode==0}:}{take all into account}
#'               \item{\code{mode==1}:}{leave out the estimate of the dth intrinsic volume}
#'               \item{\code{mode==2}:}{leave out the estimate of the 0th intrinsic volume}
#'               \item{\code{mode==3}:}{leave out both estimates of the 0th and dth intrinsic volume}
#'             }
#'
#' @return The output of \code{loglike_ivols} is the value of the normalized
#'         log-likelihood of the mixing weights \code{v} with respect to the
#'         sample data given in \code{data}
#'
#' @section See also:
#' \code{\link[conivol]{prepare_em}}, \code{\link[conivol]{estim_statdim_var}}
#'
#' Package: \code{\link[conivol]{conivol}}
#'
#' @examples
#' D <- c(5,5)
#' alpha <- c(pi/3,pi/4)
#' d <- sum(D)
#' N <- 10^5
#' v_exact <- circ_ivols( D, alpha, product=TRUE )
#'
#' # collect sample data
#' m_samp <- rbichibarsq(N,v_exact)
#' data <- prepare_em(d, m_samp)
#' est <- estim_statdim_var(d, m_samp)
#' v_estim <- list(
#'     init0 = init_ivols( d ) ,
#'     init1 = init_ivols( d, 1, delta=est$delta, var=est$var ) ,
#'     init2 = init_ivols( d, 2, delta=est$delta ) ,
#'     init3 = init_ivols( d, 3, var=est$var ) ,
#'     init4 = init_ivols( d, 4, delta=est$delta, var=est$var ) )
#'
#' # evaluate log-likelihood function
#' loglike_ivols(v_exact, data)
#' loglike_ivols(v_estim$init0, data)
#' loglike_ivols(v_estim$init1, data)
#' loglike_ivols(v_estim$init2, data)
#' loglike_ivols(v_estim$init3, data)
#' loglike_ivols(v_estim$init4, data)
#'
#' @export
#'
loglike_ivols <- function(v, data, mode=0){
    conivol:::.test_vector(v)
    d <- length(v)-1
    if (dim(data$dens)[1]!=d-1)
        stop("Wrong format.")
    if (mode==1)
        return(data$prop_pol  * log(v[1]) +
                   sum( 1/data$n * log( colSums( data$dens * v[2:d] ) ) ) )
    else if (mode==2)
        return(sum( 1/data$n * log( colSums( data$dens * v[2:d] ) ) ) +
                   data$prop_prim * log(v[d+1]))
    else if (mode==3)
        return( sum( 1/data$n * log( colSums( data$dens * v[2:d] ) ) ) )
    else
        return(data$prop_pol  * log(v[1]) +
                   sum( 1/data$n * log( colSums( data$dens * v[2:d] ) ) ) +
                   data$prop_prim * log(v[d+1]) )
}



#' Finding an initial estimate of the intrinsic volumes
#'
#' \code{init_ivols} find an initial estimate of the intrinsic volumes via
#' moment-fitting.
#'
#' @param d the dimension of the bivariate chi-bar squared distribution.
#' @param init_mode specifies the way through which the initial estimate is found:
#'             \itemize{
#'               \item \code{init_mode==0}: uniform distribution
#'               \item \code{init_mode==1}: discretized normal distribution
#'               \item \code{init_mode==2}: circular cone fitting statdim
#'               \item \code{init_mode==3}: circular cone fitting variance
#'               \item \code{init_mode==4}: circular cone fitting statdim and variance (geometric mean)
#'             }
#' @param delta an estimate of the statistical dimension of the cone.
#' @param var an estimate of the variane of the intrinsic volumes.
#'
#' @return The output of \code{init_ivols} is a \code{(d+1)}-column vector.
#'
#' @section See also:
#' \code{\link[conivol]{rbichibarsq}}, \code{\link[conivol]{circ_rbichibarsq}},
#' \code{\link[conivol]{rbichibarsq_polyh}}, \code{\link[conivol]{loglike_ivols}}
#'
#' Package: \code{\link[conivol]{conivol}}
#'
#' @examples
#' D <- c(5,5)
#' d <- sum(D)
#' alpha <- c(pi/3,pi/4)
#' v_exact <- circ_ivols(D, alpha, product=TRUE)
#' m_samp <- rbichibarsq(10^6,v)
#' est <- estim_statdim_var(d, m_samp)
#'
#' list( v_exact = v_exact , v_init0 = init_ivols( d ) ,
#'       v_init1 = init_ivols( d, 1, delta=est$delta, var=est$var ) ,
#'       v_init2 = init_ivols( d, 2, delta=est$delta ) ,
#'       v_init3 = init_ivols( d, 3, delta=est$delta, var=est$var ) ,
#'       v_init4 = init_ivols( d, 4, delta=est$delta, var=est$var ) )
#'
#' @export
#'
init_ivols <- function(d,init_mode=0,delta=d/2,var=d/4) {
    if (init_mode==1) {
        v <- sapply( 0:d, function(k){pnorm((k+0.5-delta)/sqrt(var)) - pnorm((k-0.5-delta)/sqrt(var))})
        return(v/sum(v))
    } else if (init_mode==2) {
        alpha <- asin(min(1,sqrt(delta/d)))
        return(conivol::circ_ivols(d,alpha))
    } else if (init_mode==3) {
        alpha <- asin(min(1,sqrt(2*var/(d-2))))/2
        if ((alpha<pi/4 && delta>d/2) || (alpha>pi/4 && delta<d/2))
        alpha <- pi/2-alpha;
        return(conivol::circ_ivols(d,alpha))
    } else if (init_mode==4) {
        alpha1 <- asin(min(1,sqrt(delta/d)))
        alpha2 <- asin(min(1,sqrt(2*var/(d-2))))/2
        if ((alpha2<pi/4 && delta>d/2) || (alpha2>pi/4 && delta<d/2))
            alpha2 <- pi/2-alpha2;
        v1 <- conivol::circ_ivols(d,alpha1)
        v2 <- conivol::circ_ivols(d,alpha2)
        v <- sqrt(v1*v2)
        return(v/sum(v))
    } else {
        return(1/(d+1)*rep(1,d+1))
    }
}


# create the input for mosek for EM step
#
.create_mosek_input_em <- function(const,extrap_pol,extrap_prim,selfdual) {
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


.update_mosek_input_em <- function(mos_inp,const) {
    mos_inp$scopt$opro[3, ] <- const
    return(mos_inp)
}

#' Finding the weights of the bivariate chi-bar-squared distribution using EM algorithm
#'
#' \code{estim_em} produces EM-type iterates from a two-column
#' matrix whose rows form iid samples from a bivariate chi-bar-squared
#' distribution.
#'
#' The sequence of iterates may or may not converge
#' to the maximum likelihood estimate of the mixing weights of the distribution.
#' Log-concavity of the intrinsic volumes is enforced by projecting the logarithms
#' onto the cone of log-concave sequences; this can be turned off by setting
#' \code{no_of_lcc_projections=0}.
#'
#' @param d the dimension of the bivariate chi-bar squared distribution.
#' @param m_samp two-column matrix whose rows from iid samples from a bivariate
#'               chi-bar-squared distribution.
#' @param N the number of iterates that shall be produced.
#' @param v_init the starting point for the EM iterates; if none are provided,
#'                the starting point is found in a way specified by the input \code{init_mode}.
#' @param init_mode specifies the way through which the initial estimate is found:
#'             \itemize{
#'               \item \code{init_mode==0}: uniform distribution
#'               \item \code{init_mode==1}: discretized normal distribution
#'               \item \code{init_mode==2}: circular cone fitting statdim
#'               \item \code{init_mode==3}: circular cone fitting variance
#'               \item \code{init_mode==4}: circular cone fitting statdim and variance (geometric mean)
#'             }
#'             The starting point will be returned as the first row in the
#'             output matrix.
#' @param lambda nonnegative parameters which, if positive, enforce the
#'               log-concavity inequalities. Enforcing these may have negative
#'               effects on the performance. \code{lambda} can be a scalar or vector of length \code{d-1}.
#' @param no_of_lcc_projections number of projections on the log-concavity cone
#' @param lcc_amount constant for strict log-concavity
#' @param extrapolate specifies the way the edge cases are handled:
#'             \itemize{
#'               \item \code{extrapolate==0}: extrapolate \code{v_d} if no \code{x in C} and
#'                                      extrapolate \code{v_0} if no \code{x in C°}
#'               \item \code{extrapolate==1}: extrapolate \code{v_d}, do not extrapolate \code{v_0}
#'               \item \code{extrapolate==2}: extrapolate \code{v_0}, do not extrapolate \code{v_d}
#'               \item \code{extrapolate==3}: extrapolate \code{v_0} and \code{v_d}
#'               \item \code{extrapolate<0}: neither extrapolate \code{v_0} nor \code{v_d}
#'             }
#' @param selfdual logical; if \code{TRUE}, the symmetry equations \code{v[k+1]==v[d-k+1]},
#'                 with \code{k=0,...,d}, are enforced. These equations hold for
#'                 the intrinsic volumes of self-dual cones.
#' @param data output of \code{prepare_em(d, m_samp)}; this can be called
#'              outside and passed as input to avoid re-executing this
#'              potentially time-consuming step.
#'
#' @return The output of \code{estim_em} is a list of an \code{(N+1)}-by-\code{(d+1)}
#'         matrix whose rows constitute EM-type iterates, which may or may not
#'         converge to the maximum likelihood estimate of the mixing weights of
#'         the bivariate chi-bar-squared distribution, and the corresponding values
#'         of the log-likelihood function.
#'
#' @section See also:
#' \code{\link[conivol]{rbichibarsq}}, \code{\link[conivol]{circ_rbichibarsq}},
#' \code{\link[conivol]{rbichibarsq_polyh}}, \code{\link[conivol]{prepare_em}},
#' \code{\link[conivol]{init_ivols}}, \code{\link[conivol]{loglike_ivols}}
#'
#' Package: \code{\link[conivol]{conivol}}
#'
#' @examples
#' # define cone and find sample data
#' D <- c(5,5)
#' alpha <- c(pi/3,pi/4)
#' d <- sum(D)
#' N <- 10^5
#' v_exact <- circ_ivols( D, alpha, product=TRUE )
#' m_samp <- rbichibarsq(N,v_exact)
#'
#' # prepare data and run EM algorithm twice with different inits
#' data <- prepare_em( d, m_samp )
#' est1 <- estim_em( d, m_samp, data=data )
#' est2 <- estim_em( d, m_samp, init_mode=1, data=data )
#'
#' # plot the iterates of the first EM run
#' plot(1+0:d, v_exact)
#' lines(1+0:d, v_exact, col="red")
#' lines(1+0:d, est1$iterates[1,])
#' lines(1+0:d, est1$iterates[5,])
#' lines(1+0:d, est1$iterates[10,])
#' lines(1+0:d, est1$iterates[21,])
#'
#' # plot the iterates of the second EM run
#' plot(1+0:d, v_exact)
#' lines(1+0:d, v_exact, col="red")
#' lines(1+0:d, est2$iterates[1,])
#' lines(1+0:d, est2$iterates[5,])
#' lines(1+0:d, est2$iterates[10,])
#' lines(1+0:d, est2$iterates[21,])
#'
#' @export
#'
estim_em <- function(d, m_samp, N=20, v_init=NULL, init_mode=0,
                          lambda=0, no_of_lcc_projections=1, lcc_amount=0,
                          extrapolate=0, selfdual=FALSE, data=NULL) {
    if (!requireNamespace("Rmosek", quietly = TRUE))
        stop( paste0("\n Could not find package 'Rmosek'.",
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
        data <- conivol::prepare_em(d, m_samp)

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
        est <- conivol::estim_statdim_var(d, m_samp)
        v <- conivol::init_ivols(d,init_mode,delta=est$delta,var=est$var)
    }
    out_iterates[1, ] <- v
    out_loglike[1] <- loglike_ivols(v,data)

    # prepare Mosek inputs
    mos_inp <- .create_mosek_input_em(rep(0,d+1),extrap_pol,extrap_prim,selfdual)
    # prepare log-concavity enforcing parameters
    lambda_v <- c(0,0,rep_len(lambda, d-1),0,0)
    # prepare index sets for potential normalization
    I_even <- as.logical(rep_len(1:0,d+1))
    I_odd  <- as.logical(rep_len(0:1,d+1))

    # prepare Mosek inputs for log-concavity enforcing
    A_lcc <- matrix(0,d+1,d-1)
    diag(A_lcc) <- 1
    diag(A_lcc[2:d,]) <- -2
    diag(A_lcc[3:(d+1),]) <- 1
    # mos_inp_lcc <- conivol:::.create_mosek_input_polyh_prim(A_lcc, rep(0,d+1))
    mos_inp_lcc <- conivol:::.create_mosek_input_polyh_pol(A_lcc, rep(0,d+1), -lcc_amount)

    for (i in 1:N) {
        denom <- colSums( data$dens * v[2:d] )

        if (!extrap_pol & !extrap_prim)
            const_pre <- c( data$prop_pol,
                            rowSums( sweep( 1/data$n * data$dens * v[2:d] , MARGIN=2, denom, "/") ) ,
                            data$prop_prim )
        else if (extrap_pol & !extrap_prim)
            const_pre <- c( rowSums( sweep( 1/data$n * data$dens * v[2:d] , MARGIN=2, denom, "/") ) ,
                            data$prop_prim )
        else if (!extrap_pol & extrap_prim)
            const_pre <- c( data$prop_pol,
                            rowSums( sweep( 1/data$n * data$dens * v[2:d] , MARGIN=2, denom, "/") ) )
        else if (extrap_pol & extrap_prim)
            const_pre <- rowSums( sweep( 1/data$n * data$dens * v[2:d] , MARGIN=2, denom, "/") )

        c_lambda <- (NO_TRIES_LAMBDA:0)/NO_TRIES_LAMBDA
        i_rel <- 0
        success <- FALSE
        while(!success & i_rel<=NO_TRIES_LAMBDA) {
            i_rel <- i_rel+1

            if (!extrap_pol & !extrap_prim)
                const <- const_pre + c_lambda[i_rel] * ( 2*lambda_v[2:(d+2)]-lambda_v[1:(d+1)]-lambda_v[3:(d+3)] )
            else if (extrap_pol & !extrap_prim)
                const <- const_pre + c_lambda[i_rel] * ( 2*lambda_v[3:(d+2)]-lambda_v[2:(d+1)]-lambda_v[4:(d+3)] )
            else if (!extrap_pol & extrap_prim)
                const <- const_pre + c_lambda[i_rel] * ( 2*lambda_v[2:(d+1)]-lambda_v[1:d]-lambda_v[3:(d+2)] )
            else if (extrap_pol & extrap_prim)
                const <- const_pre + c_lambda[i_rel] * ( 2*lambda_v[3:(d+1)]-lambda_v[2:d]-lambda_v[4:(d+2)] )

            mos_inp <- .update_mosek_input_em(mos_inp,const)
            mos_out <- Rmosek::mosek(mos_inp, opts)
            success <- check_mosek_success(mos_out$response$code)
        }
        if (!extrap_pol & !extrap_prim)
            v <- mos_out$sol$itr$xx
        else if (extrap_pol & !extrap_prim) {
            v[2:(d+1)] <- mos_out$sol$itr$xx
            v[1] <- exp( spline(x=1:d,y=log(v[2:(d+1)]),method="natural",xout=0)$y )
            v[I_even] <- 0.5 * v[I_even]/sum(v[I_even])
            v[I_odd]  <- 0.5 * v[I_odd] /sum(v[I_odd])
        } else if (!extrap_pol & extrap_prim) {
            v[1:d] <- mos_out$sol$itr$xx
            v[d+1] <- exp( spline(x=0:(d-1),y=log(v[1:d]),method="natural",xout=d)$y )
            v[I_even] <- 0.5 * v[I_even]/sum(v[I_even])
            v[I_odd]  <- 0.5 * v[I_odd] /sum(v[I_odd])
        } else if (extrap_pol & extrap_prim) {
            v[2:d] <- mos_out$sol$itr$xx
            v[c(1,d+1)] <- exp( spline(x=1:(d-1),y=log(v[2:d]),method="natural",xout=c(0,d+1))$y )
            v[I_even] <- 0.5 * v[I_even]/sum(v[I_even])
            v[I_odd]  <- 0.5 * v[I_odd] /sum(v[I_odd])
        }

        i_lcc <- 0
        while (i_lcc < no_of_lcc_projections) {
            i_lcc <- i_lcc+1
            # mos_inp_lcc <- conivol:::.update_mosek_input_polyh_prim(mos_inp_lcc, log(v))
            # mos_out <- Rmosek::mosek(mos_inp_lcc, opts)
            # v <- v/exp(mos_out$sol$itr$xx[(d+2):(2*d+2)])
            mos_inp_lcc <- conivol:::.update_mosek_input_polyh_pol(mos_inp_lcc, log(v))
            mos_out <- Rmosek::mosek(mos_inp_lcc, opts)
            v <- exp(mos_out$sol$itr$xx[1:(d+1)])
            v[I_even] <- 0.5 * v[I_even]/sum(v[I_even])
            v[I_odd]  <- 0.5 * v[I_odd] /sum(v[I_odd])
        }

        out_iterates[i+1, ] <- v
        out_loglike[i+1] <- loglike_ivols(v,data)
    }

    return(list(iterates=out_iterates, loglike=out_loglike))
}


#' Finding the weights of the bivariate chi-bar-squared distribution using gradient descent
#'
#' \code{estim_gd} produces gradient descent iterates from a two-column
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
#' @param data output of \code{prepare_em(d, m_samp)}; this can be called
#'              outside and passed as input to avoid re-executing this
#'              potentially time-consuming step.
#'
#' @return The output of \code{estim_gd} is a list of an \code{(N+1)}-by-\code{(d+1)}
#'         matrix whose rows constitute gradient descent-type iterates, which may or may not
#'         converge to the maximum likelihood estimate of the mixing weights of
#'         the bivariate chi-bar-squared distribution, and the corresponding values
#'         of the log-likelihood function.
#'
#' @examples
#' m_samp <- circ_rbichibarsq(10^6,c(5,5),c(pi/3,pi/4))
#' estim_gd( 10, m_samp )
#' estim_gd( 10, m_samp, init_mode=1 )
#'
#' #@export #(gradient descent doesn't seem to be working well, so unless this is fixed, it should be not exported)
#'
.estim_gd <- function(d, m_samp, N=20, v_init=NULL, init_mode=0,
                          lambda=0, step_len=1, extrapolate=0, selfdual=FALSE, data=NULL) {
    # find the values of the chi-squared densities at the sample points
    if (is.null(data))
        data <- conivol::prepare_em(d, m_samp)

    out_iterates <- matrix(0,N+1,d+1)
    out_loglike  <- vector("double", N+1)

    # set the starting point for GD
    if (length(v_init)==d+1)
        v <- v_init
    else {
        est <- conivol::estim_statdim_var(d, m_samp)
        v <- conivol::init_ivols(d,init_mode,delta=est$delta,var=est$var)
    }
    out_iterates[1, ] <- v
    out_loglike[1] <- loglike_ivols(v,data)

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
        # out_loglike[i+1] <- loglike_ivols(v,data)
    }
    return(list(iterates=out_iterates, loglike=out_loglike))
}




#' Finding the weights of the bivariate chi-bar-squared distribution using Newton's method
#'
#' \code{estim_newton} produces Newton-type iterates from a two-column
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
#' @param data output of \code{prepare_em(d, m_samp)}; this can be called
#'              outside and passed as input to avoid re-executing this
#'              potentially time-consuming step.
#'
#' @return The output of \code{estim_newton} is a list of an \code{(N+1)}-by-\code{(d+1)}
#'         matrix whose rows constitute Newton-type iterates, which may or may not
#'         converge to the maximum likelihood estimate of the mixing weights of
#'         the bivariate chi-bar-squared distribution, and the corresponding values
#'         of the log-likelihood function.
#'
#' @examples
#' m_samp <- circ_rbichibarsq(10^6,c(5,5),c(pi/3,pi/4))
#' estim_newton( 10, m_samp )
#' estim_newton( 10, m_samp, init_mode=1 )
#'
#' #@export  #(Newton's method doesn't seem to be working well, so unless this is fixed, it should be not exported)
#'
.estim_newton <- function(d, m_samp, N=20, v_init=NULL, init_mode=0,
                              lambda=0, step_len=1, extrapolate=0, selfdual=FALSE, data=NULL) {
    # find the values of the chi-squared densities at the sample points
    if (is.null(data))
        data <- conivol::prepare_em(d, m_samp)

    out_iterates <- matrix(0,N+1,d+1)
    out_loglike  <- vector("double", N+1)

    # set the starting point for Newton
    if (length(v_init)==d+1)
        v <- v_init
    else {
        est <- conivol::estim_statdim_var(d, m_samp)
        v <- conivol::init_ivols(d,init_mode,delta=est$delta,var=est$var)
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
        out_loglike[i+1] <- loglike_ivols(v,data)
    }
    return(list(iterates=out_iterates, loglike=out_loglike))
}
