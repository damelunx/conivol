# create the mosek input for the projection on {Ax|x>=0}
#
.create_mosek_input_polyh_prim <- function(A,z) {
    m <- dim(A)[1]
    n <- dim(A)[2]

    Aext <- cbind( A, matrix(0,m,2), diag(-1,m) )

    mos_inp <- list(sense = "min")
    mos_inp$c     <- c( -as.vector(t(A)%*%z), 1, rep(0,m+1) )
    mos_inp$A     <- Matrix::Matrix( Aext, sparse=TRUE )
    mos_inp$bc    <- rbind(blc = rep(0,m),
                           buc = rep(0,m))
    mos_inp$bx    <- rbind(blx = c(rep(0,n+1),1,rep(-Inf,m)),
                           bux = c(rep(Inf,n+1),1,rep(Inf,m)))
    mos_inp$cones <- matrix( list("RQUAD", c(n+1,n+2,(n+3):(m+n+2))), 2, 1 )
    rownames(mos_inp$cones) <- c("type", "sub")

    return(mos_inp)
}

.update_mosek_input_polyh_prim <- function(mos_inp,z) {
    m <- length(z)
    n <- dim(mos_inp$A)[2]-m-2
    mos_inp$c[1:n] <- -as.vector( Matrix::t(mos_inp$A[, 1:n]) %*% z )
    return(mos_inp)
}

# create the mosek input for the projection on {y|A^Ty<=c}
#
.create_mosek_input_polyh_pol <- function(A,z,c=0) {
    m <- dim(A)[1]
    n <- dim(A)[2]

    Aext <- cbind( t(A), matrix(0,n,2), diag(1,n) )

    mos_inp <- list(sense = "min")
    mos_inp$c     <- c( -z, 1, rep(0,n+1) )
    mos_inp$A     <- Matrix::Matrix( Aext, sparse=TRUE )
    mos_inp$bc    <- rbind(blc = rep_len(c,n),
                           buc = rep_len(c,n))
    mos_inp$bx    <- rbind(blx = c(rep(-Inf,m),0,1,rep(0,n)),
                           bux = c(rep(Inf,m+1),1,rep(Inf,n)))
    mos_inp$cones <- matrix( list("RQUAD", c(m+1,m+2,1:m)), 2, 1 )
    rownames(mos_inp$cones) <- c("type", "sub")

    return(mos_inp)
}

.update_mosek_input_polyh_pol <- function(mos_inp,z) {
    m <- length(z)
    mos_inp$c[1:m] <- -z
    return(mos_inp)
}



#' Find a reduced form of a polyhedral cone given by generators
#'
#' \code{polyh_reduce} takes as input a \code{n} by \code{m} matrix \code{A}
#' and returns a reduced form described by orthogonal bases for lineality space
#' and linear span, as well as a matrix generating the reduced cone.
#' See below for more details.
#'
#' @param A matrix
#' @param solver either "nnls" or "mosek"
#' @param tol tolerance (single precision machine epsilon by default)
#'
#' @return The output of \code{polyh_reduce(A)} is a list containing the following elements:
#' \itemize{
#'   \item \code{dimC}: the dimension of the linear span of \code{C},
#'   \item \code{linC}: the lineality of the cone \code{C},
#'   \item \code{QL}: an orthogonal basis of the lineality space of \code{C},
#'                    set to \code{NA} if lineality space is zero-dimensional,
#'   \item \code{QC}: an orthogonal basis of the projection of \code{C} onto
#'                    the orthogonal complement of the lineality space of \code{C},
#'                    set to \code{NA} if \code{C} is a linear space,
#'   \item \code{A_reduced}: a matrix defining the reduced cone.
#' }
#'
#' @section See also:
#'
#' Package: \code{\link[conivol]{conivol}}
#'
#' @note See \href{../doc/conic-intrinsic-volumes.html#sampling_polyh}{this vignette}
#'       for further info.
#'
#' @examples
#' A <- cbind(diag(1,3),diag(1,3)+matrix(1,ncol=3,nrow=3),c(-1,0,0))
#' A <- A[,sample(ncol(A))]
#' print(A)
#'
#' A_red <- polyh_reduce(A)$A_reduced
#' print(A_red)
#'
#' A_red <- polyh_reduce(A, solver="mosek")$A_reduced
#' print(A_red)
#'
#' @export
#'
polyh_reduce <- function(A, solver="nnls", tol=1e-7) {
    if (solver=="nnls") {
        if (!requireNamespace("nnls", quietly = TRUE))
            stop("\n Could not find package 'nnls'.")
    } else if (solver=="mosek") {
        if (!requireNamespace("Rmosek", quietly = TRUE))
            stop("\n Could not find package 'Rmosek'.")
        if (!requireNamespace("Matrix", quietly = TRUE))
            stop("\n Could not find package 'Matrix'.")
        opts <- list(verbose=0)
    } else
        stop("\n Parameter solver must have value either 'nnls' or 'mosek'.")
    if (!is.matrix(A))
        stop("\n Input is not a matrix.")

    d <- dim(A)[1]
    dimC <- qr(A)$rank
    if (dimC==0)
        return( list( dimC=0, linC=0, QL=NA, QC=NA, A_reduced=0) )
    # find lineality space L
    nn <- dim(A)[2]
    is_in_L <- vector("logical",nn)
    if (solver=="nnls") {
        for (j in 1:nn) {
            is_in_L[j] <- isTRUE(nnls::nnls(A[,-j],-A[,j])$deviance<d*tol)
        }
    } else {
        for (j in 1:nn) {
            mos_inp <- conivol:::.create_mosek_input_polyh_prim(A[ ,-j], -A[ ,j])
            is_in_L[j] <- isTRUE( sum(
                (Rmosek::mosek(mos_inp,opts)$sol$itr$xx[(nn+2):(nn+d+1)]+A[ ,j])^2 ) < d*tol )
        }
    }
    A_L <- A[,is_in_L]
    linC <- qr(A_L)$rank
    if (linC==0)
        QL = NA
    else {
        if (linC==1)
            QL = matrix( A_L[,1]/sqrt(sum(A_L[,1]^2)) )
        else
            QL = svd(A_L)$u[ , 1:linC]

        if (dimC==linC)
            return( list( dimC=dimC, linC=linC, QL=QL, QC=NA, A_reduced=cbind(QL,-QL)) )

        A <- A[ , !is_in_L ]
        A <- A - QL %*% (t(QL) %*% A)
    }
    if (dimC-linC == d)
        QC <- diag(1,d)
    else if (dimC-linC == 1)
        QC <- matrix( svd(A)$u[ , 1] )
    else
        QC <- svd(A)$u[ , 1:(dimC-linC)]

    A <- t(QC) %*% A

    # remove redundancies
    d <- dim(A)[1]
    nn <- dim(A)[2]
    is_relevant <- vector("logical", nn)
    if (solver=="nnls") {
        for (j in 1:nn) {
            is_relevant[j] <- isTRUE(nnls::nnls(A[,-j], A[,j])$deviance > d*tol)
        }
    } else {
        for (j in 1:nn) {
            mos_inp <- conivol:::.create_mosek_input_polyh_prim(A[ ,-j], A[ ,j])
            is_relevant[j] <- isTRUE( sum(
                (Rmosek::mosek(mos_inp,opts)$sol$itr$xx[(nn+2):(nn+d+1)]-A[ ,j])^2 ) > d*tol )
        }
    }
    return( list( dimC=dimC, linC=linC, QL=QL, QC=QC, A_reduced=A[,is_relevant] ) )
}



#' Sample from bivariate chi-bar-squared distribution of a polyhedral cone given by generators
#'
#' \code{polyh_rbichibarsq_gen} generates an \code{n} by \code{2} matrix
#' such that the rows form iid samples from the bivariate chi-bar-squared
#' distribution of the polyhedral cone given by generators, that is, in the
#' form \code{{Ax|x>=0}}. If \code{reduce==TRUE}, which is the default, then a
#' reduced form of the cone will be computed and the bivariate chi-bar-squared
#' distribution will correspond to the reduced form, and the output will contain
#' further elements (in form of a list), see below.
#'
#' @param n number of samples
#' @param A matrix
#' @param solver either "nnls" or "mosek"
#' @param reduce logical; if \code{TRUE}, the cone defined by \code{A} will be
#'               decomposed orthogonally w.r.t. its lineality space
#' @param tol tolerance used in the reduction step
#'            (single precision machine epsilon by default)
#'
#' @return The output of \code{polyh_rbichibarsq_gen(n,A)}, with the default value
#'         \code{reduce==TRUE}, is a list containing the following elements:
#' \itemize{
#'   \item \code{dimC}: the dimension of the linear span of \code{C},
#'   \item \code{linC}: the lineality of the cone \code{C},
#'   \item \code{QL}: an orthogonal basis of the lineality space of \code{C},
#'                    set to \code{NA} if lineality space is zero-dimensional,
#'   \item \code{QC}: an orthogonal basis of the projection of \code{C} onto
#'                    the orthogonal complement of the lineality space of \code{C},
#'                    set to \code{NA} if \code{C} is a linear space,
#'   \item \code{A_reduced}: a matrix defining the reduced cone,
#'   \item \code{samples}: an \code{n} by \code{2} matrix whose rows form
#'         iid samples from the bivariate chi-bar-squared distribution with
#'         weights given by the intrinsic volumes of the reduced cone \code{{A_reduced x|x>=0}};
#'         set to \code{NA} if \code{C} is a linear space.
#' }
#'         If \code{reduce==FALSE} then the output is only an
#'         \code{n} by \code{2} matrix such that its rows form
#'         iid samples from the bivariate chi-bar-squared distribution with
#'         weights given by the intrinsic volumes of the cone \code{{Ax|x>=0}}.
#'
#' @note See \href{../doc/conic-intrinsic-volumes.html#sampling_polyh}{this vignette}
#'       for further info.
#'
#' @section See also:
#' \code{\link[conivol]{polyh_rbichibarsq_ineq}}, \code{\link[conivol]{rbichibarsq}},
#' \code{\link[conivol]{circ_rbichibarsq}}
#'
#' Package: \code{\link[conivol]{conivol}}
#'
#' @examples
#' set.seed(1234)
#' out <- polyh_rbichibarsq_gen(20, matrix(1:12,4,3))
#' print(out)
#'
#' set.seed(1234)
#' sampmos <- polyh_rbichibarsq_gen(20, out$A_reduced, solver="mosek", reduce=FALSE)
#' print(sampmos)
#'
#' sum( (out$samples - sampmos)^2 )
#'
#' @export
#'
polyh_rbichibarsq_gen <- function(n, A, solver="nnls", reduce=TRUE, tol=1e-7) {
    if (solver=="nnls") {
        if (!requireNamespace("nnls", quietly = TRUE))
            stop("\n Could not find package 'nnls'.")
    } else if (solver=="mosek") {
        if (!requireNamespace("Rmosek", quietly = TRUE))
            stop("\n Could not find package 'Rmosek'.")
        if (!requireNamespace("Matrix", quietly = TRUE))
            stop("\n Could not find package 'Matrix'.")
        opts <- list(verbose=0)
    } else
        stop("\n Parameter solver must have value either 'nnls' or 'mosek'.")
    if (!is.matrix(A))
        stop("\n Input is not a matrix.")

    if (reduce) {
        red <- polyh_reduce(A, solver=solver, tol=tol)
        dimC <- red$dimC
        if (dimC==0)
            return( list( dimC=0, linC=0, QL=NA, QC=NA, A_reduced=0, samples=NA) )
        linC <- red$linC
        QL   <- red$QL
        if (dimC==linC)
            return( list( dimC=dimC, linC=linC, QL=QL, QC=NA, A_reduced=cbind(QL,-QL), samples=NA) )
        QC   <- red$QC
        A    <- red$A_reduced
    }

    d <- dim(A)[1]
    out <- matrix(0,n,2)
    if (solver=="nnls") {
        for (i in 1:n) {
            y <- rnorm(d)
            q <- nnls::nnls(A,y)$deviance
            out[i, ] <- c(sum(y^2)-q,q)
        }
    } else {
        nn <- dim(A)[2]
        mos_inp <- conivol:::.create_mosek_input_polyh_prim(A,rep(0,d))
        for (i in 1:n) {
            y <- rnorm(d)
            mos_inp <- conivol:::.update_mosek_input_polyh_prim(mos_inp,y)
            p <- sum( Rmosek::mosek(mos_inp,opts)$sol$itr$xx[(nn+3):(nn+d+2)]^2 )
            out[i, ] <- c(p,sum(y^2)-p)
        }
    }

    if (!reduce)
        return(out)
    else
        return( list( dimC=dimC, linC=linC, QL=QL, QC=QC, A_reduced=A, samples=out) )
}




#' Sample from bivariate chi-bar-squared distribution of a polyhedral cone given by inequalities
#'
#' \code{polyh_rbichibarsq_ineq} generates an \code{n} by \code{2} matrix
#' such that the rows form iid samples from the bivariate chi-bar-squared
#' distribution of the polyhedral cone given by inequalities, that is, in the
#' form \code{{y|A^Ty<=0}}. If \code{reduce==TRUE}, which is the default, then a
#' reduced form of the cone will be computed and the bivariate chi-bar-squared
#' distribution will correspond to the reduced form, and the output will contain
#' further elements (in form of a list), see below.
#'
#' @param n number of samples
#' @param A matrix
#' @param solver either "nnls" or "mosek"
#' @param reduce logical; if \code{TRUE}, the cone defined by \code{A} will be
#'               decomposed orthogonally w.r.t. its lineality space
#' @param tol tolerance used in the reduction step
#'            (single precision machine epsilon by default)
#'
#' @return The output of \code{polyh_rbichibarsq_ineq(n,A)}, with the default value
#'         \code{reduce==TRUE}, is a list containing the following elements:
#' \itemize{
#'   \item \code{dimC}: the dimension of the linear span of \code{C},
#'   \item \code{linC}: the lineality of the cone \code{C},
#'   \item \code{QL}: an orthogonal basis of the orthogonal complement of the
#'                    linear span of \code{C},
#'                    set to \code{NA} if \code{dim(C)==d},
#'   \item \code{QC}: an orthogonal basis of the projection of \code{C} onto
#'                    the orthogonal complement of the lineality space of \code{C},
#'                    set to \code{NA} if \code{C} is a linear space,
#'   \item \code{A_reduced}: a matrix defining the reduced cone,
#'   \item \code{samples}: an \code{n} by \code{2} matrix whose rows form
#'         iid samples from the bivariate chi-bar-squared distribution with
#'         weights given by the intrinsic volumes of the reduced cone \code{{y|A_reduced^T y<=0}}.
#' }
#'         If \code{reduce==FALSE} then the output is only an
#'         \code{n} by \code{2} matrix such that its rows form
#'         iid samples from the bivariate chi-bar-squared distribution with
#'         weights given by the intrinsic volumes of the cone \code{{y|A^Ty<=0}},
#'         set to \code{NA} if \code{C} is a linear space.
#'
#' @section See also:
#' \code{\link[conivol]{polyh_rbichibarsq_gen}}, \code{\link[conivol]{rbichibarsq}},
#' \code{\link[conivol]{circ_rbichibarsq}}
#'
#' Package: \code{\link[conivol]{conivol}}
#'
#' @note See \href{../doc/conic-intrinsic-volumes.html#sampling_polyh}{this vignette}
#'       for further info.
#'
#' @examples
#' set.seed(1234)
#' out <- polyh_rbichibarsq_ineq(20, matrix(1:12,4,3))
#' print(out)
#'
#' set.seed(1234)
#' sampmos <- polyh_rbichibarsq_ineq(20, out$A_reduced, solver="mosek", reduce=FALSE)
#' print(sampmos)
#'
#' sum( (out$samples - sampmos)^2 )
#'
#' @export
#'
polyh_rbichibarsq_ineq <- function(n, A, solver="nnls", reduce=TRUE, tol=1e-7) {
    if (reduce) {
        d <- dim(A)[1]
        out <- polyh_rbichibarsq_gen(n, A, solver=solver, reduce=TRUE, tol=tol)
        dimCpol <- out$dimC
        linCpol <- out$linC

        out$dimC <- d-linCpol
        out$linC <- d-dimCpol
        out$samples <- out$samples[ , c(2,1) ]
        return(out)
    } else
        return( polyh_rbichibarsq_gen(n, A, solver=solver, reduce=FALSE)[ , c(2,1) ] )
}



#' Sample from intrinsic volumes distribution of a polyhedral cone given by generators
#'
#' \code{polyh_rivols_gen} generates a vector of iid samples from the intrinsic
#' volumes distribution, that is, the distribution on \code{{0,1,...,d}} with the
#' probability for \code{k} given by \code{v_k(C)}, where \code{C} is the
#' polyhedral cone by the generator matrix \code{A}, that is, \code{C={Ax|x>=0}}.
#' If \code{reduce==TRUE}, which is the default, then a
#' reduced form of the cone will be computed and returned, see below;
#' however, the intrinsic volumes distribution will be that of the original
#' (non-reduced) cone.
#'
#' @param n number of samples
#' @param A matrix
#' @param solver either "nnls" or "mosek"
#' @param reduce logical; if \code{TRUE}, the cone defined by \code{A} will be
#'               decomposed orthogonally w.r.t. its lineality space
#' @param tol tolerance used in the reduction step
#'            (single precision machine epsilon by default)
#'
#' @return The output of \code{polyh_rivols_gen(n,A)}, with the default value
#'         \code{reduce==TRUE}, is a list containing the following elements:
#' \itemize{
#'   \item \code{dimC}: the dimension of the linear span of \code{C},
#'   \item \code{linC}: the lineality of the cone \code{C},
#'   \item \code{QL}: an orthogonal basis of the orthogonal complement of the
#'                    linear span of \code{C},
#'                    set to \code{NA} if \code{dim(C)==d},
#'   \item \code{QC}: an orthogonal basis of the projection of \code{C} onto
#'                    the orthogonal complement of the lineality space of \code{C},
#'                    set to \code{NA} if \code{C} is a linear space,
#'   \item \code{A_reduced}: a matrix defining the reduced cone,
#'   \item \code{samples}: an \code{n}-element vector of integers in \code{linC,...,dimC}
#'         representing iid samples from the distribution on \code{{0,1,...,d}} with the
#'         probability for \code{k} given by \code{v_k(C)}, where \code{C={Ax|x>=0}},
#'   \item \code{multsamp}: a \code{(d+1)}-element vector of integers
#'         in \code{0,...,n} that sum up to \code{n} representing the frequency
#'         table of the above categorical samples.
#' }
#'         If \code{reduce==FALSE} then the output is a list containing only
#'         the vectors \code{samples} and \code{multsamp}.
#'
#' @section See also:
#' \code{\link[conivol]{polyh_rivols_ineq}}
#'
#' Package: \code{\link[conivol]{conivol}}
#'
#' @note See \href{../doc/conic-intrinsic-volumes.html#sampling_polyh}{this vignette}
#'       for further info.
#'
#' @examples
#' set.seed(1234)
#' out <- polyh_rivols_gen(20, matrix(1:12,4,3))
#' print(out)
#'
#' set.seed(1234)
#' out$linC + polyh_rivols_gen(20, out$A_reduced, solver="mosek", reduce=FALSE)$samples
#'
#' @export
#'
polyh_rivols_gen <- function(n, A, solver="nnls", reduce=TRUE, tol=1e-7) {
    if (solver=="nnls") {
        if (!requireNamespace("nnls", quietly = TRUE))
            stop("\n Could not find package 'nnls'.")
    } else if (solver=="mosek") {
        if (!requireNamespace("Rmosek", quietly = TRUE))
            stop("\n Could not find package 'Rmosek'.")
        if (!requireNamespace("Matrix", quietly = TRUE))
            stop("\n Could not find package 'Matrix'.")
        opts <- list(verbose=0)
    } else
        stop("\n Parameter solver must have value either 'nnls' or 'mosek'.")
    if (!is.matrix(A))
        stop("\n Input is not a matrix.")

    if (reduce) {
        red <- polyh_reduce(A, solver=solver, tol=tol)
        dimC <- red$dimC
        if (dimC==0)
            return( list( dimC=0, linC=0, QL=NA, QC=NA, A_reduced=0, samples=NA) )
        linC <- red$linC
        QL   <- red$QL
        if (dimC==linC)
            return( list( dimC=dimC, linC=linC, QL=QL, QC=NA, A_reduced=cbind(QL,-QL), samples=NA) )
        QC   <- red$QC
        A_tmp    <- red$A_reduced
    } else
        A_tmp <- A

    d_tmp <- dim(A_tmp)[1]
    samples <- vector("integer",n)
    if (solver=="nnls") {
        for (i in 1:n) {
            y <- rnorm(d_tmp)
            # compute rank of submatrix corresponding to nonzero components in expression of projection
            samples[i] <- qr(A_tmp[ ,which(nnls::nnls(A_tmp,y)$x>tol)])$rank
        }
    } else {
        nn <- dim(A_tmp)[2]
        mos_inp <- conivol:::.create_mosek_input_polyh_prim(A_tmp,rep(0,d_tmp))
        for (i in 1:n) {
            y <- rnorm(d_tmp)
            mos_inp <- conivol:::.update_mosek_input_polyh_prim(mos_inp,y)
            mos_out <- Rmosek::mosek(mos_inp,opts)
            samples[i] <- qr(A_tmp[ ,which(mos_out$sol$itr$xx[1:nn]>tol)])$rank
        }
    }

    d <- dim(A)[1]
    multsamp <- rep(0,d+1)
    if (!reduce) {
        multsamp <- tabulate( 1+samples, d+1 )
        return( list (samples=samples, multsamp=multsamp ) )
    } else {
        multsamp <- tabulate( 1+samples+linC, d+1 )
        return( list( dimC=dimC, linC=linC, QL=QL, QC=QC, A_reduced=A_tmp,
                      samples=samples+linC, multsamp=multsamp) )
    }
}


#' Sample from intrinsic volumes distribution of a polyhedral cone given by inequalities
#'
#' \code{polyh_rivols_ineq} generates a vector of iid samples from the intrinsic
#' volumes distribution, that is, the distribution on \code{{0,1,...,d}} with the
#' probability for \code{k} given by \code{v_k(C)}, where \code{C} is the
#' polyhedral cone by the inequalities matrix \code{A}, that is, \code{C={y|A^Ty<=0}}.
#' If \code{reduce==TRUE}, which is the default, then a
#' reduced form of the cone will be computed and returned, see below;
#' however, the intrinsic volumes distribution will be that of the original
#' (non-reduced) cone.
#'
#' @param n number of samples
#' @param A matrix
#' @param solver either "nnls" or "mosek"
#' @param reduce logical; if \code{TRUE}, the cone defined by \code{A} will be
#'               decomposed orthogonally w.r.t. its lineality space
#' @param tol tolerance used in the reduction step
#'            (single precision machine epsilon by default)
#'
#' @return The output of \code{polyh_rivols_ineq(n,A)}, with the default value
#'         \code{reduce==TRUE}, is a list containing the following elements:
#' \itemize{
#'   \item \code{dimC}: the dimension of the linear span of \code{C},
#'   \item \code{linC}: the lineality of the cone \code{C},
#'   \item \code{QL}: an orthogonal basis of the orthogonal complement of the
#'                    linear span of \code{C},
#'                    set to \code{NA} if \code{dim(C)==d},
#'   \item \code{QC}: an orthogonal basis of the projection of \code{C} onto
#'                    the orthogonal complement of the lineality space of \code{C},
#'                    set to \code{NA} if \code{C} is a linear space,
#'   \item \code{A_reduced}: a matrix defining the reduced cone,
#'   \item \code{samples}: an \code{n}-element vector of integers in \code{linC,...,dimC}
#'         representing iid samples from the distribution on \code{{0,1,...,d}} with the
#'         probability for \code{k} given by \code{v_k(C)}, where \code{C={y|A^Ty<=0}},
#'   \item \code{multsamp}: a \code{(d+1)}-element vector of integers
#'         in \code{0,...,n} that sum up to \code{n} representing the frequency
#'         table of the above categorical samples.
#' }
#'         If \code{reduce==FALSE} then the output is a list containing only
#'         the vectors \code{samples} and \code{multsamp}.
#'
#' @section See also:
#' \code{\link[conivol]{polyh_rivols_gen}}
#'
#' Package: \code{\link[conivol]{conivol}}
#'
#' @note See \href{../doc/conic-intrinsic-volumes.html#sampling_polyh}{this vignette}
#'       for further info.
#'
#' @examples
#' set.seed(1234)
#' out <- polyh_rivols_ineq(20, matrix(1:12,4,3))
#' print(out)
#'
#' set.seed(1234)
#' out$linC + polyh_rivols_ineq(20, out$A_reduced, solver="mosek", reduce=FALSE)$samples
#'
#' @export
#'
polyh_rivols_ineq <- function(n, A, solver="nnls", reduce=TRUE, tol=1e-7) {
    d <- dim(A)[1]
    if (reduce) {
        out <- polyh_rivols_gen(n, A, solver=solver, reduce=TRUE, tol=tol)
        dimCpol <- out$dimC
        linCpol <- out$linC

        out$dimC <- d-linCpol
        out$linC <- d-dimCpol
        out$samples <- d-out$samples
        out$multsamp <- rev(out$multsamp)
        return(out)
    } else {
        out <- polyh_rivols_gen(n, A, solver=solver, reduce=FALSE, tol=tol)
        out$samples <- d-out$samples
        out$multsamp <- rev(out$multsamp)
        return(out)
    }
}


#' Bayesian posterior for samples of intrinsic volumes distribution
#'
#' \code{polyh_bayes} generates functions for computing quantiles of marginals
#' of the posterior distribution and for sampling from the posterior distribution,
#' given direct (multinomial) samples of the intrinsic volumes distribution.
#'
#' @param multsamp vector of integers representing a sample from the
#'                 multinomial intrinsic volumes distribution of a convex cone
#' @param dimC the dimension of the cone
#' @param linC the lineality of the cone
#' @param prior either "noninformative" (default) or "informative"
#' @param v_prior a prior estimate of the vector of intrinsic volumes (NA by default)
#'
#' @return The output of \code{polyh_bayes} is a list containing the following elements:
#' \itemize{
#'   \item \code{post_marg_quant}: a function that computes the quantiles of the
#'                    marginals of the posterior distribution;
#'                    \code{marg_quant(i,alpha)} returns the value \code{x}
#'                    such that \code{Prob(v_i<x)=alpha};
#'                    the index \code{i} as well as the probability \code{alpha}
#'                    may be vectors,
#'   \item \code{post_samp}: a function that returns samples of the posterior distribution;
#'                    \code{post_samp(n)} returns an \code{n}-by-\code{(dimC+1)}
#'                    matrix whose rows form a set of \code{n} independent samples
#'                    of the posterior distribution,
#'   \item \code{Dir}: a list containing the weights of the Dirichlet distributions
#'                    that make up both prior and posterior distributions.
#' }
#'
#' @section See also:
#' \code{\link[conivol]{polyh_rivols_gen}}, \code{\link[conivol]{polyh_rivols_ineq}},
#' \code{\link[conivol]{polyh_stan}}
#'
#' Package: \code{\link[conivol]{conivol}}
#'
#' @note See \href{../doc/bayesian.html#sampl_latent}{this vignette}
#'       for further info.
#'
#' @examples
#' # set parameters of cones
#' D <- c(5,7)
#' cone_types <- c("BC","BCp")
#' d <- sum(D)
#'
#' # collect matrix representation and true intrinsic volumes
#' v <- weyl_ivols(D, cone_types, product = TRUE)
#' A <- weyl_matrix(D, cone_types, product = TRUE)
#' true_data <- list( ivols=v, A=A )
#' print(true_data)
#'
#' # collect sample data from intrinsic volumes distribution
#' n <- 10^4
#' set.seed(1234)
#' out <- polyh_rivols_ineq(n,A)
#' str(out)
#'
#' # evaluate posterior distribution
#' bayes_est <- polyh_bayes( out$multsamp, out$dimC, out$linC )
#' str(bayes_est)
#'
#' # compare posterior median with true values
#' v_est_med <- bayes_est$post_marg_quant(0:sum(D),0.5)
#' v_est_med / v
#' sum( (v_est_med-v)^2 )
#'
#' # display boxplot of posterior distribution, overlayed with true values
#' data <- as.data.frame( bayes_est$post_samp(1e4) )
#' colnames(data) <- paste0(rep("V",d+1),as.character(0:d))
#' boxplot( value~key, tidyr::gather( data, factor_key=TRUE ) )
#' lines(1+0:d, v, col="red")
#' lines(1+0:d, v_est_med, col="blue")
#'
#' # display boxplot of posterior distribution of logs, overlayed with true values
#' data <- as.data.frame( log(bayes_est$post_samp(1e4)) )
#' colnames(data) <- paste0(rep("logV",d+1),as.character(0:d))
#' boxplot( value~key, tidyr::gather( data, factor_key=TRUE ) )
#' lines(1+0:d, log(v), col="red")
#' lines(1+0:d, log(v_est_med), col="blue")
#'
#' @export
#'
polyh_bayes <- function(multsamp, dimC, linC, prior="noninformative", v_prior=NA) {
    if ( !(prior %in% c("noninformative", "informative")) )
        stop("\n Parameter prior must be \"noninformative\" or \"informative\".")
    if ( linC>dimC )
        stop("\n Lineality linC must be less than dimension dimC.")
    if ( linC==dimC )
        stop("\n Lineality and dimension (linC==dimC) indicate that cone is linear subspace.")

    d <- dimC-linC
    I_even <- 2*(0:floor(d/2)) + 1          # final +1 is because of R indices start at 1
    I_odd  <- 1+2*(0:floor((d-1)/2)) + 1    # final +1 is because of R indices start at 1

    if (is.na(v_prior)) {
        v_prior_adj <- rep(0,d+1)
        v_prior_adj[I_even] <- 1/ceiling((d+1)/2) / 2
        v_prior_adj[I_odd]  <- 1/floor((d+1)/2) / 2
    } else {
        v_prior_adj <- v_prior[linC:dimC]
    }

    Dir_prior_even <- 2*v_prior_adj[I_even]
    Dir_prior_odd  <- 2*v_prior_adj[I_odd]

    if (prior=="informative"){
        Dir_prior_even <- 1+Dir_prior_even
        Dir_prior_odd  <- 1+Dir_prior_odd
    }
    update <- multsamp[linC:(dimC-linC)+1]

    Dir_post_even <- Dir_prior_even + update[I_even]
    Dir_post_odd  <- Dir_prior_odd  + update[I_odd]

    marg_quant <- function(i,alpha) {
        m <- max(length(i),length(alpha))
        i <- rep_len(i,m)
        alpha <- rep_len(alpha,m)
        X <- rep(0,m)
        for (j in 1:m) {
            if ( i[j] < linC || i[j] > dimC)
                X[j] <- 0
            else {
                if ((i[j]-linC)%%2==0)
                    X[j] <- qbeta(alpha[j], Dir_post_even[(i[j]-linC)/2+1],
                                  sum(Dir_post_even[-((i[j]-linC)/2+1)])) / 2
                else
                    X[j] <- qbeta(alpha[j], Dir_post_odd[(i[j]-linC-1)/2+1],
                                  sum(Dir_post_odd[-((i[j]-linC-1)/2+1)])) / 2
            }
        }
        return(X)
    }

    post_samp <- function(n) {
        samp <- matrix(0,n,dimC+1)
        for (j in I_even)
            samp[ ,j+linC] <- rgamma(n,Dir_post_even[(j-1)/2+1])
        for (j in I_odd)
            samp[ ,j+linC] <- rgamma(n,Dir_post_odd[(j-2)/2+1])
        samp[ ,I_even+linC] <- samp[ ,I_even+linC] / rowSums(samp[ ,I_even+linC]) / 2
        samp[ ,I_odd+linC]  <- samp[ ,I_odd+linC]  / rowSums(samp[ ,I_odd+linC])  / 2
        return(samp)
    }

    out <- list()
    out$post_marg_quant <- marg_quant
    out$post_samp  <- post_samp
    out$Dir <- list(prior=list(even=Dir_prior_even, odd=Dir_prior_odd),
                    post =list(even=Dir_post_even,  odd=Dir_post_odd))
    return(out)
}





#' Stan model creation for Bayesian posterior given direct samples, enforcing log-concavity
#'
#' \code{polyh_stan} generates inputs for Stan (data list and model string or external file)
#' for sampling from the posterior distribution,
#' given direct (multinomial) samples of the intrinsic volumes distribution.
#' The prior distribution is taken on the log-concavity parameters
#' (second iterated differences of the logarithms of the intrinsic volumes),
#' which enforces log-concavity of the intrinsic volumes.
#'
#' @param multsamp vector of integers representing a sample from the
#'                 multinomial intrinsic volumes distribution of a convex cone
#' @param dimC the dimension of the cone
#' @param linC the lineality of the cone
#' @param prior either "noninformative" (default) or "informative"
#' @param v_prior a prior estimate of the vector of intrinsic volumes (NA by default)
#' @param filename filename for output (NA by default, in which case the return is a string)
#' @param overwrite logical; determines whether the output should overwrite an existing file
#'
#' @return If \code{filename==NA} then the output of \code{polyh_stan} is a list containing the following elements:
#' \itemize{
#'   \item \code{model}: a string that forms the description of the Stan model,
#'   \item \code{data}: a data list containing the prepared data to be used
#'                    for defining a Stan model object,
#'                    ## rest still has to be adapted ##
#'   \item \code{variable.names}: the single string "V" to be used as additional
#'                    parameter when creating samples from the JAGS object to
#'                    indicate that only this vector should be tracked.
#' }
#' If \code{filename!=NA} then the model string will be written to the file with
#' the specified name and the output will only contain the elements \code{data}
#' and \code{variable.names}.
#'
#' @section See also:
#' \code{\link[conivol]{polyh_rivols_gen}}, \code{\link[conivol]{polyh_rivols_ineq}},
#' \code{\link[conivol]{polyh_bayes}}
#'
#' Package: \code{\link[conivol]{conivol}}
#'
#' @note See \href{../doc/bayesian.html#sampl_latent}{this vignette}
#'       for further info.
#'
#' @examples
#' # set parameters of cones
#' D <- c(5,7)
#' cone_types <- c("BC","BCp")
#' d <- sum(D)
#'
#' # collect matrix representation and true intrinsic volumes
#' v <- weyl_ivols(D, cone_types, product = TRUE)
#' A <- weyl_matrix(D, cone_types, product = TRUE)
#' true_data <- list( ivols=v, A=A )
#' print(true_data)
#'
#' # collect sample data from intrinsic volumes distribution
#' n <- 10^4
#' set.seed(1234)
#' out <- polyh_rivols_ineq(n,A)
#' str(out)
#'
#' # define stan model
#' filename <- "ex_stan_model.stan"
#' staninp <- polyh_stan(out$multsamp, out$dimC, out$linC, prior="informative", filename=filename)
#'
#' # run the stan model
#' stanfit <- stan( file = filename, data = staninp$data, chains = 4,
#'                  warmup = 1000, iter = 2000, cores = 2, refresh = 200 )
#' str(extract(stanfit))
#'
#' # remove stan file
#' file.remove(filename)
#'
#' # compare posterior median with true values
#' v_est_med <- apply(extract(stanfit)$V, 2, FUN = median)
#' v_est_med / v
#' sum( (v_est_med-v)^2 )
#'
#' # display boxplot of posterior distribution, overlayed with true values
#' data <- as.data.frame( extract(stanfit)$V )
#' colnames(data) <- paste0(rep("V",d+1),as.character(0:d))
#' boxplot( value~key, tidyr::gather( data, factor_key=TRUE ) )
#' lines(1+0:d, v, col="red")
#' lines(1+0:d, v_est_med, col="blue")
#'
#' # display boxplot of posterior distribution of logs, overlayed with true values
#' data <- as.data.frame( extract(stanfit)$logV_nonz )
#' colnames(data) <- paste0(rep("logV",d+1),as.character(0:d))
#' boxplot( value~key, tidyr::gather( data, factor_key=TRUE ) )
#' lines(1+0:d, log(v), col="red")
#' lines(1+0:d, log(v_est_med), col="blue")
#'
#' @export
#'
polyh_stan <- function(multsamp, dimC, linC, prior="noninformative", v_prior=NA, filename=NA, overwrite=FALSE) {
    if ( !(prior %in% c("noninformative", "informative")) )
        stop("\n Parameter prior must be \"noninformative\" or \"informative\".")
    if ( linC>dimC )
        stop("\n Lineality linC must be less than dimension dimC.")
    if ( linC==dimC )
        stop("\n Lineality and dimension (linC==dimC) indicate that cone is linear subspace.")

    if (is.na(v_prior)) {
        v_nonz_prior <- rep(1,dimC-linC+1)/(dimC-linC+1)
    } else {
        v_nonz_prior <- v_prior[linC:dimC]
    }

    if (prior=="informative"){
        alpha <- v_nonz_prior / 2
        beta  <- rep(1/2,dimC-linC+1)
    } else {
        alpha <- rep(1,dimC-linC+1)
        beta  <- 1 / v_nonz_prior
    }

    T <- matrix( rep(0,(dimC-linC+1)^2), dimC-linC+1, dimC-linC+1 )
    for (i in 1:(dimC-linC-1)) {
        T[i,i] = 1
        T[i,i+1] = -2
        T[i,i+2] = 1
    }
    if ((dimC-linC)%%2==0) {
        for (i in 0:((dimC-linC)/2-1)) {
            T[dimC-linC,  2*i+1] = 1
            T[dimC-linC+1,2*i+2] = 1
        }
        T[dimC-linC,  dimC-linC+1] = 1
        T[dimC-linC+1,1] = 1
    } else {
        for (i in 0:((dimC-linC-1)/2)) {
            T[dimC-linC,  2*i+1] = 1
            T[dimC-linC+1,2*i+2] = 1
        }
    }

    data_list <- list(
        d           = length(multsamp)-1 ,
        dimC        = dimC ,
        linC        = linC ,
        multsamp    = multsamp ,
        alpha       = alpha ,
        beta        = beta ,
        T           = T
    )
        model_string <-
"// Log-concavity enforcing model for estimating conic intrinsic volumes
// given a multinomial sample from the intrinsic volumes distribution.
// Model is created with the R-method 'polyh_stan' from the 'conivol' package.
// See 'https://github.com/damelunx/conivol' for more information.

data {
    int<lower=1> d ;                                        // ambient dimension
    int<lower=1, upper=d> dimC ;                            // dimension of linear span of C
    int<lower=0, upper=dimC-1> linC ;                       // lineality
    int<lower=0> multsamp[d+1] ;                            // sample of multinomial distribution (even and odd)
    vector<lower=0>[dimC-linC+1] alpha ;                    // prior values for hyperparameters alpha
    vector<lower=0>[dimC-linC+1] beta ;                     // prior values for hyperparameters beta
    matrix[dimC-linC+1,dimC-linC+1] T ;                     // transformation matrix for u ~> t
}

transformed data {
    int<lower=0> samp_nonz[dimC-linC+1] ;
    samp_nonz = multsamp[(linC+1):(dimC+1)] ;
}

parameters {
    vector<lower=0>[dimC-linC+1] t ;
}

transformed parameters {
    vector<lower=0>[dimC-linC+1] v_nonz ;
    v_nonz = exp(- T \\ t) ;
}

model {
    for (k in 0:(dimC-linC)) {
        t[k+1] ~ gamma(alpha[k+1], beta[k+1]) ;
    }
    samp_nonz ~ multinomial(v_nonz/sum(v_nonz)) ;
}

generated quantities {
    vector[d+1] V ;
    vector[dimC-linC+1] logV_nonz ;
    if (linC>0) V[1:linC] = rep_vector(0,linC) ;
    if (dimC<d) V[(dimC+2):(d+1)] = rep_vector(0,d-dimC) ;
    V[(linC+1):(dimC+1)] = v_nonz / sum(v_nonz) ;
    logV_nonz = log(v_nonz) - log(sum(v_nonz)) ;
}"

    out                <- list()
    out$data           <- data_list
    out$variable.names <- c("V","logV_nonz","t")

    if ( !is.na(filename) && file.exists(filename) && overwrite==FALSE )
        stop("\n File with given filename exists and overwrite==FALSE.")
    else if ( !is.na(filename) ) {
        if (file.exists(filename))
            file.remove(filename)
        file_conn<-file(filename)
        writeLines(model_string, file_conn)
        close(file_conn)
    } else {
        out$model <- model_string
    }
    return(out)

}
