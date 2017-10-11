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
#'   \item \code{dim}: the dimension of the linear span of \code{C},
#'   \item \code{lin}: the lineality of the cone \code{C},
#'   \item \code{QL}: an orthogonal basis of the lineality space of \code{C},
#'                    set to \code{NA} if lineality space is zero-dimensional,
#'   \item \code{QC}: an orthogonal basis of the projection of \code{C} onto
#'                    the orthogonal complement of the lineality space of \code{C},
#'                    set to \code{NA} if \code{C} is a linear space,
#'   \item \code{A_reduced}: a matrix defining the reduced cone.
#' }
#'
#'
#' @note See \href{../doc/conic-intrinsic-volumes.html#sampling_polyh}{this vignette}
#'       for further info.
#'
#' Package: \code{\link[conivol]{conivol}}
#'
#' @examples
#' A <- cbind(diag(1,3),diag(1,3)+matrix(1,ncol=3,nrow=3))
#' A <- A[,sample(ncol(A))]
#' print(A)
#' A_red <- polyh_reduce(A)$A_reduced
#' print(A_red)
#'
#' @export
#'
polyh_reduce <- function(A, solver="nnls", tol=1e-8) {
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
        return( list( dim=0, lin=0, QL=NA, QC=NA, A_reduced=0) )
    # find lineality space L
    nn <- dim(A)[2]
    is_in_L <- vector("logical",nn)
    if (solver=="nnls") {
        for (j in 1:nn) {
            is_in_L[j] <- isTRUE(nnls::nnls(A[,-j],-A[,j])$deviance<d*tol)
        }
    } else {
        for (j in 1:nn) {
            mos_inp <- .create_mosek_input_polyh_prim(A[ ,-j], -A[ ,j])
            # is_in_L[j] <- all.equal( Rmosek::mosek(mos_inp,opts)$sol$itr$xx[(nn+2):(nn+d+1)] , A[ ,j] ) == TRUE
            is_in_L[j] <- isTRUE( sum( Rmosek::mosek(mos_inp,opts)$sol$itr$xx[(nn+2):(nn+d+1)]^2 ) < d*tol )
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
            return( list( dim=dimC, lin=linC, QL=QL, QC=NA, A_reduced=cbind(QL,-QL)) )

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
            mos_inp <- .create_mosek_input_polyh_prim(A[ ,-j], A[ ,j])
            is_relevant[j] <- isTRUE( sum( Rmosek::mosek(mos_inp,opts)$sol$itr$xx[(nn+2):(nn+d+1)]^2 ) > d*tol )
        }
    }
    A <- A[,is_relevant]
    return( list( dim=dimC, lin=linC, QL=QL, QC=QC, A_reduced=A ) )
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
#'   \item \code{dim}: the dimension of the linear span of \code{C},
#'   \item \code{lin}: the lineality of the cone \code{C},
#'   \item \code{QL}: an orthogonal basis of the lineality space of \code{C},
#'                    set to \code{NA} if lineality space is zero-dimensional,
#'   \item \code{QC}: an orthogonal basis of the projection of \code{C} onto
#'                    the orthogonal complement of the lineality space of \code{C},
#'                    set to \code{NA} if \code{C} is a linear space,
#'   \item \code{A_reduced}: a matrix defining the reduced cone,
#'   \item \code{samples}: an \code{n} by \code{2} matrix rows of the matrix form
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
#' polyh_rbichibarsq_gen_nnls(20,matrix(1:12,4,3))
#'
#' @export
#'
polyh_rbichibarsq_gen <- function(n, A, solver="nnls", reduce=TRUE, tol=1e-8) {
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
        dimC <- red$dim
        if (dimC==0)
            return( list( dim=0, lin=0, QL=NA, QC=NA, A_reduced=0, samples=NA) )
        linC <- red$lin
        QL   <- red$QL
        if (dimC==linC)
            return( list( dim=dimC, lin=linC, QL=QL, QC=NA, A_reduced=cbind(QL,-QL), samples=NA) )
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
        mos_inp <- .create_mosek_input_polyh_prim(A,rep(0,d))
        for (i in 1:n) {
            y <- rnorm(d)
            mos_inp <- .update_mosek_input_polyh_prim(mos_inp,y)
            p <- sum( Rmosek::mosek(mos_inp,opts)$sol$itr$xx[(nn+3):(nn+d+2)]^2 )
            out[i, ] <- c(p,sum(y^2)-p)
        }
    }

    if (!reduce)
        return(out)
    else
        return( list( dim=dimC, lin=linC, QL=QL, QC=QC, A_reduced=A, samples=out) )
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
#'   \item \code{dim}: the dimension of the linear span of \code{C},
#'   \item \code{lin}: the lineality of the cone \code{C},
#'   \item \code{QL}: an orthogonal basis of the orthogonal complement of the
#'                    linear span of \code{C},
#'                    set to \code{NA} if \code{dim(C)==d},
#'   \item \code{QC}: an orthogonal basis of the projection of \code{C} onto
#'                    the orthogonal complement of the lineality space of \code{C},
#'                    set to \code{NA} if \code{C} is a linear space,
#'   \item \code{A_reduced}: a matrix defining the reduced cone,
#'   \item \code{samples}: an \code{n} by \code{2} matrix rows of the matrix form
#'         iid samples from the bivariate chi-bar-squared distribution with
#'         weights given by the intrinsic volumes of the reduced cone \code{{y|A_reduced^T y<=0}}.
#' }
#'         If \code{reduce==FALSE} then the output is only an
#'         \code{n} by \code{2} matrix such that its rows form
#'         iid samples from the bivariate chi-bar-squared distribution with
#'         weights given by the intrinsic volumes of the cone \code{{y|A^Ty<=0}},
#'         set to \code{NA} if \code{C} is a linear space.
#'
#' @note See \href{../doc/conic-intrinsic-volumes.html#sampling_polyh}{this vignette}
#'       for further info.
#'
#' @section See also:
#' \code{\link[conivol]{polyh_rbichibarsq_gen}}, \code{\link[conivol]{rbichibarsq}},
#' \code{\link[conivol]{circ_rbichibarsq}}
#'
#' Package: \code{\link[conivol]{conivol}}
#'
#' @examples
#' polyh_rbichibarsq_ineq(20,matrix(1:12,4,3)
#'
#' @export
#'
polyh_rbichibarsq_ineq <- function(n, A, solver="nnls", reduce=TRUE, tol=1e-8) {
    if (reduce) {
        d <- dim(A)[1]
        out <- polyh_rbichibarsq_gen(n, A, solver=solver, reduce=TRUE, tol=tol)
        dimCpol <- out$dim
        linCpol <- out$lin

        out$dim <- d-linCpol
        out$lin <- d-dimCpol
        out$samples <- out$samples[ , c(2,1) ]
        return(out)
    } else
        return( polyh_rbichibarsq_gen(n, A, solver=solver, reduce=FALSE)[ , c(2,1) ] )
}



#' Sample from intrinsic volumes distribution of a polyhedral cone given by generators
#'
#' \code{polyh_samp_ivol_gen} generates a vector of iid samples from the intrinsic
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
#' @return The output of \code{polyh_samp_ivol_gen(n,A)}, with the default value
#'         \code{reduce==TRUE}, is a list containing the following elements:
#' \itemize{
#'   \item \code{dim}: the dimension of the linear span of \code{C},
#'   \item \code{lin}: the lineality of the cone \code{C},
#'   \item \code{QL}: an orthogonal basis of the orthogonal complement of the
#'                    linear span of \code{C},
#'                    set to \code{NA} if \code{dim(C)==d},
#'   \item \code{QC}: an orthogonal basis of the projection of \code{C} onto
#'                    the orthogonal complement of the lineality space of \code{C},
#'                    set to \code{NA} if \code{C} is a linear space,
#'   \item \code{A_reduced}: a matrix defining the reduced cone,
#'   \item \code{samples}: an \code{n}-element vector whose elements form
#'         iid samples from the distribution on \code{{0,1,...,d}} with the
#'         probability for \code{k} given by \code{v_k(C)}, where \code{C={Ax|x>=0}}.
#' }
#'         If \code{reduce==FALSE} then the output is only the above vector of samples.
#'
#' @note See \href{../doc/conic-intrinsic-volumes.html#sampling_polyh}{this vignette}
#'       for further info.
#'
#' @section See also:
#' \code{\link[conivol]{polyh_samp_ivol_ineq}}
#'
#' Package: \code{\link[conivol]{conivol}}
#'
#' @examples
#' polyh_samp_ivol_gen(20,matrix(1:12,4,3)
#'
#' @export
#'
polyh_samp_ivol_gen <- function(n, A, solver="nnls", reduce=TRUE, tol=1e-8) {
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
        dimC <- red$dim
        if (dimC==0)
            return( list( dim=0, lin=0, QL=NA, QC=NA, A_reduced=0, samples=NA) )
        linC <- red$lin
        QL   <- red$QL
        if (dimC==linC)
            return( list( dim=dimC, lin=linC, QL=QL, QC=NA, A_reduced=cbind(QL,-QL), samples=NA) )
        QC   <- red$QC
        A    <- red$A_reduced
    }

    d <- dim(A)[1]
    out <- vector("integer",n)
    if (solver=="nnls") {
        for (i in 1:n) {
            y <- rnorm(d)
            # compute rank of submatrix corresponding to nonzero components in expression of projection
            out[i] <- qr(A[ ,which(nnls::nnls(A,y)$x>tol)])$rank
        }
    } else {
        nn <- dim(A)[2]
        mos_inp <- .create_mosek_input_polyh_prim(A,rep(0,d))
        for (i in 1:n) {
            y <- rnorm(d)
            mos_inp <- .update_mosek_input_polyh_prim(mos_inp,y)
            mos_out <- Rmosek::mosek(mos_inp,opts)
            out[i] <- qr(A[ ,which(mos_out$sol$itr$xx[1:nn]>tol)])$rank
        }
    }

    if (!reduce)
        return(out)
    else
        return( list( dim=dimC, lin=linC, QL=QL, QC=QC, A_reduced=A, samples=out+linC) )
}


#' Sample from intrinsic volumes distribution of a polyhedral cone given by inequalities
#'
#' \code{polyh_samp_ivol_ineq} generates a vector of iid samples from the intrinsic
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
#' @return The output of \code{polyh_samp_ivol_ineq(n,A)}, with the default value
#'         \code{reduce==TRUE}, is a list containing the following elements:
#' \itemize{
#'   \item \code{dim}: the dimension of the linear span of \code{C},
#'   \item \code{lin}: the lineality of the cone \code{C},
#'   \item \code{QL}: an orthogonal basis of the orthogonal complement of the
#'                    linear span of \code{C},
#'                    set to \code{NA} if \code{dim(C)==d},
#'   \item \code{QC}: an orthogonal basis of the projection of \code{C} onto
#'                    the orthogonal complement of the lineality space of \code{C},
#'                    set to \code{NA} if \code{C} is a linear space,
#'   \item \code{A_reduced}: a matrix defining the reduced cone,
#'   \item \code{samples}: an \code{n}-element vector whose elements form
#'         iid samples from the distribution on \code{{0,1,...,d}} with the
#'         probability for \code{k} given by \code{v_k(C)}, where \code{C={y|A^Ty<=0}}.
#' }
#'         If \code{reduce==FALSE} then the output is only the above vector of samples.
#'
#' @note See \href{../doc/conic-intrinsic-volumes.html#sampling_polyh}{this vignette}
#'       for further info.
#'
#' @section See also:
#' \code{\link[conivol]{polyh_samp_ivol_gen}}
#'
#' Package: \code{\link[conivol]{conivol}}
#'
#' @examples
#' polyh_samp_ivol_ineq(20,matrix(1:12,4,3)
#'
#' @export
#'
polyh_samp_ivol_ineq <- function(n, A, solver="nnls", reduce=TRUE, tol=1e-8) {
    d <- dim(A)[1]
    if (reduce) {
        out <- polyh_samp_ivol_gen(n, A, solver=solver, reduce=TRUE, tol=tol)
        dimCpol <- out$dim
        linCpol <- out$lin

        out$dim <- d-linCpol
        out$lin <- d-dimCpol
        out$samples <- d-out$samples
        return(out)
    } else
        return( d-polyh_samp_ivol_gen(n, A, solver=solver, reduce=FALSE, tol=tol) )
}



#' Bayesian posterior for samples of intrinsic volumes distribution
#'
#' \code{polyh_ivol_Bayes} generates functions for computing quantiles of marginals
#' of the posterior distribution and for sampling from the posterior distribution,
#' given samples of the intrinsic volumes distribution.
#'
#' @param samples vector of integers representing independent samples from the
#'                intrinsic volumes distribution of a convex cone
#' @param dim the dimension of the cone
#' @param lin the lineality of the cone
#' @param prior either "noninformative" (default) or "informative"
#' @param v_prior a prior estimate of the vector of intrinsic volumes (NA by default)
#'
#' @return The output of \code{polyh_ivol_Bayes} is a list containing the following elements:
#' \itemize{
#'   \item \code{post_marg_quant}: a function that computes the quantiles of the
#'                    marginals of the posterior distribution;
#'                    \code{marg_quant(i,alpha)} returns the value \code{x}
#'                    such that \code{Prob(v_i<x)=alpha};
#'                    the index \code{i} as well as the probability \code{alpha}
#'                    may be vectors,
#'   \item \code{post_samp}: a function that returns samples of the posterior distribution;
#'                    \code{post_samp(n)} returns an \code{n}-by-\code{(dim+1)}
#'                    matrix whose rows form a set of \code{n} independent samples
#'                    of the posterior distribution,
#'   \item \code{Dir}: a list containing the weights of the Dirichlet distributions
#'                    that make up both prior and posterior distributions.
#' }
#'
#' @note See \href{../doc/bayesian.html#sampl_latent}{this vignette}
#'       for further info.
#'
#' @section See also:
#' \code{\link[conivol]{polyh_samp_ivol_gen}}, \code{\link[conivol]{polyh_samp_ivol_ineq}}
#'
#' Package: \code{\link[conivol]{conivol}}
#'
#' @examples
#' samp <- polyh_samp_ivol_gen(20,matrix(1:12,4,3))
#' polyh_ivol_Bayes(samp$samples, samp$dim, samp$lin)
#'
#' @export
#'
polyh_ivol_Bayes <- function(samples, dim, lin, prior="noninformative", v_prior=NA) {
    # check whether prior is "noninformative" or "informative"
    # check whether lin==dim

    d <- dim-lin
    I_even <- 2*(0:floor(d/2)) + 1          # final +1 is because of R indices start at 1
    I_odd  <- 1+2*(0:floor((d-1)/2)) + 1    # final +1 is because of R indices start at 1

    if (is.na(v_prior)) {
        v_prior_adj <- rep(0,d+1)
        v_prior_adj[I_even] <- 1/ceiling((d+1)/2) / 2
        v_prior_adj[I_odd]  <- 1/floor((d+1)/2) / 2
    } else {
        v_prior_adj <- v_prior[lin:dim]
    }

    Dir_prior_even <- 2*v_prior_adj[I_even]
    Dir_prior_odd  <- 2*v_prior_adj[I_odd]

    if (prior=="informative"){
        Dir_prior_even <- 1+Dir_prior_even
        Dir_prior_odd  <- 1+Dir_prior_odd
    }
    update <- tabulate(1+samples-lin, 1+dim-lin)

    Dir_post_even <- Dir_prior_even + update[I_even]
    Dir_post_odd  <- Dir_prior_odd  + update[I_odd]

    marg_quant <- function(i,alpha) {
        m <- max(length(i),length(alpha))
        i <- rep_len(i,m)
        alpha <- rep_len(alpha,m)
        X <- rep(0,m)
        for (j in 1:m) {
            if ( i[j] < lin || i[j] > dim)
                X[j] <- 0
            else {
                if ((i[j]-lin)%%2==0)
                    X[j] <- qbeta(alpha[j], Dir_post_even[(i[j]-lin)/2+1],
                                  sum(Dir_post_even[-((i[j]-lin)/2+1)])) / 2
                else
                    X[j] <- qbeta(alpha[j], Dir_post_odd[(i[j]-lin-1)/2+1],
                                  sum(Dir_post_odd[-((i[j]-lin-1)/2+1)])) / 2
            }
        }
        return(X)
    }

    post_samp <- function(n) {
        samp <- matrix(0,n,dim+1)
        for (j in I_even)
            samp[ ,j+lin] <- rgamma(n,Dir_post_even[(j-1)/2+1])
        for (j in I_odd)
            samp[ ,j+lin] <- rgamma(n,Dir_post_odd[(j-2)/2+1])
        samp[ ,I_even+lin] <- samp[ ,I_even+lin] / rowSums(samp[ ,I_even+lin]) / 2
        samp[ ,I_odd+lin]  <- samp[ ,I_odd+lin]  / rowSums(samp[ ,I_odd+lin])  / 2
        return(samp)
    }

    out <- list()
    out$post_marg_quant <- marg_quant
    out$post_samp  <- post_samp
    out$Dir <- list(prior=list(even=Dir_prior_even, odd=Dir_prior_odd),
                    post =list(even=Dir_post_even,  odd=Dir_post_odd))
    return(out)
}





