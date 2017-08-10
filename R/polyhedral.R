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


#' Sample from bivariate chi-bar-squared distribution of a polyhedral cone given by generators
#'
#' \code{rbichibarsq_polyh_gen} generates an \code{n} by \code{2} matrix
#' such that the rows form iid samples from the bivariate chi-bar-squared
#' distribution of the polyhedral cone given by generators, that is, in the
#' form \code{{Ax|x>=0}}. If \code{reduce==TRUE}, which is the default, then a
#' reduced form of the cone will be computed and the bivariate chi-bar-squared
#' distribution will correspond to the reduced form, and the output will contain
#' further elements (in form of a list), see below.
#'
#' @param n number of samples
#' @param A matrix
#' @param reduce logical; if \code{TRUE}, the cone defined by \code{A} will be
#'               decomposed orthogonally w.r.t. its lineality space
#'
#' @return The output of \code{rbichibarsq_polyh_gen(n,A)}, with the default value
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
#' \code{\link[conivol]{rbichibarsq_polyh_ineq}}, \code{\link[conivol]{rbichibarsq}},
#' \code{\link[conivol]{rbichibarsq_circ}}
#'
#' Package: \code{\link[conivol]{conivol}}
#'
#' @examples
#' rbichibarsq_polyh_gen(20,matrix(1:12,4,3)
#'
#' @export
#'
rbichibarsq_polyh_gen <- function(n, A, reduce=TRUE) {
    if (!requireNamespace("Rmosek", quietly = TRUE))
        stop("\n Could not find package 'Rmosek'.")
    if (!requireNamespace("Matrix", quietly = TRUE))
        stop("\n Could not find package 'Matrix'.")
    opts <- list(verbose=0)
    if (reduce) {
        d <- dim(A)[1]
        dimC <- qr(A)$rank
        if (dimC==0)
            return( list( dim=0, lin=0, QL=NA, QC=NA, A_reduced=0, samples=NA) )
        else if (dimC==d)
            return( list( dim=d, lin=d, QL=diag(1,d), QC=NA, A_reduced=diag(1,d), samples=NA) )

        # find lineality space L
        nn <- dim(A)[2]
        is_in_L <- vector("logical",nn)
        for (j in 1:nn) {
            mos_inp <- .create_mosek_input_polyh_prim(A[ ,-j], -A[ ,j])
            # is_in_L[j] <- sum( Rmosek::mosek(mos_inp,opts)$sol$itr$xx[(nn+2):(nn+d+1)]^2 ) == sum(A[ ,j]^2)
            is_in_L[j] <- all.equal( Rmosek::mosek(mos_inp,opts)$sol$itr$xx[(nn+2):(nn+d+1)] , A[ ,j] ) == TRUE
        }
        linC <- qr(A[ , is_in_L ])$rank
        if (linC==0)
            QL = NA
        else {
            if (linC==1)
                QL = matrix( svd(A[ , is_in_L ])$u[ , 1] )
            else
                QL = svd(A[ , is_in_L ])$u[ , 1:linC]

            if (dimC==linC)
                return( list( dim=dimC, lin=linC, QL=QL, QC=NA, A_reduced=cbind(QL,-QL), samples=NA) )

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

        # test whether columns of A lie in the cone spanned by the other columns
        e <- dim(A)[1]
        m <- dim(A)[2]
        is_superfluous <- vector("logical",m)
        for (j in 1:m) {
            mos_inp <- .create_mosek_input_polyh_prim(A[ ,-j],A[ ,j])
            # is_superfluous[j] <- sum( Rmosek::mosek(mos_inp,opts)$sol$itr$xx[(m+2):(m+e+1)]^2 ) == sum(A[ ,j]^2)
            is_superfluous[j] <- all.equal( Rmosek::mosek(mos_inp,opts)$sol$itr$xx[(m+2):(m+e+1)] , A[ ,j] ) == TRUE
        }
        A <- A[ ,!is_superfluous]
    }

    d <- dim(A)[1]
    nn <- dim(A)[2]
    mos_inp <- .create_mosek_input_polyh_prim(A,rep(0,d))
    out <- matrix(0,n,2)
    for (i in 1:n) {
        y <- rnorm(d)
        mos_inp <- .update_mosek_input_polyh_prim(mos_inp,y)
        p <- sum( Rmosek::mosek(mos_inp,opts)$sol$itr$xx[(nn+3):(nn+d+2)]^2 )
        out[i, ] <- c(p,sum(y^2)-p)
    }
    if (!reduce)
        return(out)
    else
        return( list( dim=dimC, lin=linC, QL=QL, QC=QC, A_reduced=A, samples=out) )
}



#' Sample from bivariate chi-bar-squared distribution of a polyhedral cone given by inequalities
#'
#' \code{rbichibarsq_polyh_ineq} generates an \code{n} by \code{2} matrix
#' such that the rows form iid samples from the bivariate chi-bar-squared
#' distribution of the polyhedral cone given by inequalities, that is, in the
#' form \code{{y|A^Ty<=0}}. If \code{reduce==TRUE}, which is the default, then a
#' reduced form of the cone will be computed and the bivariate chi-bar-squared
#' distribution will correspond to the reduced form, and the output will contain
#' further elements (in form of a list), see below.
#'
#' @param n number of samples
#' @param A matrix
#' @param reduce logical; if \code{TRUE}, the cone defined by \code{A} will be
#'               decomposed orthogonally w.r.t. its lineality space
#'
#' @return The output of \code{rbichibarsq_polyh_ineq(n,A)}, with the default value
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
#' \code{\link[conivol]{rbichibarsq_polyh_gen}}, \code{\link[conivol]{rbichibarsq}},
#' \code{\link[conivol]{rbichibarsq_circ}}
#'
#' Package: \code{\link[conivol]{conivol}}
#'
#' @examples
#' rbichibarsq_polyh_ineq(20,matrix(1:12,4,3)
#'
#' @export
#'
rbichibarsq_polyh_ineq <- function(n, A, reduce=TRUE) {
    if (reduce) {
        d <- dim(A)[1]
        out <- rbichibarsq_polyh_gen(n, A, TRUE)
        dimCpol <- out$dim
        linCpol <- out$lin

        out$dim <- d-linCpol
        out$lin <- d-dimCpol
        out$samples <- out$samples[ , c(2,1) ]
        return(out)
    } else
        return( rbichibarsq_polyh_gen(n, A, FALSE)[ , c(2,1) ] )
}

