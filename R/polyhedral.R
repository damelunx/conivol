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


#' Sample from bivariate chi-bar-squared distribution of polyhedral cones
#'
#' \code{rbichibarsq_polyh} generates an \code{n} by \code{2} matrix
#' such that the rows form iid samples from the bivariate chi-bar-squared
#' distribution of the polyhedral cone \code{{Ax|x>=0}}. If input parameter
#' \code{reduce==TRUE} then output will contain further elements (in form of a list),
#' see below.
#'
#' @param n number of samples
#' @param A matrix
#' @param reduce logical; if \code{TRUE}, \code{A} will be replaced by a smaller
#'               matrix if possible
#'
#' @return If \code{reduce==TRUE} then the output of \code{rbichibarsq_polyh(n,A)}
#'         is a list containing an \code{n} by \code{2} matrix \code{samples},
#'         a logical \code{replaced}, and, if \code{replaced==TRUE}, another
#'         matrix \code{A_reduced}; if \code{reduce==FALSE} then the output
#'         is only the \code{n} by \code{2} matrix. The rows of the matrix form
#'         iid samples from the bivariate chi-bar-squared distribution with
#'         weights given by the intrinsic volumes of the polyhedral cone \code{{Ax|x>=0}}.
#'
#' @note See \href{link.to.vignette/TBD}{this vignette}
#'       for further info; in particular, info about how
#'       to work with polyhedral cones given in the form \code{{x|Ax>=0}}.
#'
#' @section See also:
#' \code{\link[conivol]{rbichibarsq}}, \code{\link[conivol]{rbichibarsq_circ}},
#'
#' @examples
#' rbichibarsq_polyh(20,matrix(1:12,4,3)
#'
#' @export
#'
rbichibarsq_polyh <- function(n, A, reduce=TRUE) {
    if (!requireNamespace("Rmosek", quietly = TRUE))
        stop("\n Could not find package 'Rmosek'.")
    if (!requireNamespace("Matrix", quietly = TRUE))
        stop("\n Could not find package 'Matrix'.")
    opts <- list(verbose=0)
    replaced <- FALSE
    if (reduce) {
        # test if cone lies in lower-dimensional space
        QR <- qr(A)
        d <- QR$rank
        if ( d<dim(A)[1] ) {
            A <- qr.R(QR)[1:d, ]
            replaced <- TRUE
        }
        # test whether columns of A lie in the cone spanned by the other columns
        m <- dim(A)[1]
        nn <- dim(A)[2]
        is_superfluous <- rep(0,nn)
        for (j in 1:nn) {
            mos_inp <- .create_mosek_input_polyh_prim(A[ ,-j],A[ ,j])
            is_superfluous[j] <- sum( Rmosek::mosek(mos_inp,opts)$sol$itr$xx[(nn+2):(nn+m+1)]^2 ) == sum(A[ ,j]^2)
        }
        if (any(is_superfluous>0)) {
            A <- A[ ,-is_superfluous]
            replaced = TRUE
        }
    }

    m <- dim(A)[1]
    nn <- dim(A)[2]
    mos_inp <- .create_mosek_input_polyh_prim(A,rep(0,m))
    out <- matrix(0,n,2)
    for (i in 1:n) {
        y <- rnorm(m)
        mos_inp <- .update_mosek_input_polyh_prim(mos_inp,y)
        p <- sum( Rmosek::mosek(mos_inp,opts)$sol$itr$xx[(nn+3):(nn+m+2)]^2 )
        out[i, ] <- c(p,sum(y^2)-p)
    }
    if (!reduce)
        return(out)
    else {
        if (replaced)
            return(list(replaced=TRUE,  A_reduced=A, samples=out))
        else
            return(list(replaced=FALSE, samples=out))
    }
}







