#' The conic intrinsic volumes of (products of) circular cones
#'
#' \code{circ_ivols} computes the conic intrinsic volumes of circular cones,
#' whose dimensions and angles are given in the vectors \code{d} and
#' \code{alpha} (vectors must be of same lengths); if the length of the vectors
#' is one, a single vector is returned; if the length of the vectors is greater
#' than one and \code{product==FALSE}, a list of vectors is returned;
#' if the length of the vectors is greater than one and \code{product==TRUE},
#' a single vector with the intrinsic volumes of the product cone is returned.
#'
#' @param d vector of dimensions; must be same length as \code{alpha}.
#' @param alpha vector of angles; must be same length as \code{d}.
#' @param product logical; if \code{TRUE}, intrinsic volumes of product cone are returned.
#'
#' @return If \code{length(d)==1} or \code{(length(d)>1 & product==TRUE)}
#'         then a single vectors will be returned.
#'         If \code{(length(d)>1 & product==FALSE)} then a list of
#'         vectors will be returned.
#'
#' @section See also:
#' Package: \code{\link[conivol]{conivol}}
#'
#' @examples
#' circ_ivols(5, pi/4)
#' circ_ivols(c(5,5), c(pi/4,pi/8))
#' circ_ivols(c(5,5), c(pi/4,pi/8), product = TRUE)
#'
#' @export
#'
circ_ivols <- function(d, alpha, product = FALSE) {
    if (length(d)!=length(alpha))
        stop("Inputs d and alpha must be of same length.")

    V <- list()
    for (i in 1:length(d)) {
        if (d[i]<=1 || alpha[i]<0 || alpha[i]>pi/2)
            v <- NA
        else if (alpha[i]==0)
            v <- c(0.5, 0.5, rep(0,d[i]-2))
        else if (alpha[i]==pi/2)
            v <- c(rep(0,d[i]-2), 0.5, 0.5)
        else {
            v <- rep(0,d[i]+1)
            # v[1] <- exp( lgamma(d[i]/2)-lgamma((d[i]+1)/2)-lgamma(1/2)+log(d[i]-1)-log(2)+
            #                  log( integrate(function(x){sin(x)^(d[i]-2)},0,pi/2-alpha[i])$value ) )
            v[1] <- 1/2 * pbeta(cos(alpha[i])^2, (d[i]-1)/2, 1/2)

            # v[d[i]+1] <- exp( lgamma(d[i]/2)-lgamma((d[i]+1)/2)-lgamma(1/2)+log(d[i]-1)-log(2)+
            #                       log( integrate(function(x){sin(x)^(d[i]-2)},0,alpha[i])$value ) )
            v[d[i]+1] <- 1/2 * pbeta(sin(alpha[i])^2, (d[i]-1)/2, 1/2)

            k <- 1:(d[i]-1)
            v[2:d[i]] <- exp( lgamma(d[i]/2)-lgamma((k+1)/2)-lgamma((d[i]-k+1)/2)+
                                  (k-1)*log(sin(alpha[i]))+(d[i]-k-1)*log(cos(alpha[i]))-log(2) )
        }
        V[[i]] <- v
    }

    if (length(d)==1)
        return(V[[1]])
    else if (product)
        return(conivol::prod_ivols(V))
    else return(V)
}


#' Computes the semiaxes of an ellipsoidal cone given as the linear image of the Lorentz cone
#'
#' \code{ellips_semiax} takes as input an invertible matrix \code{A}
#' and returns the lengths of the semiaxes of the ellipsoidal cone given by
#' the image of the Lorentz cone under \code{A}.
#'
#' @param A matrix
#'
#' @return The output of \code{ellips_semiax(A)}, \code{A in Gl_d}, is a positive
#' vector \code{a in R^(d-1)} such that the cone \code{A*L}, where \code{L=\{x in R^d | x_d >= ||x||\}},
#' is isometric to the cone \code{\{x in R^d | x_d >= sum_(j=1)^(d-1) x_j^2/a_j^2\}}.
#'
#' @section See also:
#' \code{\link[conivol]{ellips_rbichibarsq}}
#'
#' Package: \code{\link[conivol]{conivol}}
#'
#' @note See \href{../doc/conic-intrinsic-volumes.html#ellips_cone}{this vignette}
#'       for further info.
#'
#' @examples
#' A <- matrix(c(2,3,5,7,11,13,17,19,23),3,3)
#' ellips_semiax(A)
#'
#' @export
#'
ellips_semiax <- function(A) {
    if (!is.matrix(A))
        stop("\n Input is not a matrix.")
    d <- dim(A)[1]
    if (dim(A)[2] != d)
        stop("\n Input matrix is not a square matrix.")
    if (qr(A)$rank != d)
        stop("\n Input matrix is not invertible.")

    eigs_AJAt <- sort(eigen( A %*% (c(rep(-1,d-1),1) * t(A)) )$values)
    return( sqrt(-eigs_AJAt[1:(d-1)]/eigs_AJAt[d]) )
}



# create the mosek input for the projection on ellipsoidal cone
#
.create_mosek_input_ellips <- function(A,z) {
    d <- dim(A)[1]

    Aext <- cbind( A, matrix(0,d,2), diag(-1,d) )

    mos_inp <- list(sense = "min")
    mos_inp$c     <- c( -z, 1, rep(0,d+1) )
    mos_inp$A     <- Matrix::Matrix( Aext, sparse=TRUE )
    mos_inp$bc    <- rbind(blc = rep(0,d),
                           buc = rep(0,d))
    mos_inp$bx    <- rbind(blx = c(rep(-Inf,d),0,1,rep(-Inf,d)),
                           bux = c(rep(Inf,d+1),1,rep(Inf,d)))
    mos_inp$cones <- matrix( list(
        "QUAD", c(d,1:(d-1)) ,
        "RQUAD", c(d+1,d+2,(d+3):(2*d+2))
    ), 2, 2 )
    rownames(mos_inp$cones) <- c("type", "sub")

    return(mos_inp)
}

.update_mosek_input_ellips <- function(mos_inp,z) {
    d <- length(z)
    mos_inp$c[1:d] <- -z
    return(mos_inp)
}


#' Sample from bivariate chi-bar-squared distribution of an ellipsoidal cone
#'
#' \code{ellips_rbichibarsq} generates an \code{n} by \code{2} matrix
#' such that the rows form iid samples from the bivariate chi-bar-squared
#' distribution of an ellipsoidal cone.
#'
#' @param n number of samples
#' @param A invertible matrix or positive vector
#' @param semiax logical; if \code{TRUE}, the algorithm will first compute the
#'               lengths of the semiaxes of the cone
#'
#' @return If \code{A} is an invertible matrix and \code{reduce==TRUE} (default),
#'         then the output of \code{ellips_rbichibarsq_gen(n,A)} is a list
#'         containing the following elements:
#' \itemize{
#'   \item \code{semiax}: a vector of the lengths of the semiaxes of the ellipsoidal
#'         cone defined by the matrix \code{A},
#'   \item \code{samples}: an \code{n} by \code{2} matrix whose rows form
#'         iid samples from the bivariate chi-bar-squared distribution with
#'         weights given by the intrinsic volumes of the ellipsoidal cone \code{A*L},
#'         where \code{L=\{x in R^d | x_d >= ||x||\}}.
#' }
#'         If \code{A} is a positive vector or \code{reduce==FALSE} then the output is only an
#'         \code{n} by \code{2} matrix whose rows form
#'         iid samples from the bivariate chi-bar-squared distribution with
#'         weights given by the intrinsic volumes of the ellipsoidal cone \code{A*L},
#'         or \code{diag(c(A,1))*L}, depending on whether \code{A} is an invertible
#'         matrix or a positive vector.
#'
#' @section See also:
#' \code{\link[conivol]{ellips_semiax}}
#'
#' @note See \href{../doc/conic-intrinsic-volumes.html#ellips_cone}{this vignette}
#'       for further info.
#'
#' @section See also:
#' Package: \code{\link[conivol]{conivol}}
#'
#' @examples
#' A <- matrix(c(2,3,5,7,11,13,17,19,23),3,3)
#' ellips_rbichibarsq(10,A)
#'
#' # the lengths of the semiaxes (directly):
#' ellips_semiax(A)
#'
#' # sampling without computing semiaxes first:
#' ellips_rbichibarsq(10,A,FALSE)
#'
#' @export
#'
ellips_rbichibarsq <- function(n,A, semiax = TRUE) {
    if (!is.matrix(A) && !is.vector(A))
        stop("\n A must be a vector or a matrix.")
    if ( is.null( dim(A) ) ) {         # A is a vector
        if ( !all( A>0 ) )
            stop("\n If A is a vector then it must be positive.")
        d <- length(A)+1
        A_ell <- diag( c(A,1) )
    } else {                           # A is a matrix
        d <- dim(A)[1]
        if (dim(A)[2] != d || qr(A)$rank != d)
            stop("\n If A is a matrix then it must be invertible.")
        if (semiax) {
            alpha <- conivol::ellips_semiax(A)
            A_ell <- diag( c(alpha,1) )
        } else {
            A_ell <- A
        }
    }
    if (!requireNamespace("Rmosek", quietly = TRUE))
        stop("\n Could not find package 'Rmosek'.")
    if (!requireNamespace("Matrix", quietly = TRUE))
        stop("\n Could not find package 'Matrix'.")

    opts <- list(verbose=0)
    mos_inp <- conivol:::.create_mosek_input_ellips(A_ell,rep(0,d))
    out <- matrix(0,n,2)

    for (i in 1:n) {
        z <- rnorm(d)
        mos_inp <- conivol:::.update_mosek_input_ellips(mos_inp,z)
        nrmprojsq <- sum(Rmosek::mosek(mos_inp,opts)$sol$itr$xx[(d+3):(2*d+2)]^2)
        out[i,1] <- nrmprojsq
        out[i,2] <- sum(z^2)-nrmprojsq
    }
    if (semiax)
        return( list( semiax=alpha, samples=out ) )
    else
        return(out)
}



#' Matrix representation of (products/duals of) Weyl chambers
#'
#' \code{weyl_matrix} computes a matrix representation of the (polars of)
#' Weyl chambers of finite reflection groups of type A, BC, and D.
#'
#' The dimensions and types are given
#' in the vectors \code{d} and \code{cone_type} (vectors must be of same lengths,
#' entries of \code{conetype} must be 'A', 'BC', 'D', 'Ap', 'BCp', or 'Dp').
#' If the length of the vectors is one, a single matrix is returned; if the
#' length of the vectors is greater
#' than one and \code{product==FALSE}, a list of matrices is returned;
#' if the length of the vectors is greater than one and \code{product==TRUE},
#' a single matrix representing the product cone is returned.
#'
#' @param d vector of dimensions; must be same length as \code{cone_type}.
#' @param cone_type vector of cone types; must be same length as \code{d}, the
#'             available types are as follows:
#'             \describe{
#'               \item{\code{"A"}:}{chamber of type A}
#'               \item{\code{"A_red"}:}{chamber of type A, reduced form}
#'               \item{\code{"BC"}:}{chamber of type BC}
#'               \item{\code{"D"}:}{chamber of type D}
#'               \item{\code{"Ap"}:}{polar of chamber of type A}
#'               \item{\code{"Ap_red"}:}{polar of chamber of type A, reduced form}
#'               \item{\code{"BCp"}:}{polar of chamber of type BC}
#'               \item{\code{"Dp"}:}{polar of chamber of type D}
#'             }
#' @param product logical; if \code{TRUE}, a representation of the product cone is returned.
#'
#' @return If \code{length(d)==1} or \code{(length(d)>1 & product==TRUE)}
#'         then a single matrix will be returned.
#'         If \code{(length(d)>1 & product==FALSE)} then a list of
#'         matrices will be returned.
#'
#' @section See also:
#' \code{\link[conivol]{weyl_ivols}}
#'
#' Package: \code{\link[conivol]{conivol}}
#'
#' @examples
#' weyl_matrix(5, "BC")
#' weyl_matrix(c(5,5), c("BC","Ap"))
#' weyl_matrix(c(5,5), c("BC","Ap"), product = TRUE)
#'
#' # testing the reduced cones
#' d <- 6
#' A <- weyl_matrix(d, "A")
#' A_red <- weyl_matrix(d, "A_red")
#' t(A) %*% A
#' round(t(A_red) %*% A_red, digits=14)
#'
#' Ap <- weyl_matrix(d, "Ap")
#' Ap_red <- weyl_matrix(d, "Ap_red")
#' t(Ap[,-c(1,2)]) %*% Ap[,-c(1,2)]
#' t(Ap_red) %*% Ap_red
#'
#' @export
#'
weyl_matrix <- function(d, cone_type, product = FALSE) {
    if (!requireNamespace("Matrix", quietly = TRUE))
        stop("\n Could not find package 'Matrix'.")
    if (length(d)!=length(cone_type))
        stop("Inputs d and cone_type must be of same length.")
    if (!all(cone_type %in% c("A","A_red","BC","D","Ap","Ap_red","BCp","Dp")))
        stop("Input cone_type must hav values 'A', 'A_red', 'BC', 'D', 'Ap', 'Ap_red', 'BCp', or 'Dp'.")
    if (!all(d[which(cone_type %in% c("A","A_red","Ap","Ap_red", "D","Dp"))] > 1))
        stop("Chambers of type 'D' and 'Dp' must be of dimension >1.")

    M <- list()
    for (i in 1:length(d)) {
        if (cone_type[i]=="A") {
            A <- rbind(0,diag(rep(-1,d[i])))
            diag(A) <- 1
        } else if (cone_type[i]=="A_red") {
            A <- rbind(0,diag(rep(-1,d[i])))
            diag(A) <- 1
            Q <- svd(A)$u
            A <- t(Q) %*% A
        } else if (cone_type[i]=="BC") {
            A <- diag(rep(-1,d[i]))
            diag(A[,-1]) <- 1
        } else if (cone_type[i]=="D") {
            A <- diag(rep(-1,d[i]))
            diag(A[,-1]) <- 1
            A[2,1] <- -1
        } else if (cone_type[i]=="Ap") {
            # A <- matrix(0,d[i]+1,d[i]+1)
            # A[lower.tri(A,diag=TRUE)] <- 1
            # A <- cbind(rep(-1,d[i]+1),A)

            A <- matrix(0,d[i]+1,d[i])
            I <- lower.tri(A)
            A[I] <- matrix(rep(1:d[i],d[i]+1),d[i]+1,d[i],byrow=TRUE)[I]
            I <- upper.tri(A,diag=TRUE)
            A[I] <- matrix(rep(-(d[i]:1),d[i]+1),d[i]+1,d[i],byrow=TRUE)[I]
            A <- cbind(rep(1,d[i]+1),A)
            A <- cbind(rep(-1,d[i]+1),A)
        } else if (cone_type[i]=="Ap_red") {
            A <- matrix(0,d[i]+1,d[i])
            I <- lower.tri(A)
            A[I] <- matrix(rep(1:d[i],d[i]+1),d[i]+1,d[i],byrow=TRUE)[I]
            I <- upper.tri(A,diag=TRUE)
            A[I] <- matrix(rep(-(d[i]:1),d[i]+1),d[i]+1,d[i],byrow=TRUE)[I]
            Q <- svd(A)$u
            A <- t(Q) %*% A
        } else if (cone_type[i]=="BCp") {
            A <- matrix(0,d[i],d[i])
            A[lower.tri(A,TRUE)] <- 1
        } else if (cone_type[i]=="Dp") {
            # A <- matrix(0,d[i],d[i])
            # A[lower.tri(A,TRUE)] <- 1
            # A[ ,c(1,2)] <- 1/2
            # A[1,2] <- -1/2
            A <- matrix(0,d[i],d[i])
            A[lower.tri(A,TRUE)] <- 1
            A[1,2] <- -1
        }
        M[[i]] <- A
    }
    if (length(d)==1)
        return(M[[1]])
    else if (product)
        return(as.matrix(Matrix::bdiag(M)))
    else return(M)
}


#' The conic intrinsic volumes of (products/duals of) Weyl chambers
#'
#' \code{weyl_ivols} computes the conic intrinsic volumes of (polars of)
#' Weyl chambers of finite reflection groups of type A, BC, D.
#'
#' The dimensions and types are given
#' in the vectors \code{d} and \code{cone_type} (vectors must be of same lengths,
#' entries of \code{conetype} must be 'A', 'BC', 'D', 'Ap', 'BCp', or 'Dp').
#' If the length of the vectors is one, a single vector is returned; if the
#' length of the vectors is greater
#' than one and \code{product==FALSE}, a list of vectors is returned;
#' if the length of the vectors is greater than one and \code{product==TRUE},
#' a single vector with the intrinsic volumes of the product cone is returned.
#'
#' @param d vector of dimensions; must be same length as \code{cone_type}.
#' @param cone_type vector of cone types; must be same length as \code{d}, the
#'             available types are as follows:
#'             \describe{
#'               \item{\code{"A"}:}{chamber of type A}
#'               \item{\code{"A_red"}:}{chamber of type A, reduced form}
#'               \item{\code{"BC"}:}{chamber of type BC}
#'               \item{\code{"D"}:}{chamber of type D}
#'               \item{\code{"Ap"}:}{polar of chamber of type A}
#'               \item{\code{"Ap_red"}:}{polar of chamber of type A, reduced form}
#'               \item{\code{"BCp"}:}{polar of chamber of type BC}
#'               \item{\code{"Dp"}:}{polar of chamber of type D}
#'             }
#' @param product logical; if \code{TRUE}, intrinsic volumes of product cone are returned.
#'
#' @return If \code{length(d)==1} or \code{(length(d)>1 & product==TRUE)}
#'         then a single vectors will be returned.
#'         If \code{(length(d)>1 & product==FALSE)} then a list of
#'         vectors will be returned.
#'
#' @section See also:
#' \code{\link[conivol]{weyl_matrix}}
#'
#' Package: \code{\link[conivol]{conivol}}
#'
#' @examples
#' weyl_ivols(5, "BC")
#' weyl_ivols(c(5,5), c("BC","Ap"))
#' weyl_ivols(c(5,5), c("BC","Ap"), product = TRUE)
#'
#' # testing the reduced cones
#' d <- 6
#' weyl_ivols(d, "A")
#' weyl_ivols(d, "A_red")
#'
#' weyl_ivols(d, "Ap")
#' weyl_ivols(d, "Ap_red")
#'
#' @export
#'
weyl_ivols <- function(d, cone_type, product = FALSE) {
    if (!requireNamespace("polynom", quietly = TRUE))
        stop("\n Could not find package 'polynom'.")
    if (length(d)!=length(cone_type))
        stop("Inputs d and cone_type must be of same length.")
    if (!all(cone_type %in% c("A","A_red","BC","D","Ap","Ap_red","BCp","Dp")))
        stop("Input cone_type must hav values 'A', 'A_red', 'BC', 'D', 'Ap', 'Ap_red', 'BCp', or 'Dp'.")
    if (!all(d[which(cone_type %in% c("D","Dp"))] > 1))
        stop("Chambers of type 'D' and 'Dp' must be of dimension >1.")

    V <- list()
    for (i in 1:length(d)) {
        if (cone_type[i] %in% c("A","A_red","Ap","Ap_red")) {
            v <- as.vector(polynom::poly.calc( -(1:d[i]) ))
            if (cone_type[i] %in% c("A","Ap"))
                v <- c(0,v)
        } else if (cone_type[i] %in% c("BC","BCp")) {
            v <- as.vector(polynom::poly.calc( -(2*(1:d[i])-1) ))
        } else if (cone_type[i] %in% c("D","Dp")) {
            v <- as.vector(polynom::poly.calc( -c( d[i]-1, 2*(1:(d[i]-1))-1 ) ))
        }
        if (cone_type[i] %in% c("Ap","Ap_red","BCp","Dp"))
            v <- rev(v)
        V[[i]] <- v/sum(v)
    }
    if (length(d)==1)
        return(V[[1]])
    else if (product)
        return(conivol::prod_ivols(V))
    else return(V)
}
