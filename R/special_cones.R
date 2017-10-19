#' The conic intrinsic volumes of (products of) circular cones.
#'
#' \code{circ_ivols} computes the conic intrinsic volumes of circular cones,
#' whose dimensions and angles are given in the vectors \code{d} and
#' \code{alpha} (vectors must be of same lengths); if the length of the vectors
#' is one, a single vector is returned; if the length of the vectors is greater
#' than one and \code{prcoduct==FALSE}, a list of vectors is returned;
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
#' \code{\link[conivol]{circ_rbichibarsq}}
#'
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
            v[1] <- exp( lgamma(d[i]/2)-lgamma((d[i]+1)/2)-lgamma(1/2)+log(d[i]-1)-log(2)+
                             log( integrate(function(x){sin(x)^(d[i]-2)},0,pi/2-alpha[i])$value ) )

            v[d[i]+1] <- exp( lgamma(d[i]/2)-lgamma((d[i]+1)/2)-lgamma(1/2)+log(d[i]-1)-log(2)+
                                  log( integrate(function(x){sin(x)^(d[i]-2)},0,alpha[i])$value ) )

            k <- 1:(d[i]-1)
            v[2:d[i]] <- exp( lgamma(d[i]/2)-lgamma((k+1)/2)-lgamma((d[i]-k+1)/2)+
                                  (k-1)*log(sin(alpha[i]))+(d[i]-k-1)*log(cos(alpha[i]))-log(2) )
        }
        V[[i]] <- v
    }

    if (length(d)==1)
        return(V[[1]])
    else if (product)
        return(conivol::comp_ivols_product(V))
    else return(V)
}



#' Matrix representation of (products/duals of) Weyl chambers.
#'
#' \code{weyl_matrix} computes a matrix representation of the (polars of)
#' Weyl chambers of finite reflection groups of type A, BC, and D.
#' The dimensions and types are given
#' in the vectors \code{d} and \code{cone_type} (vectors must be of same lengths,
#' entries of \code{conetype} must be 'A', 'BC', 'D', 'Ap', 'BCp', or 'Dp').
#' If the length of the vectors is one, a single matrix is returned; if the
#' length of the vectors is greater
#' than one and \code{prcoduct==FALSE}, a list of matrices is returned;
#' if the length of the vectors is greater than one and \code{product==TRUE},
#' a single matrix representing the product cone is returned.
#'
#' @param d vector of dimensions; must be same length as \code{cone_type}.
#' @param cone_type vector of cone types; must be same length as \code{d},
#'             \describe{
#'               \item{\code{cone_type=="A"}:}{chamber of type A}
#'               \item{\code{cone_type=="BC"}:}{chamber of type BC}
#'               \item{\code{cone_type=="D"}:}{chamber of type D}
#'               \item{\code{cone_type=="Ap"}:}{polar of chamber of type A}
#'               \item{\code{cone_type=="BCp"}:}{polar of chamber of type BC}
#'               \item{\code{cone_type=="Dp"}:}{polar of chamber of type D}
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
#' @export
#'
weyl_matrix <- function(d, cone_type, product = FALSE) {
    if (!requireNamespace("Matrix", quietly = TRUE))
        stop("\n Could not find package 'Matrix'.")
    if (length(d)!=length(cone_type))
        stop("Inputs d and cone_type must be of same length.")
    if (!all(cone_type %in% c("A","BC","D","Ap","BCp","Dp")))
        stop("Input cone_type must hav values 'A', 'BC', 'D', 'Ap', 'BCp', or 'Dp'.")
    if (!all(d[which(cone_type %in% c("D","Dp"))] > 1))
        stop("Chambers of type 'D' and 'Dp' must be of dimension >1.")

    M <- list()
    for (i in 1:length(d)) {
        if (cone_type[i]=="A") {
            A <- rbind(0,diag(rep(-1,d[i])))
            diag(A) <- 1
        } else if (cone_type[i]=="BC") {
            A <- diag(rep(-1,d[i]))
            diag(A[,-1]) <- 1
        } else if (cone_type[i]=="D") {
            A <- diag(rep(-1,d[i]))
            diag(A[,-1]) <- 1
            A[2,1] <- -1
        } else if (cone_type[i]=="Ap") {
            A <- matrix(0,d[i]+1,d[i]+1)
            A[lower.tri(A)] <- 1
        } else if (cone_type[i]=="BCp") {
            A <- matrix(0,d[i],d[i])
            A[lower.tri(A,TRUE)] <- 1
        } else if (cone_type[i]=="Dp") {
            A <- matrix(0,d[i],d[i])
            A[lower.tri(A,TRUE)] <- 1
            A[ ,c(1,2)] <- 1/2
            A[1,2] <- -1/2
        }
        if (cone_type[i] %in% c("A","Ap")) {
            tmp <- rbind(0,diag(rep(-1,d)))
            diag(tmp) <- 1
            Q <- svd(tmp)$u
            A <- t(Q) %*% A
        }
        M[[i]] <- A
    }
    if (length(d)==1)
        return(M[[1]])
    else if (product)
        return(as.matrix(Matrix::bdiag(M)))
    else return(V)
}


#' The conic intrinsic volumes of (products/duals of) Weyl chambers.
#'
#' \code{weyl_ivols} computes the conic intrinsic volumes of (polars of)
#' Weyl chambers of finite reflection groups of type A, BC, D.
#' The dimensions and types are given
#' in the vectors \code{d} and \code{cone_type} (vectors must be of same lengths,
#' entries of \code{conetype} must be 'A', 'BC', 'D', 'Ap', 'BCp', or 'Dp').
#' If the length of the vectors is one, a single vector is returned; if the
#' length of the vectors is greater
#' than one and \code{prcoduct==FALSE}, a list of vectors is returned;
#' if the length of the vectors is greater than one and \code{product==TRUE},
#' a single vector with the intrinsic volumes of the product cone is returned.
#'
#' @param d vector of dimensions; must be same length as \code{cone_type}.
#' @param cone_type vector of cone types; must be same length as \code{d},
#'             \describe{
#'               \item{\code{cone_type=="A"}:}{chamber of type A}
#'               \item{\code{cone_type=="BC"}:}{chamber of type BC}
#'               \item{\code{cone_type=="D"}:}{chamber of type D}
#'               \item{\code{cone_type=="Ap"}:}{polar of chamber of type A}
#'               \item{\code{cone_type=="BCp"}:}{polar of chamber of type BC}
#'               \item{\code{cone_type=="Dp"}:}{polar of chamber of type D}
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
#' @export
#'
weyl_ivols <- function(d, cone_type, product = FALSE) {
    if (!requireNamespace("polynom", quietly = TRUE))
        stop("\n Could not find package 'polynom'.")
    if (length(d)!=length(cone_type))
        stop("Inputs d and cone_type must be of same length.")
    if (!all(cone_type %in% c("A","BC","D","Ap","BCp","Dp")))
        stop("Input cone_type must hav values 'A', 'BC', 'D', 'Ap', 'BCp', or 'Dp'.")
    if (!all(d[which(cone_type %in% c("D","Dp"))] > 1))
        stop("Chambers of type 'D' and 'Dp' must be of dimension >1.")

    V <- list()
    for (i in 1:length(d)) {
        if (cone_type[i] %in% c("A","Ap")) {
            v <- as.vector(polynom::poly.calc( -(1:d[i]) ))
        } else if (cone_type[i] %in% c("BC","BCp")) {
            v <- as.vector(polynom::poly.calc( -(2*(1:d[i])-1) ))
        } else if (cone_type[i] %in% c("D","Dp")) {
            v <- as.vector(polynom::poly.calc( -c( d-1, 2*(1:d[i])-3 ) ))
        }
        if (cone_type[i] %in% c("Ap","BCp","Dp"))
            v <- rev(v)
        V[[i]] <- v/sum(v)
    }
    if (length(d)==1)
        return(V[[1]])
    else if (product)
        return(conivol::comp_ivols_product(V))
    else return(V)
}
