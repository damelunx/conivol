#' Compute intrinsic volumes of product cone from intrinsic volumes of components.
#'
#' \code{comp_ivols_product} computes the intrinsic volumes of a product cone
#' from the intrinsic volumes of ites components. That is, it returns the
#' convolution of the vectors of intrinsic volumes given in \code{V}.
#'
#' @param V list of vectors of intrinsic volumes
#'
#' @return The output of \code{comp_ivols_product(V)} is the vector that is
#'         obtained from convolving the components of \code{V}.
#'
#' @examples
#' comp_ivols_product(list( c(0.5,0.5), c(0.1,0.4,0.5) ))
#' comp_ivols_product(list( c(0.5,0.5), c(0.5,0.5), c(0.5, 0.5) ))
#'
comp_ivols_product <- function(V) {
    lapply(V, .test_vector)
    for (i in 2:length(V))
        V[[1]] <- convolve(V[[1]],rev(V[[i]]),type="o")
    # test which one are numerically zero, then set them equal to zero
    # (to avoid negative entries)
    I <- which(sapply(V[[1]], function(t){isTRUE(all.equal(t,0,tolerance=.adj_tol))}))
    V[[1]][I]=0
    return(V[[1]])
}



