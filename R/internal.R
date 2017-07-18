
# adjust tolerance level for some numeric computations
#
.adj_tol <- 10*.Machine$double.eps



# just a small function to use in the other function for testing whether an
# input is a genuine vector of probabilities
#
.test_vector <- function(v){
    if (any(v<0 & sapply(v, function(t){!isTRUE(all.equal(t,0,tolerance=.adj_tol))}) )
        | any(v>1) | !isTRUE(all.equal(sum(v),1)))
        stop("Vector v must be a discrete probability distribution.")
}


# computing the log-likelihood (normalized by sample size)
#
.comp_loglike <- function(v, D){
    d <- length(v)-1
    return(
        D$prop_pol  * log(v[1]) +
        sum( 1/D$n * log( colSums( D$dens * v[2:d] ) ) ) +
        D$prop_prim * log(v[d+1])
    )
}







