
# adjust tolerance level for some numeric computations
#
.adj_tol <- 10*.Machine$double.eps



# just a small function to use in the other function for testing whether an
# input is a genuine vector of probabilities
#
.test_vector <- function(v){
    if (any(v<0 & sapply(v, function(t){!isTRUE(all.equal(t,0,tolerance=conivol:::.adj_tol))}) )
        | any(v>1) | !isTRUE(all.equal(sum(v),1)))
        stop("Vector v must be a discrete probability distribution.")
}




