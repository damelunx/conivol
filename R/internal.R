
# adjust tolerance level for some numeric computations
#
.conivol_adj_tol <- 10*.Machine$double.eps



# just a small function to use in the other function for testing whether an
# input is a genuine vector of probabilities
#
.conivol_test_vector <- function(v){
    if (any(v<0 & sapply(v, function(t){!isTRUE(all.equal(t,0,tolerance=.conivol_adj_tol))}) )
        | any(v>1) | !isTRUE(all.equal(sum(v),1)))
        stop("Vector v must be a discrete probability distribution.")
}


# computing the log-likelihood (normalized by sample size)
#
.conivol_comp_loglike <- function(v, D){
    d <- length(v)-1
    return(
        D$prop_pol  * log(v[1]) +
        sum( 1/D$n * log( colSums( D$dens * v[2:d] ) ) ) +
        D$prop_prim * log(v[d+1])
    )
}



# create the mosek input for the projection
#
.conivol_create_mosek_input_polyh <- function(A,y) {
    m <- dim(A)[1]
    n <- dim(A)[2]

    Aext <- cbind( A, matrix(0,m,2), diag(-1,m,m) )

    mos_inp <- list(sense = "min")
    mos_inp$c     <- c( -as.vector(t(A)%*%y), 1, rep(0,m+1) )
    mos_inp$A     <- Matrix::Matrix( as.vector(Aext), nrow=m, byrow=FALSE, sparse=TRUE )
    mos_inp$bc    <- rbind(blc = rep(0,m),
                           buc = rep(0,m))
    mos_inp$bx    <- rbind(blx = c(rep(0,n+1),1,rep(-Inf,m)),
                           bux = c(rep(Inf,n+1),1,rep(Inf,m)))
    mos_inp$cones <- matrix( list("RQUAD", c(n+1,n+2,(n+3):(m+n+2))), 2, 1 )
    rownames(mos_inp$cones) <- c("type", "sub")

    return(mos_inp)
}

.conivol_update_mosek_input_polyh <- function(mos_inp,y) {
    m <- length(y)
    n <- dim(mos_inp$A)[2]-m-2
    mos_inp$c[1:n] <- -as.vector( t(mos_inp$A[1:m,1:n]) %*% y )
    return(mos_inp)
}




