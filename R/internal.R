
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




# create the input for mosek for EM step
#
.conivol_create_mosek_input_EM <- function(const,extrap_pol,extrap_prim,selfdual) {
    d <- length(const)-1

    mos_inp <- list()
    mos_inp$sense <- "max"

    # setting optimizer
    if (!extrap_pol & !extrap_prim) {
        mos_inp$c <- c(rep(0,d+1))
        opro <- matrix(list(), nrow=5, ncol=d+1)
        rownames(opro) <- c("type","j","f","g","h")
        for (i in 1:(d+1)) {
            opro[ ,i] <- list("LOG", i, const[i], 1.0, 0.0)
        }
    } else if ((extrap_pol & !extrap_prim) | (!extrap_pol & extrap_prim)) {
        mos_inp$c <- c(rep(0,d))
        opro <- matrix(list(), nrow=5, ncol=d)
        rownames(opro) <- c("type","j","f","g","h")
        for (i in 1:d) {
            opro[ ,i] <- list("LOG", i, const[i], 1.0, 0.0)
        }
    } else if (extrap_pol & extrap_prim) {
        mos_inp$c <- c(rep(0,d-1))
        opro <- matrix(list(), nrow=5, ncol=d-1)
        rownames(opro) <- c("type","j","f","g","h")
        for (i in 1:(d-1)) {
            opro[ ,i] <- list("LOG", i, const[i], 1.0, 0.0)
        }
    }
    mos_inp$scopt <- list(opro=opro)

    # variable constraints
    if (!extrap_pol & !extrap_prim) {
        blx <- rep(0, d+1)            # v_i >= 0
        bux <- rep(.5, d+1)           # v_i <= 0.5
    } else if ((extrap_pol & !extrap_prim) | (!extrap_pol & extrap_prim)) {
        blx <- rep(0, d)              # v_i >= 0
        bux <- rep(.5, d)             # v_i <= 0.5
    } else if (extrap_pol & extrap_prim) {
        blx <- rep(0, d-1)            # v_i >= 0
        bux <- rep(.5, d-1)           # v_i <= 0.5
    }
    mos_inp$bx <- rbind(blx, bux)

    # constraint matrix:
    if (!selfdual) {
        nrow <- 2
    } else {
        if (!extrap_pol & !extrap_prim)
            nrow <- 2+floor((d+1)/2)
        else
            nrow <- 2+floor((d-1)/2)
    }
    if (!extrap_pol & !extrap_prim) {
        A <- Matrix::Matrix(c( rep_len(1:0,d+1), rep_len(0:1,d+1) , rep(0,(nrow-2)*(d+1)) ), nrow=nrow, byrow=TRUE, sparse=TRUE)
        if (nrow>2)
            for (i in 1:(nrow-2))
                A[2+i,c(i,d+2-i)] <- c(1,-1)
    } else if (extrap_pol & !extrap_prim) {
        A <- Matrix::Matrix(c( rep_len(1:0,d), rep_len(0:1,d) , rep(0,(nrow-2)*d) ), nrow=nrow, byrow=TRUE, sparse=TRUE)
        if (nrow>2)
            for (i in 1:(nrow-2))
                A[2+i,c(1+i,d+2-i)] <- c(1,-1)
    } else if (!extrap_pol & extrap_prim) {
        A <- Matrix::Matrix(c( rep_len(1:0,d), rep_len(0:1,d) , rep(0,(nrow-2)*d) ), nrow=nrow, byrow=TRUE, sparse=TRUE)
        if (nrow>2)
            for (i in 1:(nrow-2))
                A[2+i,c(i,d+2-i-1)] <- c(1,-1)
    } else if (extrap_pol & extrap_prim) {
        A <- Matrix::Matrix(c( rep_len(1:0,d-1), rep_len(0:1,d-1) , rep(0,(nrow-2)*(d-1)) ), nrow=nrow, byrow=TRUE, sparse=TRUE)
        if (nrow>2)
            for (i in 1:(nrow-2))
                A[2+i,c(i,d+2-i)] <- c(1,-1)
    }
    mos_inp$A <- A

    # constraint rhs:
    buc <- c( .5, .5, rep(0,nrow-2) )
    if (!extrap_pol & !extrap_prim)
        blc <- buc
    else if (extrap_pol & !extrap_prim)
        blc <- c( 0, .5, rep(0,nrow-2) )
    else if (!extrap_pol & extrap_prim) {
        if (d%%2 == 0)
            blc <- c( 0, .5, rep(0,nrow-2) )
        else
            blc <- c( .5, 0, rep(0,nrow-2) )
    } else if (extrap_pol & extrap_prim) {
        if (d%%2 == 0)
            blc <- c( 0, .5, rep(0,nrow-2) )
        else
            blc <- c( 0, 0, rep(0,nrow-2) )
    }

    mos_inp$bc <- rbind(blc, buc);

    return( mos_inp )
}


.conivol_update_mosek_input_EM <- function(mos_inp,const) {
    mos_inp$scopt$opro[3, ] <- const
    return(mos_inp)
}


