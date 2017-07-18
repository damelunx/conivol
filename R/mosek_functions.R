
# create the input for mosek for EM step
#
.create_mosek_input_EM <- function(const,extrap_pol,extrap_prim,selfdual) {
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
        A <- Matrix(c( rep_len(1:0,d+1), rep_len(0:1,d+1) , rep(0,(nrow-2)*(d+1)) ), nrow=nrow, byrow=TRUE, sparse=TRUE)
        if (nrow>2)
            for (i in 1:(nrow-2))
                A[2+i,c(i,d+2-i)] <- c(1,-1)
    } else if (extrap_pol & !extrap_prim) {
        A <- Matrix(c( rep_len(1:0,d), rep_len(0:1,d) , rep(0,(nrow-2)*d) ), nrow=nrow, byrow=TRUE, sparse=TRUE)
        if (nrow>2)
            for (i in 1:(nrow-2))
                A[2+i,c(1+i,d+2-i)] <- c(1,-1)
    } else if (!extrap_pol & extrap_prim) {
        A <- Matrix(c( rep_len(1:0,d), rep_len(0:1,d) , rep(0,(nrow-2)*d) ), nrow=nrow, byrow=TRUE, sparse=TRUE)
        if (nrow>2)
            for (i in 1:(nrow-2))
                A[2+i,c(i,d+2-i-1)] <- c(1,-1)
    } else if (extrap_pol & extrap_prim) {
        A <- Matrix(c( rep_len(1:0,d-1), rep_len(0:1,d-1) , rep(0,(nrow-2)*(d-1)) ), nrow=nrow, byrow=TRUE, sparse=TRUE)
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


.update_mosek_input_EM <- function(mos_inp,const) {
    mos_inp$scopt$opro[3, ] <- const
    return(mos_inp)
}


