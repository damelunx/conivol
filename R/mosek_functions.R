
# create the input for mosek
#
.create_mosek_input <- function(d, const, v0, vd, selfdual) {
    mos_inp0 <- list()
    mos_inp1 <- list()
    mos_inp2 <- list()

    mos_inp0$sense <- "max"
    mos_inp1$sense <- "max"
    mos_inp2$sense <- "max"

    # setting optimizer
    mos_inp0$c <- c(rep(0,d-1))
    mos_inp1$c <- c(rep(0,d))
    mos_inp2$c <- c(rep(0,d))
    opro0 <- matrix(list(), nrow=5, ncol=d-1)
    opro1 <- matrix(list(), nrow=5, ncol=d)
    opro2 <- matrix(list(), nrow=5, ncol=d)
    rownames(opro0) <- c("type","j","f","g","h")
    rownames(opro1) <- c("type","j","f","g","h")
    rownames(opro2) <- c("type","j","f","g","h")
    for (i in 1:(d-1)) {
        opro0[ ,i] <- list("LOG", i, const$c0[i], 1.0, 0.0)
    }
    for (i in 1:d) {
        opro1[ ,i] <- list("LOG", i, const$c1[i], 1.0, 0.0)
        opro2[ ,i] <- list("LOG", i, const$c2[i], 1.0, 0.0)
    }
    mos_inp0$scopt <- list(opro=opro0)
    mos_inp1$scopt <- list(opro=opro1)
    mos_inp2$scopt <- list(opro=opro2)


    # variable constraints
    blx0 <- rep(0, d-1)               # v_i >= 0
    blx1 <- rep(0, d)
    blx2 <- rep(0, d)
    bux0 <- rep(.5/(1-v0-vd), d-1)    # v_i <= ...
    bux1 <- rep(.5/(1-v0), d)
    bux2 <- rep(.5/(1-vd), d)
    mos_inp0$bx <- rbind(blx0, bux0)
    mos_inp1$bx <- rbind(blx1, bux1)
    mos_inp2$bx <- rbind(blx2, bux2)

    # constraint matrix:
    if (!selfdual) {
        nrow <- 2
    } else {
        nrow <- 2+floor((d-1)/2)
    }
    A0 <- Matrix(c( rep_len(1:0,d-1), rep_len(0:1,d-1) , rep(0,(nrow-2)*(d-1)) ), nrow=nrow, byrow=TRUE, sparse=TRUE)
    A1 <- Matrix(c( rep_len(1:0,d),   rep_len(0:1,d)   , rep(0,(nrow-2)*d)     ), nrow=nrow, byrow=TRUE, sparse=TRUE)
    A2 <- Matrix(c( rep_len(1:0,d),   rep_len(0:1,d)   , rep(0,(nrow-2)*d)     ), nrow=nrow, byrow=TRUE, sparse=TRUE)
    if (nrow>2)
        for (i in 1:(nrow-2)) {
            A0[2+i,c(i,d-i)] <- c(1,-1)
            A1[2+i,c(i,d-i)] <- c(1,-1)
            A2[2+i,c(i+1,d-i+1)] <- c(1,-1)
        }
    mos_inp0$A <- A0
    mos_inp1$A <- A1
    mos_inp2$A <- A2

    # constraint rhs:
    if (d%%2) {
        blc0 <- c( (.5-vd)/(1-v0-vd) ,    (.5-v0)/(1-v0-vd) , rep(0,nrow-2) )
        blc1 <- c(      .5/(1-v0)    ,    (.5-v0)/(1-v0)    , rep(0,nrow-2) )
        blc2 <- c(      .5/(1-vd)    ,    (.5-vd)/(1-vd)    , rep(0,nrow-2) )
    } else {
        blc0 <- c(      .5/(1-v0-vd) , (.5-v0-vd)/(1-v0-vd) , rep(0,nrow-2) )
        blc1 <- c(      .5/(1-v0)    ,    (.5-v0)/(1-v0)    , rep(0,nrow-2) )
        blc2 <- c(   .5-vd/(1-vd)    ,         .5/(1-vd)    , rep(0,nrow-2) )
    }
    # since all constraints are equalities, the upper constraints match the lower ones
    buc0 <- blc0
    buc1 <- blc0
    buc2 <- blc0
    mos_inp0$bc <- rbind(blc0, buc0);
    mos_inp1$bc <- rbind(blc1, buc1);
    mos_inp2$bc <- rbind(blc2, buc2);

    return( list(mode0=mos_inp0, mode1=mos_inp1, mode2=mos_inp2) )
}


.update_mosek_input <- function(mos_inp, d, c0, c1, c2, v0, vd) {
    # update optimizer
    mos_inp$mode0$scopt$opro[3, ] <- c0
    mos_inp$mode1$scopt$opro[3, ] <- c1
    mos_inp$mode2$scopt$opro[3, ] <- c2

    # update variable constraints
    mos_inp$mode0$bx[2, ] <- rep(.5/(1-v0-vd), d-1)
    mos_inp$mode1$bx[2, ] <- rep(.5/(1-v0), d)
    mos_inp$mode2$bx[2, ] <- rep(.5/(1-vd), d)

    # update constraint rhs:
    if (d%%2) {
        mos_inp$mode0$bc[1,c(1,2)] <- c( (.5-vd)/(1-v0-vd) ,    (.5-v0)/(1-v0-vd) )
        mos_inp$mode1$bc[1,c(1,2)] <- c(      .5/(1-v0)    ,    (.5-v0)/(1-v0)    )
        mos_inp$mode2$bc[1,c(1,2)] <- c(      .5/(1-vd)    ,    (.5-vd)/(1-vd)    )
    } else {
        mos_inp$mode0$bc[1,c(1,2)] <- c(      .5/(1-v0-vd) , (.5-v0-vd)/(1-v0-vd) )
        mos_inp$mode1$bc[1,c(1,2)] <- c(      .5/(1-v0)    ,    (.5-v0)/(1-v0)    )
        mos_inp$mode2$bc[1,c(1,2)] <- c(   .5-vd/(1-vd)    ,         .5/(1-vd)    )
    }
    # match the upper constraints
    mos_inp$mode0$bc[2,c(1,2)] <- mos_inp$mode0$bc[1,c(1,2)]
    mos_inp$mode1$bc[2,c(1,2)] <- mos_inp$mode1$bc[1,c(1,2)]
    mos_inp$mode2$bc[2,c(1,2)] <- mos_inp$mode2$bc[1,c(1,2)]

    return(mos_inp)
}


