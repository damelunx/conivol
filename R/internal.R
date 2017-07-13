
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



# takes the two-column matrix of samples of the bivariate chi-bar-squared
# distribution and outputs a list of three elements:
# (1) an integer vector with values in {0,1,2} that indicates whether
#     x\not\in (C\cup C°) (0) , x\in C (1) , x\in C° (2)
# (2) a matrix with the density values of the primal projection
# (3) a matrix with the density values of the primal projection
#
.prepare_proj_data <- function(d,m_input,log=FALSE) {
    # d is the dimension of the normal distribution used in the sampling
    # m_input should be a two-column matrix with sample from the bivariate chi-bar-squared distribution

    n <- dim(m_input)[1]
    ind <- rep(0,n)

    I1 <- which( sapply(m_input[ ,1], function(t){isTRUE(all.equal(t,0,tolerance=.adj_tol))}) )
    I2 <- which( sapply(m_input[ ,2], function(t){isTRUE(all.equal(t,0,tolerance=.adj_tol))}) )
    ind[I1] <- 2
    ind[I2] <- 1

    distrib_prim <- matrix( sapply( m_input[ ,1], function(x){dchisq(x,1:d,log=log)} ), n, d, byrow=TRUE)
    distrib_pol  <- matrix( sapply( m_input[ ,2], function(x){dchisq(x,1:d,log=log)} ), n, d, byrow=TRUE)

    return(list(ind=ind, prim=distrib_prim, pol=distrib_pol))
}



# creating the starting point for the EM algorithm
#
.init_v <- function(d,mode=0,delta=d/2,var=d/4) {
    # mode==0: uniform distribution
    # mode==1: discretized normal distribution
    # mode==2: circular cone fitting statdim
    # mode==3: circular cone fitting variance
    # mode==4: circular cone fitting statdim and variance (geometric mean)
    if (mode==1) {
        v <- sapply( 0:d, function(k){pnorm((k+0.5-delta)/sqrt(var)) - pnorm((k-0.5-delta)/sqrt(var))})
        return(v/sum(v))
    } else if (mode==2) {
        alpha <- asin(sqrt(delta/d))
        return(circ_ivol(d,alpha))
    } else if (mode==3) {
        alpha <- asin(sqrt(2*var/(d-2)))/2
        if ((alpha<pi/4 && delta>d/2) || (alpha>pi/4 && delta<d/2))
            alpha <- pi/2-alpha;
        return(circ_ivol(d,alpha))
    } else if (mode==4) {
        alpha1 <- asin(sqrt(delta/d))
        alpha2 <- asin(sqrt(2*var/(d-2)))/2
        if ((alpha2<pi/4 && delta>d/2) || (alpha2>pi/4 && delta<d/2))
            alpha2 <- pi/2-alpha2;
        v1 <- circ_ivol(d,alpha1)
        v2 <- circ_ivol(d,alpha2)
        v <- sqrt(v1*v2)
        return(v/sum(v))
    } else {
        return(1/(d+1)*rep(1,d+1))
    }
}



# computing the log-likelihoods
#
.comp_loglike <- function(v, D){
    ind <- D$ind
    prim <- D$prim
    pol <- D$pol
    d <- dim(prim)[2]

    I0 <- which(D$ind==0)
    I1 <- which(D$ind==0 | D$ind==1)
    I2 <- which(D$ind==0 | D$ind==2)

    n <- dim(v)[1]
    loglike0 <- vector("double", n)
    loglike1 <- vector("double", n)
    loglike2 <- vector("double", n)

    for (i in 1:n) {
        v0 <- v[i,2:d]/(1-v[i,1]-v[i,d+1])
        v1 <- v[i,2:(d+1)]/(1-v[i,1])
        v2 <- v[i,1:d]/(1-v[i,d+1])

        denom0 <- rowSums( sweep( D$prim[I0,1:(d-1)] * D$pol[I0,rev(1:(d-1))] , MARGIN=2,v0,"*") )
        denom1 <- rowSums( sweep( D$prim[I1,1:d] , MARGIN=2,v1,"*") )
        denom2 <- rowSums( sweep( D$pol[I2,1:d]  , MARGIN=2,v2,"*") )

        loglike0[i] <- sum(log(denom0))
        loglike1[i] <- sum(log(denom1))
        loglike2[i] <- sum(log(denom2))
    }

    return(list( loglike0=loglike0, loglike1=loglike1, loglike2=loglike2 ))
}







