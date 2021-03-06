## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "conic-intrinsic-volumes_figures/"
)

## ----load-pkgs, include=FALSE--------------------------------------------
library(conivol)
library(tidyverse)
library(logcondiscr)
library(partitions)
library(Rmisc)

## ----polyh-red-----------------------------------------------------------
A <- matrix(c(-(1:4),1:24),4,7); A
polyh_reduce_gen(A)
polyh_reduce_ineq(A)

## ----weyl-comps1---------------------------------------------------------
A <- weyl_matrix(5,"A")
A_red <- weyl_matrix(5,"A_red")
list( A=A, A_red=A_red )

## ----weyl-comps2---------------------------------------------------------
t(A) %*% A
round( t(A_red) %*% A_red ,digits=14)

## ----ellips-semiax1------------------------------------------------------
d <- 5
ellips_semiax( diag(d:1) )

## ----ellips-semiax2------------------------------------------------------
Q <- svd( matrix(rnorm(d^2),d,d) )$u  # find random rotation
round( t(Q) %*% Q, 14 )               # test orthogonality
ellips_semiax( Q %*% diag(d:1) )      # compute semiaxes of rotated cone

## ----prod-ivols----------------------------------------------------------
v_halfline <- c(1/2,1/2); v_halfline
2^4 * prod_ivols( list(v_halfline, v_halfline, v_halfline, v_halfline) )

## ----read-canc-data, warning=FALSE, message=FALSE------------------------
address <- "http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt"
A_canc <- readr::read_tsv(address)[ , 1:33]
A <- as.matrix(A_canc[ , 4:33])

## ----sum-canc-data-------------------------------------------------------
dim(A)
A[1:5,1:5]
all(A>=0)
colSums(A)

## ----samp-canc-data------------------------------------------------------
n <- 1e5
S <- polyh_rivols_gen(n,A)         # sampling from intrinsic volumes distribution

str(S)
linC <- S$linC
dimC <- S$dimC
msamp <- S$multsamp

## ----plot-canc-samp, fig.width = 7---------------------------------------
tib_plot <- tibble( k=linC:dimC, value=msamp[1+linC:dimC]/n )
ggplot(tib_plot, aes(x=k, y=value))      + geom_line() + theme_bw()
ggplot(tib_plot, aes(x=k, y=log(value))) + geom_line() + theme_bw()

## ----canc-bayes, fig.width = 7-------------------------------------------
bayes_est <- polyh_bayes( msamp, dimC, linC )
tib_plot_bay <- bayes_est$post_samp(1e4) %>%
    as_tibble() %>%
    `colnames<-`(paste0(rep("V",dimC-linC+1),as.character(linC:dimC))) %>%
    gather(factor_key=TRUE)
ggplot(tib_plot_bay, aes(x=key, y=value)) +
    geom_boxplot() + theme_bw() +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank())

## ----samp-canc-data-small, fig.width = 7---------------------------------
n_sm <- 3e2
set.seed(1111)
S_sm <- polyh_rivols_gen(n_sm,A)
msamp_sm <- S_sm$multsamp

# point estimate:
tib_plot_sm <- tibble( k=linC:dimC, value=msamp_sm[1+linC:dimC]/n_sm )
ggplot(tib_plot_sm, aes(x=k, y=value)) + geom_line() + theme_bw() +
    geom_line(data=tib_plot, aes(x=k, y=value), linetype="dashed", color="red")
ggplot(tib_plot_sm, aes(x=k, y=log(value))) + geom_line() + theme_bw() +
    geom_line(data=tib_plot, aes(x=k, y=log(value)), linetype="dashed", color="red")
# Bayes posterior:
bayes_est_sm <- polyh_bayes( msamp_sm, dimC, linC )
tib_plot_bay_sm <- bayes_est_sm$post_samp(1e4) %>%
    as_tibble() %>%
    `colnames<-`(paste0(rep("V",dimC-linC+1),as.character(linC:dimC))) %>%
    gather(factor_key=TRUE)
ggplot(tib_plot_bay_sm, aes(x=key, y=value)) +
    geom_boxplot() + theme_bw() +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
    geom_line(data=tib_plot, aes(x=k+1, y=value), linetype="dashed", color="red")

## ----samp-canc-data-logconc, fig.width = 7-------------------------------
logconc_MLE <- logConDiscrMLE(S_sm$samples, output=FALSE)
log_est <- rep(0,dimC-linC+1)
log_est[logconc_MLE$x+linC+1] <- exp(logconc_MLE$psi)

# log-concavity improved point estimate:
tib_plot_logc <- tibble( k=linC:dimC, value=log_est )
ggplot(tib_plot_logc, aes(x=k, y=value)) + geom_line() + theme_bw() +
    geom_line(data=tib_plot_sm, aes(x=k, y=value), linetype="dashed", color="blue") +
    geom_line(data=tib_plot, aes(x=k, y=value), linetype="dashed", color="red")
ggplot(tib_plot_logc, aes(x=k, y=log(value))) + geom_line() + theme_bw() +
    geom_line(data=tib_plot_sm, aes(x=k, y=log(value)), linetype="dashed", color="blue") +
    geom_line(data=tib_plot, aes(x=k, y=log(value)), linetype="dashed", color="red")

## ----weyl-ivols1---------------------------------------------------------
factorial(6) * weyl_ivols(5,"A")
factorial(6) * weyl_ivols(5,"A_red")   # intrinsic volumes of reduced cone are just shifted
2^5*factorial(5) * weyl_ivols(5,"BC")
2^4*factorial(5) * weyl_ivols(5,"D")

## ----weyl-ivols2---------------------------------------------------------
factorial(6)*2^5*factorial(5) * weyl_ivols( c(5,5) , c("A","BC"), product=TRUE )

v_list <- weyl_ivols( c(5,5) , c("A","BC") )
factorial(6)*2^5*factorial(5) * prod_ivols(v_list)

## ----circ-ivols----------------------------------------------------------
v <- circ_ivols(10,pi/5)
v
sum(v)
sum( v[ 2*(1:5) ] )       # odd index intrinsic volumes add up to 1/2

2 * v[ 2*(1:5) ]          # comparing with binomial distribution
dbinom(0:4,4,sin(pi/5)^2)

## ----biv-chi-bar-sq, fig.width = 5, fig.height=5-------------------------
# sample from the bivariate chi-bar-squared distribution of a product of circular cones
D <- c(7,17)
alpha <- c(0.7*pi/2, 0.6*pi/2)
v_true <- circ_ivols( D, alpha, product=TRUE)
m_samp <- rbichibarsq(1e5, v_true)
d <- sum(D)

# scatter plot of the sample
ggplot(as_tibble(m_samp), aes(V1,V2)) + geom_point(alpha=.02) +
    theme_bw() +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())

# estimate moments, compare with true values
est <- estim_statdim_var(d, m_samp); est
list( statdim_true=sum((0:d)*v_true),
      var_true=sum((0:d)^2*v_true)-sum((0:d)*v_true)^2 )

## ----statdim-var-circ, fig.width = 7-------------------------------------
d <- 9
N <- 1e3
alpha <- (0:N)/N * pi/2

Sdim <- matrix(0,d-1,N+1)
Var  <- matrix(0,d-1,N+1)
for (k in 2:d) {
    V <- circ_ivols( rep(k,N+1) , alpha)
    Sdim[k-1,] <- sapply(V, function(v) return(sum((0:k) * v)))
    Var[k-1,]  <- sapply(V, function(v) return(sum((0:k)^2 * v)))-Sdim[k-1,]^2
}

G <- ggplot()
for (k in 2:d) {
    G <- G + geom_line(data=tibble(sdim=Sdim[k-1,],var=Var[k-1,]), aes(sdim,var))
}
G <- G + theme_bw()
G

## ----statdim-var-part----------------------------------------------------
P <- parts(d); P

## ----statdim-var-prod-curv, fig.width = 7--------------------------------
sing <- which(colSums(P>1)==1)

G <- ggplot()
for (k in sing) {
    n_ray <- sum(P[,k]==1)
    sdim <- Sdim[ P[1,k]-1, ] + n_ray/2
    var  <- Var[  P[1,k]-1, ] + n_ray/4
    G <- G + geom_line(data=tibble(sdim,var), aes(sdim,var))
}
G <- G + theme_bw(); G

## ----statdim-var-prod-scatter, fig.width = 7, fig.height=10, echo=FALSE----
# line of d-dimensional circular cones
Gd <- geom_line(data=tibble(sdim=Sdim[d-1,],var=Var[d-1,]), aes(sdim,var))

# recalculate statdim and variance with less points
N_sm <- 20
alpha_sm <- (0:N_sm)/N_sm * pi/2
Sdim_sm <- matrix(0,d-1,N_sm+1)
Var_sm  <- matrix(0,d-1,N_sm+1)
for (k in 2:d) {
    V <- circ_ivols( rep(k,N_sm+1) , alpha_sm )
    Sdim_sm[k-1,] <- sapply(V, function(v) return(sum((0:k) * v)))
    Var_sm[k-1,]  <- sapply(V, function(v) return(sum((0:k)^2 * v)))-Sdim_sm[k-1,]^2
}

# indices of nontrivial cones
nonsing <- which(colSums(P>1)>1)

# function that produces scatterplots
plotSV <- function(k) {
    n_ray <- sum(P[,k]==1)
    n_full <- sum(P[,k]>1)
    sdim <- Sdim_sm[ P[1,k]-1, ]
    var  <- Var_sm[  P[1,k]-1, ]
    for (l in 2:n_full) {
        sdim <- c(rep(1,N_sm+1) %*% t(sdim) + Sdim_sm[ P[l,k]-1, ])
        var  <- c(rep(1,N_sm+1) %*% t(var)  + Var_sm[  P[l,k]-1, ])
    }
    sdim <- sdim + n_ray/2
    var  <- var  + n_ray/4
    return(geom_point(data=tibble(sdim,var), aes(sdim,var), size=1, alpha=0.05))
}

plotsSV <- list()
for (i in 1:length(nonsing)) {
    plotsSV[[i]] <- ggplot() + Gd + plotSV(nonsing[i]) + theme_bw() +
        theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
              axis.text.x =element_blank(), axis.text.y =element_blank(),
              axis.ticks.x=element_blank(), axis.ticks.y=element_blank())
}

multiplot(plotlist = plotsSV, cols = 3)

