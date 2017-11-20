## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "conic-intrinsic-volumes_figures/"
)

## ----load-pkgs, include=FALSE--------------------------------------------
library(conivol)
library(tidyverse)

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

## ----sum-canc-data-------------------------------------------------------
dim(A)
A[1:5,1:5]
all(A>=0)
colSums(A)

## ----plot-canc-samp, fig.width = 7---------------------------------------
tib_plot <- as_tibble(msamp[1+linC:dimC]/n) %>%
    add_column( k=linC:dimC, .before=1)
ggplot(tib_plot, aes(x=k, y=value))      + geom_line() + theme_bw()
ggplot(tib_plot, aes(x=k, y=log(value))) + geom_line() + theme_bw()

## ----canc-bayes, fig.width = 7-------------------------------------------
bayes_est <- polyh_bayes( msamp, dimC, linC )
tib_plot <- bayes_est$post_samp(1e4) %>%
    as_tibble() %>%
    `colnames<-`(paste0(rep("V",dimC-linC+1),as.character(linC:dimC))) %>%
    gather(factor_key=TRUE)
ggplot(tib_plot, aes(x=key, y=value)) +
    geom_boxplot() + theme_bw() +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank())

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

