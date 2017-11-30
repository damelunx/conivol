## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "conic-intrinsic-volumes_figures/"
)

## ----load-pkgs, include=FALSE--------------------------------------------
library(conivol)
library(tidyverse)
library(knitr)
library(png)
library(rstan)

## ----load-imgs, include=FALSE--------------------------------------------
img_paths <- list( dl="bayes_diagrams/bayes_direct_logconc.png",
                   dnl="bayes_diagrams/bayes_direct_nologconc.png",
                   idle="bayes_diagrams/bayes_indirect_logconc_even.png",
                   idlo="bayes_diagrams/bayes_indirect_logconc_odd.png",
                   idnle="bayes_diagrams/bayes_indirect_nologconc_even.png",
                   idnlo="bayes_diagrams/bayes_indirect_nologconc_odd.png" )
img_dpi <- 450

## ----define-cone---------------------------------------------------------
d <- 17
num_gen <- 50
set.seed(1324)
A <- matrix( c(rep(1,num_gen), rnorm((d-1)*num_gen)), d, num_gen, byrow = TRUE )
out <- polyh_reduce_gen(A)
str(out)
dimC <- out$dimC
linC <- out$linC
A_red <- out$A_reduced

## ----coll-samples--------------------------------------------------------
samp_iv_sm <- polyh_rivols_gen(1e2, A_red, reduce=FALSE)$multsamp
samp_iv_la <- polyh_rivols_gen(1e4, A_red, reduce=FALSE)$multsamp
samp_bcb_sm <- polyh_rbichibarsq_gen(1e2, A_red, reduce=FALSE)
samp_bcb_la <- polyh_rbichibarsq_gen(1e4, A_red, reduce=FALSE)

## ----start-V-------------------------------------------------------------
v0_iv_sm <- init_ivols( dimC, sum(0:dimC * samp_iv_sm/sum(samp_iv_sm)), init_mode=1 )
v0_iv_la <- init_ivols( dimC, sum(0:dimC * samp_iv_la/sum(samp_iv_la)), init_mode=1 )
v0_bcb_sm <- init_ivols( dimC, estim_statdim_var(dimC,samp_bcb_sm)$delta, init_mode=1)
v0_bcb_la <- init_ivols( dimC, estim_statdim_var(dimC,samp_bcb_la)$delta, init_mode=1)

## ----start-V-disp, fig.width = 7, echo=FALSE-----------------------------
tib_plot <- tibble(iv_small = v0_iv_sm,
                   iv_large = v0_iv_la,
                   bcb_small = v0_bcb_sm,
                   bcb_large = v0_bcb_la) %>%
    add_column(k=0:dimC,.before=1) %>% gather(sampling,v,2:5,factor_key=TRUE)
ggplot(tib_plot,aes(x=k,y=v,color=sampling)) +
    geom_line() + theme_bw() +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank())

## ----prior-iv-noninf, fig.width=7, echo=FALSE----------------------------
# we use the built-in function; by providing zero weights the posterior will equal the prior
prior_iv <- polyh_bayes(rep(0,dimC+1),dimC,linC,v_prior=v0_iv_sm)
N <- 10
tib_plot <- t(prior_iv$post_samp(N)) %>%
    as_tibble() %>% add_column(k=linC:dimC, .before=1) %>%
    gather(key,value,2:(N+1))
ggplot(tib_plot, aes(x=k, y=(value), color=key)) +
    geom_line() + theme_bw() +
    theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank())

## ----prior-iv-inf, fig.width=7, echo=FALSE-------------------------------
# we use the built-in function; by providing zero weights the posterior will equal the prior
prior_iv <- polyh_bayes(rep(0,dimC+1),dimC,linC,v_prior=v0_iv_sm,prior_sample_size=d)
N <- 10
tib_plot <- t(prior_iv$post_samp(N)) %>%
    as_tibble() %>% add_column(k=linC:dimC, .before=1) %>%
    gather(key,value,2:(N+1))
ggplot(tib_plot, aes(x=k, y=(value), color=key)) +
    geom_line() + theme_bw() +
    theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank())

## ----prior-iv-noninf-logc, fig.width=7, echo=FALSE-----------------------
# construct matrix T
T <- matrix(0,d+1,d+1)
for (k in 0:(d-2)) {
    T[k+1,k+1] <- 1
    T[k+1,k+2] <- -2
    T[k+1,k+3] <- 1
}
T[d,] <- rep_len(c(1,0),d+1)
T[d+1,] <- rep_len(c(0,1),d+1)
if (d%%2==0) T[d+1,1] <- 1

# find s^(0) and sample from prior
s0 <- T %*% (-log(v0_iv_sm))

S_samp <- matrix(0,d+1,N)
for (k in 0:d)
    S_samp[k+1,] <- rgamma(N,shape=2,scale=s0[k+1]/2)

# compute back to sample for V
V_samp_unnorm <- t(exp(-solve(T,S_samp)))
V_samp <- t(V_samp_unnorm/rowSums(V_samp_unnorm))

tib_plot <- V_samp %>%
    as_tibble() %>% add_column(k=linC:dimC, .before=1) %>%
    gather(key,value,2:(N+1))
ggplot(tib_plot, aes(x=k, y=(value), color=key)) +
    geom_line() + theme_bw() +
    theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank())

## ----prior-iv-inf-logc, fig.width=7, echo=FALSE--------------------------
S_samp <- matrix(0,d+1,N)
for (k in 0:d)
    S_samp[k+1,] <- rgamma(N,shape=s0[k+1],scale=1)

# compute back to sample for V
V_samp_unnorm <- t(exp(-solve(T,S_samp)))
V_samp <- t(V_samp_unnorm/rowSums(V_samp_unnorm))

tib_plot <- V_samp %>%
    as_tibble() %>% add_column(k=linC:dimC, .before=1) %>%
    gather(key,value,2:(N+1))
ggplot(tib_plot, aes(x=k, y=(value), color=key)) +
    geom_line() + theme_bw() +
    theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank())

## ----disp-dnl, echo=FALSE------------------------------------------------
include_graphics(img_paths$dnl, dpi=img_dpi)

## ----dir-enf-par-sm-noninf-comp------------------------------------------
post_iv <- polyh_bayes(samp_iv_sm,dimC,linC,v_prior=v0_iv_sm,prior_sample_size=1)
str(post_iv)

## ----dir-enf-par-sm-noninf, fig.width=7, echo=FALSE----------------------
tib_plot <- t(post_iv$post_samp(N)) %>%
    as_tibble() %>% add_column(k=linC:dimC, .before=1) %>%
    gather(key,value,2:(N+1))
ggplot(tib_plot, aes(x=k, y=(value), color=key)) +
    geom_line() + theme_bw() +
    theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank())

## ----dir-enf-par-sm-inf, fig.width=7, echo=FALSE-------------------------
post_iv <- polyh_bayes(samp_iv_sm,dimC,linC,v_prior=v0_iv_sm,prior_sample_size=d)
tib_plot <- t(post_iv$post_samp(N)) %>%
    as_tibble() %>% add_column(k=linC:dimC, .before=1) %>%
    gather(key,value,2:(N+1))
ggplot(tib_plot, aes(x=k, y=(value), color=key)) +
    geom_line() + theme_bw() +
    theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank())

## ----dir-enf-par-la-noninf, fig.width=7, echo=FALSE----------------------
post_iv <- polyh_bayes(samp_iv_la,dimC,linC,v_prior=v0_iv_la)
tib_plot <- t(post_iv$post_samp(N)) %>%
    as_tibble() %>% add_column(k=linC:dimC, .before=1) %>%
    gather(key,value,2:(N+1))
ggplot(tib_plot, aes(x=k, y=(value), color=key)) +
    geom_line() + theme_bw() +
    theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank())

## ----dir-enf-par-la-inf, fig.width=7, echo=FALSE-------------------------
post_iv <- polyh_bayes(samp_iv_la,dimC,linC,v_prior=v0_iv_la,prior_sample_size=d)
tib_plot <- t(post_iv$post_samp(N)) %>%
    as_tibble() %>% add_column(k=linC:dimC, .before=1) %>%
    gather(key,value,2:(N+1))
ggplot(tib_plot, aes(x=k, y=(value), color=key)) +
    geom_line() + theme_bw() +
    theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank())

## ----disp-dl, echo=FALSE-------------------------------------------------
include_graphics(img_paths$dl, dpi=img_dpi)

## ----dir-enf-logc--------------------------------------------------------
# asdf

## ----disp-idnle, echo=FALSE----------------------------------------------
include_graphics(img_paths$idnle, dpi=img_dpi)

## ----disp-idnlo, echo=FALSE----------------------------------------------
include_graphics(img_paths$idnlo, dpi=img_dpi)

## ----indir-enf-par-------------------------------------------------------
# asdf

## ----disp-idle, echo=FALSE-----------------------------------------------
include_graphics(img_paths$idle, dpi=img_dpi)

## ----disp-idlo, echo=FALSE-----------------------------------------------
include_graphics(img_paths$idlo, dpi=img_dpi)

## ----indir-enf-logc------------------------------------------------------
# asdf

