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
library(rjags)
options(mc.cores = parallel::detectCores())

## ----load-imgs, include=FALSE--------------------------------------------
img_paths <- list( d="bayes_diagrams/bayes_direct.png",
                   dp="bayes_diagrams/bayes_direct_par.png",
                   dl="bayes_diagrams/bayes_direct_logconc.png",
                   id="bayes_diagrams/bayes_indirect.png",
                   idpe="bayes_diagrams/bayes_indirect_par_even.png",
                   idpo="bayes_diagrams/bayes_indirect_par_odd.png",
                   idle="bayes_diagrams/bayes_indirect_logconc_even.png",
                   idlo="bayes_diagrams/bayes_indirect_logconc_odd.png" )
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
samp_iv_la <- polyh_rivols_gen(2e3, A_red, reduce=FALSE)$multsamp
samp_bcb_sm <- polyh_rbichibarsq_gen(1e2, A_red, reduce=FALSE)
samp_bcb_la <- polyh_rbichibarsq_gen(2e3, A_red, reduce=FALSE)

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

## ----disp-d, echo=FALSE--------------------------------------------------
include_graphics(img_paths$d, dpi=img_dpi)

## ----disp-dp, echo=FALSE-------------------------------------------------
include_graphics(img_paths$dp, dpi=img_dpi)

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

## ----dir-enf-logc-comp, warning=FALSE------------------------------------
filename <- "ex_stan_model_iv.stan"
staninp <- polyh_stan(samp_iv_sm, dimC, linC,
                      v_prior=v0_iv_sm, prior="noninformative", filename=filename)

## ----include=FALSE-------------------------------------------------------
# this is so that some compiler warnings do not appear in the markdown
# (those warnings are only shown when running the Stan model the first time)
stanfit <- stan( file = filename, data = staninp$data, chains = 1,
                 warmup = 1, iter = 2, cores = 1, refresh = 10 )

## ----dir-enf-logc-comp2, fig.width=7, warning=FALSE----------------------
stanfit <- stan( file = filename, data = staninp$data, chains = 1,
                 warmup = 1000, iter = 2000, cores = 2, refresh = 1000 )

str(rstan::extract(stanfit))

post_iv_logc <- rstan::extract(stanfit)$V[100*(1:10),]

## ----dir-enf-logc-sm-noninf, fig.width=7, echo=FALSE---------------------
tib_plot <- t(post_iv_logc) %>%
    as_tibble() %>% add_column(k=linC:dimC, .before=1) %>%
    gather(key,value,2:(N+1))
ggplot(tib_plot, aes(x=k, y=(value), color=key)) +
    geom_line() + theme_bw() +
    theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank())

## ----dir-enf-logc-sm-inf, fig.width=7, echo=FALSE, warning=FALSE---------
# obtain Stan input (model already defined)
staninp <- polyh_stan(samp_iv_sm, dimC, linC,
                      v_prior=v0_iv_sm, prior="informative")
# run Stan model
stanfit <- stan( file = filename, data = staninp$data, chains = 1,
                 warmup = 1000, iter = 2000, cores = 2, refresh = 1000 )

post_iv_logc <- rstan::extract(stanfit)$V[100*(1:10),]

tib_plot <- t(post_iv_logc) %>%
    as_tibble() %>% add_column(k=linC:dimC, .before=1) %>%
    gather(key,value,2:(N+1))
ggplot(tib_plot, aes(x=k, y=(value), color=key)) +
    geom_line() + theme_bw() +
    theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank())

## ----dir-enf-logc-la-noninf, fig.width=7, echo=FALSE, warning=FALSE------
# obtain Stan input (model already defined)
staninp <- polyh_stan(samp_iv_la, dimC, linC,
                      v_prior=v0_iv_la, prior="noninformative")
# run Stan model
stanfit <- stan( file = filename, data = staninp$data, chains = 1,
                 warmup = 1000, iter = 2000, cores = 2, refresh = 1000 )

post_iv_logc <- rstan::extract(stanfit)$V[100*(1:10),]

tib_plot <- t(post_iv_logc) %>%
    as_tibble() %>% add_column(k=linC:dimC, .before=1) %>%
    gather(key,value,2:(N+1))
ggplot(tib_plot, aes(x=k, y=(value), color=key)) +
    geom_line() + theme_bw() +
    theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank())

## ----dir-enf-logc-la-inf, fig.width=7, echo=FALSE, warning=FALSE---------
# obtain Stan input (model already defined)
staninp <- polyh_stan(samp_iv_la, dimC, linC,
                      v_prior=v0_iv_la, prior="informative")
# run Stan model
stanfit <- stan( file = filename, data = staninp$data, chains = 1,
                 warmup = 1000, iter = 2000, cores = 2, refresh = 1000 )
# remove Stan input file
file.remove(filename)

post_iv_logc <- rstan::extract(stanfit)$V[100*(1:10),]

tib_plot <- t(post_iv_logc) %>%
    as_tibble() %>% add_column(k=linC:dimC, .before=1) %>%
    gather(key,value,2:(N+1))
ggplot(tib_plot, aes(x=k, y=(value), color=key)) +
    geom_line() + theme_bw() +
    theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank())

## ----disp-id, echo=FALSE-------------------------------------------------
include_graphics(img_paths$id, dpi=img_dpi)

## ----disp-idpe, echo=FALSE-----------------------------------------------
include_graphics(img_paths$idpe, dpi=img_dpi)

## ----disp-idpo, echo=FALSE-----------------------------------------------
include_graphics(img_paths$idpo, dpi=img_dpi)

## ----indir-JAGS-enf-par-sm-ninf, warning=FALSE, fig.width=7--------------
# obtain input data for JAGS model
in_jags <- estim_jags(samp_bcb_sm, d, v_prior=v0_bcb_sm)

# create JAGS model
model_connection <- textConnection(in_jags$model)
mod <- jags.model(model_connection ,
                  data = in_jags$data ,
                  n.chains = 1 ,
                  n.adapt = 1000)
close(model_connection)

# sample posterior distribution, take every thousandth sample
mod_sim <- coda.samples(model=mod, variable.names=in_jags$variable.names, n.iter=1e4)
mod_csim <- as.mcmc(do.call(rbind, mod_sim))
post_bcb_pa <- mod_csim[1e3*(1:10),]

tib_plot <- t(post_bcb_pa) %>%
    as_tibble() %>% add_column(k=linC:dimC, .before=1) %>%
    gather(key,value,2:(N+1))
ggplot(tib_plot, aes(x=k, y=(value), color=key)) +
    geom_line() + theme_bw() +
    theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank())

## ----indir-JAGS-enf-par-sm-inf, warning=FALSE, fig.width=7, echo=FALSE, warning=FALSE, results=FALSE----
# obtain input data for JAGS model
in_jags <- estim_jags(samp_bcb_sm, d, v_prior=v0_bcb_sm, prior_sample_size=d)

# create JAGS model
model_connection <- textConnection(in_jags$model)
mod <- jags.model(model_connection ,
                  data = in_jags$data ,
                  n.chains = 1 ,
                  n.adapt = 1000)
close(model_connection)

# sample posterior distribution, take every thousandth sample
mod_sim <- coda.samples(model=mod, variable.names=in_jags$variable.names, n.iter=1e4)
mod_csim <- as.mcmc(do.call(rbind, mod_sim))
post_bcb_pa <- mod_csim[1e3*(1:10),]

tib_plot <- t(post_bcb_pa) %>%
    as_tibble() %>% add_column(k=linC:dimC, .before=1) %>%
    gather(key,value,2:(N+1))
ggplot(tib_plot, aes(x=k, y=(value), color=key)) +
    geom_line() + theme_bw() +
    theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank())

## ----indir-JAGS-enf-par-la-ninf, warning=FALSE, fig.width=7, echo=FALSE, warning=FALSE, results=FALSE----
# obtain input data for JAGS model
in_jags <- estim_jags(samp_bcb_la, d, v_prior=v0_bcb_la)

# create JAGS model
model_connection <- textConnection(in_jags$model)
mod <- jags.model(model_connection ,
                  data = in_jags$data ,
                  n.chains = 1 ,
                  n.adapt = 1000)
close(model_connection)

# sample posterior distribution, take every thousandth sample
mod_sim <- coda.samples(model=mod, variable.names=in_jags$variable.names, n.iter=1e4)
mod_csim <- as.mcmc(do.call(rbind, mod_sim))
post_bcb_pa <- mod_csim[1000*(1:10),]

tib_plot <- t(post_bcb_pa) %>%
    as_tibble() %>% add_column(k=linC:dimC, .before=1) %>%
    gather(key,value,2:(N+1))
ggplot(tib_plot, aes(x=k, y=(value), color=key)) +
    geom_line() + theme_bw() +
    theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank())

## ----indir-JAGS-enf-par-la-inf, warning=FALSE, fig.width=7, echo=FALSE, warning=FALSE, results=FALSE----
# obtain input data for JAGS model
in_jags <- estim_jags(samp_bcb_la, d, v_prior=v0_bcb_la, prior_sample_size=d)

# create JAGS model
model_connection <- textConnection(in_jags$model)
mod <- jags.model(model_connection ,
                  data = in_jags$data ,
                  n.chains = 1 ,
                  n.adapt = 1000)
close(model_connection)

# sample posterior distribution, take every thousandth sample
mod_sim <- coda.samples(model=mod, variable.names=in_jags$variable.names, n.iter=1e4)
mod_csim <- as.mcmc(do.call(rbind, mod_sim))
post_bcb_pa <- mod_csim[1000*(1:10),]

tib_plot <- t(post_bcb_pa) %>%
    as_tibble() %>% add_column(k=linC:dimC, .before=1) %>%
    gather(key,value,2:(N+1))
ggplot(tib_plot, aes(x=k, y=(value), color=key)) +
    geom_line() + theme_bw() +
    theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank())

## ----indir-Stan-enf-par-sm-ninf-comp1, warning=FALSE, fig.width=7--------
filename <- "ex_stan_model_bcb.stan"
staninp <- estim_stan(samp_bcb_sm, d, dimC, linC,
                      v_prior=v0_bcb_sm, filename=filename)

## ----include=FALSE-------------------------------------------------------
# this is so that some compiler warnings do not appear in the markdown
# (those warnings are only shown when running the Stan model the first time)
stanfit <- stan( file = filename, data = staninp$data, chains = 1,
                 warmup = 1, iter = 2, cores = 1, refresh = 10 )

## ----indir-Stan-enf-par-sm-ninf-comp2, warning=FALSE, fig.width=7--------
stanfit <- stan( file = filename, data = staninp$data, chains = 1,
                 warmup = 1000, iter = 2000, cores = 2, refresh = 1000 )

str(rstan::extract(stanfit))

post_bcb_pa <- rstan::extract(stanfit)$V[100*(1:10),]

## ----indir-Stan-enf-par-sm-ninf, fig.width=7, echo=FALSE-----------------
tib_plot <- t(post_bcb_pa) %>%
    as_tibble() %>% add_column(k=linC:dimC, .before=1) %>%
    gather(key,value,2:(N+1))
ggplot(tib_plot, aes(x=k, y=(value), color=key)) +
    geom_line() + theme_bw() +
    theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank())

## ----indir-Stan-enf-par-sm-inf, warning=FALSE, fig.width=7, echo=FALSE----
# obtain Stan input (model already defined)
staninp <- estim_stan(samp_bcb_sm, d, dimC, linC,prior_sample_size=d,
                      v_prior=v0_bcb_sm)
stanfit <- stan( file = filename, data = staninp$data, chains = 1,
                 warmup = 1000, iter = 2000, cores = 2, refresh = 1000 )

post_bcb_pa <- rstan::extract(stanfit)$V[100*(1:10),]

tib_plot <- t(post_bcb_pa) %>%
    as_tibble() %>% add_column(k=linC:dimC, .before=1) %>%
    gather(key,value,2:(N+1))
ggplot(tib_plot, aes(x=k, y=(value), color=key)) +
    geom_line() + theme_bw() +
    theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank())

## ----indir-Stan-enf-par-la-ninf, warning=FALSE, fig.width=7, echo=FALSE----
# obtain Stan input (model already defined)
staninp <- estim_stan(samp_bcb_la, d, dimC, linC,
                      v_prior=v0_bcb_la)
stanfit <- stan( file = filename, data = staninp$data, chains = 1,
                 warmup = 1000, iter = 2000, cores = 2, refresh = 1000 )

post_bcb_pa <- rstan::extract(stanfit)$V[100*(1:10),]

tib_plot <- t(post_bcb_pa) %>%
    as_tibble() %>% add_column(k=linC:dimC, .before=1) %>%
    gather(key,value,2:(N+1))
ggplot(tib_plot, aes(x=k, y=(value), color=key)) +
    geom_line() + theme_bw() +
    theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank())

## ----indir-Stan-enf-par-la-inf, warning=FALSE, fig.width=7, echo=FALSE----
# obtain Stan input (model already defined)
staninp <- estim_stan(samp_bcb_la, d, dimC, linC,prior_sample_size=d,
                      v_prior=v0_bcb_la)
stanfit <- stan( file = filename, data = staninp$data, chains = 1,
                 warmup = 1000, iter = 2000, cores = 2, refresh = 1000 )
# remove Stan input file
file.remove(filename)

post_bcb_pa <- rstan::extract(stanfit)$V[100*(1:10),]

tib_plot <- t(post_bcb_pa) %>%
    as_tibble() %>% add_column(k=linC:dimC, .before=1) %>%
    gather(key,value,2:(N+1))
ggplot(tib_plot, aes(x=k, y=(value), color=key)) +
    geom_line() + theme_bw() +
    theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank())

## ----disp-idle, echo=FALSE-----------------------------------------------
include_graphics(img_paths$idle, dpi=img_dpi)

## ----disp-idlo, echo=FALSE-----------------------------------------------
include_graphics(img_paths$idlo, dpi=img_dpi)

## ----indir-Stan-enf-logc-sm-ninf-comp1, warning=FALSE, fig.width=7-------
filename <- "ex_stan_model_bcb_logc.stan"
staninp <- estim_stan(samp_bcb_sm, d, dimC, linC, enforce_logconc=TRUE,
                      v_prior=v0_bcb_sm, filename=filename)

## ----include=FALSE-------------------------------------------------------
# this is so that some compiler warnings do not appear in the markdown
# (those warnings are only shown when running the Stan model the first time)
stanfit <- stan( file = filename, data = staninp$data, chains = 1,
                 warmup = 1, iter = 2, cores = 1, refresh = 10 )

## ----indir-Stan-enf-logc-sm-ninf-comp2, warning=FALSE, fig.width=7-------
stanfit <- stan( file = filename, data = staninp$data, chains = 1,
                 warmup = 1000, iter = 2000, cores = 2, refresh = 1000 )

str(rstan::extract(stanfit))

post_bcb_logc <- rstan::extract(stanfit)$V[100*(1:10),]

## ----indir-Stan-enf-logc-sm-ninf, fig.width=7, echo=FALSE----------------
tib_plot <- t(post_bcb_logc) %>%
    as_tibble() %>% add_column(k=linC:dimC, .before=1) %>%
    gather(key,value,2:(N+1))
ggplot(tib_plot, aes(x=k, y=(value), color=key)) +
    geom_line() + theme_bw() +
    theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank())

## ----indir-Stan-enf-logc-sm-inf, warning=FALSE, fig.width=7, echo=FALSE----
# obtain Stan input (model already defined)
staninp <- estim_stan(samp_bcb_sm, d, dimC, linC, enforce_logconc=TRUE,
                      prior="informative", v_prior=v0_bcb_sm)
stanfit <- stan( file = filename, data = staninp$data, chains = 1,
                 warmup = 1000, iter = 2000, cores = 2, refresh = 1000 )

post_bcb_logc <- rstan::extract(stanfit)$V[100*(1:10),]

tib_plot <- t(post_bcb_logc) %>%
    as_tibble() %>% add_column(k=linC:dimC, .before=1) %>%
    gather(key,value,2:(N+1))
ggplot(tib_plot, aes(x=k, y=(value), color=key)) +
    geom_line() + theme_bw() +
    theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank())

## ----indir-Stan-enf-logc-la-ninf, warning=FALSE, fig.width=7, echo=FALSE----
# obtain Stan input (model already defined)
staninp <- estim_stan(samp_bcb_la, d, dimC, linC, enforce_logconc=TRUE,
                      v_prior=v0_bcb_la)
stanfit <- stan( file = filename, data = staninp$data, chains = 1,
                 warmup = 1000, iter = 2000, cores = 2, refresh = 1000 )

post_bcb_logc <- rstan::extract(stanfit)$V[100*(1:10),]

tib_plot <- t(post_bcb_logc) %>%
    as_tibble() %>% add_column(k=linC:dimC, .before=1) %>%
    gather(key,value,2:(N+1))
ggplot(tib_plot, aes(x=k, y=(value), color=key)) +
    geom_line() + theme_bw() +
    theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank())

## ----indir-Stan-enf-logc-la-inf, warning=FALSE, fig.width=7, echo=FALSE----
# obtain Stan input (model already defined)
staninp <- estim_stan(samp_bcb_la, d, dimC, linC, enforce_logconc=TRUE,
                      prior="informative", v_prior=v0_bcb_la)
stanfit <- stan( file = filename, data = staninp$data, chains = 1,
                 warmup = 1000, iter = 2000, cores = 2, refresh = 1000 )
# remove Stan input file
file.remove(filename)

post_bcb_logc <- rstan::extract(stanfit)$V[100*(1:10),]

tib_plot <- t(post_bcb_logc) %>%
    as_tibble() %>% add_column(k=linC:dimC, .before=1) %>%
    gather(key,value,2:(N+1))
ggplot(tib_plot, aes(x=k, y=(value), color=key)) +
    geom_line() + theme_bw() +
    theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank())

