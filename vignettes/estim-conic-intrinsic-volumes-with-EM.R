## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "estim-conic-intrinsic-volumes-with-EM_figures/"
)

## ----load-pkgs, include=FALSE--------------------------------------------
library(conivol)
library(tidyverse)

## ----ellips-def----------------------------------------------------------
d <-12
set.seed(1234)
A <- matrix(rnorm(d^2),d,d)
ellips_semiax(A)

## ----ellips-samp, fig.width = 7------------------------------------------
N <- 1e4
E <- ellips_rbichibarsq(N,A)
str(E)
mom <- estim_statdim_var(d, E$samples); mom
v_est <- list( mode0 = init_ivols(d, 0),
               mode1 = init_ivols(d, 1, mom$delta, mom$var),
               mode2 = init_ivols(d, 2, mom$delta),
               mode3 = init_ivols(d, 3, mom$delta, mom$var),
               mode4 = init_ivols(d, 4, mom$delta, mom$var) )
tib_plot <- as_tibble(v_est) %>%
    add_column(k=0:d,.before=1) %>% gather(mode,v,2:6)
ggplot(tib_plot,aes(x=k,y=v,color=mode)) +
    geom_line() + theme_bw()

## ----ellips-prep-data----------------------------------------------------
out_prep <- prepare_em(d, E$samples)          

## ----ellips-loglike------------------------------------------------------
#asdf

