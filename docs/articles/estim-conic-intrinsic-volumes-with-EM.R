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
N <- 1e5
E <- ellips_rbichibarsq(N,A)
str(E)

# scatter plot of the sample
ggplot(as_tibble(E$samples), aes(V1,V2)) + geom_point(alpha=.02) +
    theme_bw() +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())
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
loglike_ivols( v_est$mode0 , out_prep )
loglike_ivols( v_est$mode1 , out_prep )
loglike_ivols( v_est$mode2 , out_prep )
loglike_ivols( v_est$mode3 , out_prep )
loglike_ivols( v_est$mode4 , out_prep )

## ----ellips-em, fig.width = 7--------------------------------------------
em0 <- estim_em( d, E$samples, N=200, v_init=v_est$mode0, data=out_prep )
em1 <- estim_em( d, E$samples, N=200, v_init=v_est$mode1, data=out_prep )
em2 <- estim_em( d, E$samples, N=200, v_init=v_est$mode2, data=out_prep )
em3 <- estim_em( d, E$samples, N=200, v_init=v_est$mode3, data=out_prep )
em4 <- estim_em( d, E$samples, N=200, v_init=v_est$mode4, data=out_prep )

# prepare plots
tib_plot0 <- as_tibble( t(em0$iterates[1+25*(0:8), ]) ) %>%
    `colnames<-`(paste0("s_",25*(0:8))) %>%
    add_column(k=0:d,.before=1) %>% gather(step,value,2:10)
tib_plot1 <- as_tibble( t(em1$iterates[1+25*(0:8), ]) ) %>%
    `colnames<-`(paste0("s_",25*(0:8))) %>%
    add_column(k=0:d,.before=1) %>% gather(step,value,2:10)
tib_plot2 <- as_tibble( t(em2$iterates[1+25*(0:8), ]) ) %>%
    `colnames<-`(paste0("s_",25*(0:8))) %>%
    add_column(k=0:d,.before=1) %>% gather(step,value,2:10)
tib_plot3 <- as_tibble( t(em3$iterates[1+25*(0:8), ]) ) %>%
    `colnames<-`(paste0("s_",25*(0:8))) %>%
    add_column(k=0:d,.before=1) %>% gather(step,value,2:10)
tib_plot4 <- as_tibble( t(em4$iterates[1+25*(0:8), ]) ) %>%
    `colnames<-`(paste0("s_",25*(0:8))) %>%
    add_column(k=0:d,.before=1) %>% gather(step,value,2:10)

tib_plot0$step <- factor(tib_plot0$step, levels = paste0("s_",25*(0:8)))
tib_plot1$step <- factor(tib_plot1$step, levels = paste0("s_",25*(0:8)))
tib_plot2$step <- factor(tib_plot2$step, levels = paste0("s_",25*(0:8)))
tib_plot3$step <- factor(tib_plot3$step, levels = paste0("s_",25*(0:8)))
tib_plot4$step <- factor(tib_plot4$step, levels = paste0("s_",25*(0:8)))

# plot mode 0
ggplot(tib_plot0,aes(x=k,y=value,color=step)) +
    geom_line() + theme_bw() +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank())
# plot mode 1
ggplot(tib_plot1,aes(x=k,y=value,color=step)) +
    geom_line() + theme_bw() +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank())
# plot mode 2
ggplot(tib_plot2,aes(x=k,y=value,color=step)) +
    geom_line() + theme_bw() +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank())
# plot mode 3
ggplot(tib_plot3,aes(x=k,y=value,color=step)) +
    geom_line() + theme_bw() +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank())
# plot mode 4
ggplot(tib_plot4,aes(x=k,y=value,color=step)) +
    geom_line() + theme_bw() +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank())

## ----ellips-em-comp, fig.width = 7---------------------------------------
# prepare plot
tib_plot_comp <- tibble(  mode0=em0$iterates[201, ],
                          mode1=em1$iterates[201, ],
                          mode2=em2$iterates[201, ],
                          mode3=em3$iterates[201, ],
                          mode4=em4$iterates[201, ] ) %>%
    add_column(k=0:d,.before=1) %>% gather(mode,value,2:6)

# plot 200th iterate of different runs
ggplot(tib_plot_comp,aes(x=k,y=value,color=mode)) +
    geom_line() + theme_bw() +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank())

## ----ellips-em-no-logconc, fig.width = 7---------------------------------
em_no_logconc <- estim_em( d, E$samples, N=200, v_init=v_est$mode0, data=out_prep,
                           no_of_lcc_projections=0 )

# prepare plot
tib_plot_no_logconc <- as_tibble( t(em_no_logconc$iterates[1+25*(0:8), ]) ) %>%
    `colnames<-`(paste0("s_",25*(0:8))) %>%
    add_column(k=0:d,.before=1) %>% gather(step,value,2:10)

tib_plot_no_logconc$step <- factor(tib_plot_no_logconc$step, levels = paste0("s_",25*(0:8)))

ggplot(tib_plot_no_logconc,aes(x=k,y=value,color=step)) +
    geom_line() + theme_bw() +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank())

