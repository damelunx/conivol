## ----results='hide', message=FALSE, warning=FALSE, echo=FALSE------------
### loading some packages
library(conivol)
library(tidyverse)
library(Rmosek)

## ------------------------------------------------------------------------
D <- c(5,8)
alpha <- c( asin(sqrt(0.9)) , asin(sqrt(0.8)))
v <- circ_ivol(D, alpha, product = TRUE)

## ------------------------------------------------------------------------
d <- sum(D)
delta <- sum(v*(0:d))
var <- sum(v*(0:d)^2) - delta^2
vn <- sapply( 0:d, 
     function(k){ pnorm((k+0.5-delta)/sqrt(var)) - pnorm((k-0.5-delta)/sqrt(var)) })
vn <- vn/sum(vn)

## ------------------------------------------------------------------------
search_alpha <- pi/2*(1:999)/1000
ivols <- circ_ivol(rep(d,length(search_alpha)),search_alpha)
i_fit <- which.min( ( sapply(ivols, function(v){sum(v*(0:d))}) - delta )^2 )
alpha_fit <- search_alpha[i_fit]
vc <- circ_ivol(d, alpha_fit)

## ---- echo=FALSE, fig.align='center', fig.width=7, eval=TRUE-------------
tib_plot     <- gather(tibble(k=0:d,exact=v,     est_norm=vn,     est_circ=vc),    type,value,2:4)
tib_plot_log <- gather(tibble(k=0:d,exact=log(v),est_norm=log(vn),est_circ=log(vc)),type,value,2:4)

ggplot(tib_plot, aes(x=k, y=value, color=type)) +
    geom_line() +
    scale_colour_manual(values=c("red","blue","black")) +
    theme(legend.position="bottom",
          axis.title.x=element_blank(),
          axis.title.y=element_blank())

ggplot(tib_plot_log, aes(x=k, y=value, color=type)) +
    geom_line() +
    scale_colour_manual(values=c("red","blue","black")) +
    theme(legend.position="bottom",
          axis.title.x=element_blank(),
          axis.title.y=element_blank())

## ---- fig.align='center', fig.width=7, fig.height=6, eval=TRUE-----------
n <- 10^5
set.seed(1234)
m_samp <- rbichibarsq_circ(n,D,alpha)
ggplot(as_tibble(m_samp), aes(V1,V2)) + geom_point(alpha=.02) +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())

## ---- eval=FALSE---------------------------------------------------------
#  EM_iterates_mode0 <- bichibarsq_find_weights( m_samp, d, N=1000, mode=0)
#  EM_iterates_mode1 <- bichibarsq_find_weights( m_samp, d, N=1000, mode=1)
#  EM_iterates_mode2 <- bichibarsq_find_weights( m_samp, d, N=1000, mode=2)

## ---- echo=FALSE, fig.align='center', fig.width=7, eval=TRUE-------------
tib_iter <- as_tibble( t(EM_iterates_mode0$iterates[1:8, ]) ) %>%
                add_column(k=0:d,.before=1)
tib_plot <- gather(tib_iter,step,value,2:dim(tib_iter)[2])
tib_exact <- tibble( k=0:d, exact=circ_ivol(D,alpha,TRUE))
ggplot(tib_plot,aes(x=k,y=value,color=step)) +
    geom_line() +
    geom_line(data=tib_exact,aes(x=k,y=exact),colour="black") + 
    theme(legend.position="none",
          axis.title.x=element_blank(),
          axis.title.y=element_blank()) +
    scale_colour_brewer()

tib_iter_log <- as_tibble( t(log(EM_iterates_mode0$iterates[1:8, ])) ) %>%
                add_column(k=0:d,.before=1)
tib_plot_log <- gather(tib_iter_log,step,value,2:dim(tib_iter_log)[2])
tib_exact_log <- tibble( k=0:d, exact=log(circ_ivol(D,alpha,TRUE)))
ggplot(tib_plot_log,aes(x=k,y=value,color=step)) +
    geom_line() +
    geom_line(data=tib_exact_log,aes(x=k,y=exact),colour="black") + 
    theme(legend.position="none",
          axis.title.x=element_blank(),
          axis.title.y=element_blank()) +
    scale_colour_brewer()

## ---- echo=FALSE, fig.align='center', fig.width=7, eval=TRUE-------------
tib_iter <- as_tibble( t(EM_iterates_mode1$iterates[1:8, ]) ) %>%
                add_column(k=0:d,.before=1)
tib_plot <- gather(tib_iter,step,value,2:dim(tib_iter)[2])
tib_exact <- tibble( k=0:d, exact=circ_ivol(D,alpha,TRUE))
ggplot(tib_plot,aes(x=k,y=value,color=step)) +
    geom_line() +
    geom_line(data=tib_exact,aes(x=k,y=exact),colour="black") + 
    theme(legend.position="none",
          axis.title.x=element_blank(),
          axis.title.y=element_blank()) +
    scale_colour_brewer()

tib_iter_log <- as_tibble( t(log(EM_iterates_mode1$iterates[1:8, ])) ) %>%
                add_column(k=0:d,.before=1)
tib_plot_log <- gather(tib_iter_log,step,value,2:dim(tib_iter_log)[2])
tib_exact_log <- tibble( k=0:d, exact=log(circ_ivol(D,alpha,TRUE)))
ggplot(tib_plot_log,aes(x=k,y=value,color=step)) +
    geom_line() +
    geom_line(data=tib_exact_log,aes(x=k,y=exact),colour="black") + 
    theme(legend.position="none",
          axis.title.x=element_blank(),
          axis.title.y=element_blank()) +
    scale_colour_brewer()

## ---- echo=FALSE, fig.align='center', fig.width=7, eval=TRUE-------------
tib_iter <- as_tibble( t(EM_iterates_mode2$iterates[1:8, ]) ) %>%
                add_column(k=0:d,.before=1)
tib_plot <- gather(tib_iter,step,value,2:dim(tib_iter)[2])
tib_exact <- tibble( k=0:d, exact=circ_ivol(D,alpha,TRUE))
ggplot(tib_plot,aes(x=k,y=value,color=step)) +
    geom_line() +
    geom_line(data=tib_exact,aes(x=k,y=exact),colour="black") + 
    theme(legend.position="none",
          axis.title.x=element_blank(),
          axis.title.y=element_blank()) +
    scale_colour_brewer()

tib_iter_log <- as_tibble( t(log(EM_iterates_mode2$iterates[1:8, ])) ) %>%
                add_column(k=0:d,.before=1)
tib_plot_log <- gather(tib_iter_log,step,value,2:dim(tib_iter_log)[2])
tib_exact_log <- tibble( k=0:d, exact=log(circ_ivol(D,alpha,TRUE)))
ggplot(tib_plot_log,aes(x=k,y=value,color=step)) +
    geom_line() +
    geom_line(data=tib_exact_log,aes(x=k,y=exact),colour="black") + 
    theme(legend.position="none",
          axis.title.x=element_blank(),
          axis.title.y=element_blank()) +
    scale_colour_brewer()

## ---- echo=FALSE, fig.align='center', fig.width=7, eval=TRUE-------------
tib_iter <- as_tibble( t(EM_iterates_mode0$iterates[c(10,20,100,500,1000,2000,4000,10000), ]) ) %>%
                add_column(k=0:d,.before=1)
tib_plot <- gather(tib_iter,step,value,2:dim(tib_iter)[2])
tib_exact <- tibble( k=0:d, exact=circ_ivol(D,alpha,TRUE))
ggplot(tib_plot,aes(x=k,y=value,color=step)) +
    geom_line() +
    geom_line(data=tib_exact,aes(x=k,y=exact),colour="black") + 
    theme(legend.position="none",
          axis.title.x=element_blank(),
          axis.title.y=element_blank()) +
    scale_colour_brewer()

tib_iter_log <- as_tibble( t(log(EM_iterates_mode0$iterates[c(10,20,100,500,1000,2000,4000,10000), ])) ) %>%
                add_column(k=0:d,.before=1)
tib_plot_log <- gather(tib_iter_log,step,value,2:dim(tib_iter_log)[2])
tib_exact_log <- tibble( k=0:d, exact=log(circ_ivol(D,alpha,TRUE)))
ggplot(tib_plot_log,aes(x=k,y=value,color=step)) +
    geom_line() +
    geom_line(data=tib_exact_log,aes(x=k,y=exact),colour="black") + 
    theme(legend.position="none",
          axis.title.x=element_blank(),
          axis.title.y=element_blank()) +
    scale_colour_brewer()

## ---- echo=FALSE, fig.align='center', fig.width=7, eval=TRUE-------------
tib_iter <- as_tibble( t(EM_iterates_mode1$iterates[c(10,20,100,500,1000,2000,4000,10000), ]) ) %>%
                add_column(k=0:d,.before=1)
tib_plot <- gather(tib_iter,step,value,2:dim(tib_iter)[2])
tib_exact <- tibble( k=0:d, exact=circ_ivol(D,alpha,TRUE))
ggplot(tib_plot,aes(x=k,y=value,color=step)) +
    geom_line() +
    geom_line(data=tib_exact,aes(x=k,y=exact),colour="black") + 
    theme(legend.position="none",
          axis.title.x=element_blank(),
          axis.title.y=element_blank()) +
    scale_colour_brewer()

tib_iter_log <- as_tibble( t(log(EM_iterates_mode1$iterates[c(10,20,100,500,1000,2000,4000,10000), ])) ) %>%
                add_column(k=0:d,.before=1)
tib_plot_log <- gather(tib_iter_log,step,value,2:dim(tib_iter_log)[2])
tib_exact_log <- tibble( k=0:d, exact=log(circ_ivol(D,alpha,TRUE)))
ggplot(tib_plot_log,aes(x=k,y=value,color=step)) +
    geom_line() +
    geom_line(data=tib_exact_log,aes(x=k,y=exact),colour="black") + 
    theme(legend.position="none",
          axis.title.x=element_blank(),
          axis.title.y=element_blank()) +
    scale_colour_brewer()

## ---- echo=FALSE, fig.align='center', fig.width=7, eval=TRUE-------------
tib_iter <- as_tibble( t(EM_iterates_mode2$iterates[c(10,20,100,500,1000,2000,4000,10000), ]) ) %>%
                add_column(k=0:d,.before=1)
tib_plot <- gather(tib_iter,step,value,2:dim(tib_iter)[2])
tib_exact <- tibble( k=0:d, exact=circ_ivol(D,alpha,TRUE))
ggplot(tib_plot,aes(x=k,y=value,color=step)) +
    geom_line() +
    geom_line(data=tib_exact,aes(x=k,y=exact),colour="black") + 
    theme(legend.position="none",
          axis.title.x=element_blank(),
          axis.title.y=element_blank()) +
    scale_colour_brewer()

tib_iter_log <- as_tibble( t(log(EM_iterates_mode2$iterates[c(10,20,100,500,1000,2000,4000,10000), ])) ) %>%
                add_column(k=0:d,.before=1)
tib_plot_log <- gather(tib_iter_log,step,value,2:dim(tib_iter_log)[2])
tib_exact_log <- tibble( k=0:d, exact=log(circ_ivol(D,alpha,TRUE)))
ggplot(tib_plot_log,aes(x=k,y=value,color=step)) +
    geom_line() +
    geom_line(data=tib_exact_log,aes(x=k,y=exact),colour="black") + 
    theme(legend.position="none",
          axis.title.x=element_blank(),
          axis.title.y=element_blank()) +
    scale_colour_brewer()

## ---- echo=FALSE, fig.align='center', fig.width=7, eval=TRUE-------------
#N_plot <- 100
tib_prog <- as_tibble( sweep( EM_iterates_mode0$iterates[101:1001, ], MARGIN=2, circ_ivol(D,alpha,TRUE), "/") ) %>%
                add_column(iterate=100:1000, .before=1)
tib_prog_plot <- gather(tib_prog, vol, value, 2:dim(tib_prog)[2])
ggplot(tib_prog_plot, aes(x=iterate, y=value, color=vol)) +
    geom_line() +
    scale_colour_brewer(palette = "Set3")

## ---- echo=FALSE, fig.align='center', fig.width=7, eval=TRUE-------------


## ---- echo=FALSE, fig.align='center', fig.width=7, eval=TRUE-------------


