#' JAGS model creation for Bayesian posterior given samples of bivariate chi-bar-squared distribution
#'
#' \code{find_ivols_jags} generates inputs for JAGS (data list and model string)
#' for sampling from the posterior distribution,
#' given samples of the bivariate chi-bar-squared distribution.
#'
#' @param samples N-by-2 matrix representing independent samples from the
#'                bivariate chi-bar-squared distribution of a convex cone
#' @param d the dimension of the ambient space
#' @param dim the dimension of the cone
#' @param lin the lineality of the cone
#' @param prior either "noninformative" (default) or "informative"
#' @param v_prior a prior estimate of the vector of intrinsic volumes (NA by default)
#'
#' @return The output of \code{find_ivols_jags} is a list containing the following elements:
#' \itemize{
#'   \item \code{model}: a string that forms the description of the JAGS model,
#'                    which can be directly used as input (via text connection
#'                    or external file) for creating a JAGS model object;
#'                    \code{post_samp(n)} returns an \code{n}-by-\code{(dim+1)}
#'                    matrix whose rows form a set of \code{n} independent samples
#'                    of the posterior distribution,
#'   \item \code{data}: a data list containing the prepared data to be used
#'                    with the model string for defining a JAGS model object,
#'   \item \code{variable.names}: the single string "V" to be used as additional
#'                    parameter when creating samples from the JAGS object to
#'                    indicate that only this vector should be tracked.
#' }
#'
#'
#'
#' @note See the example below for details on how to us the outputs with rjags
#'       functions; see \href{../doc/bayesian.html}{this vignette}
#'       for further info.
#'
#' @section See also:
#' \code{\link[conivol]{find_ivols_EM}}
#'
#' Package: \code{\link[conivol]{conivol}}
#'
#' @examples
#'
#' library(rjags)
#' options(mc.cores = parallel::detectCores())
#'
#' # defining the example cone
#' D <- c(5,8)
#' alpha <- c( asin(sqrt(0.9)) , asin(sqrt(0.8)))
#' v <- circ_ivol(D, alpha, product = TRUE)
#' d <- sum(D)
#'
#' # getting the sample data
#' N <- 10^3
#' set.seed(1234)
#' samples <- circ_rbichibarsq(N,D,alpha)
#'
#' # compute initial guess
#' est <- estimate_statdim_var(d, m_samp)
#' v0 <- init_v(d,init_mode=1,delta=est$delta,var=est$var)
#'
#' # obtain input data for JAGS model; use v0 as prior
#' in_jags <- find_ivols_jags(samples, d, prior="informative", v_prior=v0)
#'
#' # create JAGS model
#' model_connection <- textConnection(in_jags$model)
#' mod <- jags.model(model_connection ,
#'                   data = in_jags$data ,
#'                   n.chains = 4 ,
#'                   n.adapt = 500)
#' close(model_connection)
#' update(mod, 1e3)
#'
#' # simulate posterior distribution and display trace plots and summaries
#' mod_sim <- coda.samples(model=mod, variable.names=in_jags$variable.names, n.iter=1e4)\
#' plot(mod_sim, ask=TRUE)
#' mod_csim <- as.mcmc(do.call(rbind, mod_sim))
#' tmp <- summary(mod_csim)
#'
#' # plot true values, the estimate v0, and marginals of the posterior samples
#' library(tidyverse)
#' ggplot(data=tibble(x=0:d,
#'                    y_true=v,
#'                    y_est=v0,
#'                    y_bayes_mean=tmp$statistics[ ,'Mean'],
#'                    y_bayes_quant25=tmp$quantiles[ ,'25%'],
#'                    y_bayes_quant50=tmp$quantiles[ ,'50%'],
#'                    y_bayes_quant75=tmp$quantiles[ ,'75%'])) +
#'      geom_line(aes(x=x, y=y_true)) +
#'      geom_line(aes(x=x, y=y_est), color="red") +
#'      geom_line(aes(x=x, y=y_bayes_mean), color="blue") +
#'      geom_line(aes(x=x, y=y_bayes_quant25), color="green") +
#'      geom_line(aes(x=x, y=y_bayes_quant50), color="green") +
#'      geom_line(aes(x=x, y=y_bayes_quant75), color="green")
#'
#' @export
#'
find_ivols_jags <- function(samples, d, dim=d, lin=0, prior="noninformative", v_prior=NA) {

    I_pol  <- which(samples[ ,1]==0)
    I_prim <- which(samples[ ,2]==0)
    samples_bulk <- samples[ -c(I_pol,I_prim) , ]
    N_pol  <- length(I_pol)
    N_prim <- length(I_prim)
    N_bulk <- N-N_pol-N_prim

    I_0 <- 2*(0:floor((dim-lin)/2)) + 1       # final +1 is because of R indices start at 1
    I_1 <- 1+2*(0:floor((dim-lin-1)/2)) + 1   # final +1 is because of R indices start at 1

    if (is.na(v_prior)) {
        v_prior_adj      <- rep(0,dim-lin+1)
        v_prior_adj[I_0] <- 1/ceiling((dim-lin+1)/2) / 2
        v_prior_adj[I_1] <- 1/floor((dim-lin+1)/2) / 2
    } else {
        v_prior_adj <- v_prior
    }

    Dir_prior_0 <- 2*v_prior_adj[I_0]
    Dir_prior_1 <- 2*v_prior_adj[I_1]

    if (prior=="informative"){
        Dir_prior_0 <- 1+Dir_prior_0
        Dir_prior_1 <- 1+Dir_prior_1
    }

    data_list <- list(
        d           = d ,
        dim         = dim ,
        lin         = lin ,
        d_0         = floor((dim-lin)/2) ,
        d_1         = ceiling((dim-lin)/2)-1 ,
        N           = N ,
        N_bulk      = N_bulk ,
        count_bulk  = c(N_V0, N_bulk, N_Vd) ,
        X           = samples_bulk[ ,1] ,
        Y           = samples_bulk[ ,2] ,
        Dir_prior_0 = Dir_prior_0 ,
        Dir_prior_1 = Dir_prior_1
    )

    if (dim==d && lin==0) {
        if (d%%2==1) {
            model_string <- "model {
                V_0 ~ ddirch(Dir_prior_0)
                V_1 ~ ddirch(Dir_prior_1)
                V_extreme[1] <- V_0[1]
                V_extreme[2] <- sum(V_0)+sum(V_1)-V_0[1]-V_1[d_1+1]
                V_extreme[3] <- V_1[d_1+1]
                count_bulk ~ dmulti(V_extreme, N)
                V[1]   <- V_extreme[1]
                V[d+1] <- V_extreme[3]
                for (i in 1:d_1) {
                    V_bulk[2*i-1] <- V_1[i]
                    V[2*i]        <- V_1[i]
                }
                for (i in 2:(d_0+1)) {
                    V_bulk[2*(i-1)] <- V_0[i]
                    V[2*i-1]        <- V_0[i]
                }
                for (i in 1:N_bulk) {
                    K[i] ~ dcat(V_bulk)
                    X[i] ~ dchisqr(K[i])
                    Y[i] ~ dchisqr(d-K[i])
                }
            }"
        } else {
            model_string <- "model {
                V_0 ~ ddirch(Dir_prior_0)
                V_1 ~ ddirch(Dir_prior_1)
                V_extreme[1] <- V_0[1]
                V_extreme[2] <- sum(V_0)+sum(V_1)-V_0[1]-V_0[d_0+1]
                V_extreme[3] <- V_0[d_0+1]
                count_bulk ~ dmulti(V_extreme, N)
                V[1]   <- V_extreme[1]
                V[d+1] <- V_extreme[3]
                for (i in 1:(d_1+1)) {
                    V_bulk[2*i-1] <- V_1[i]
                    V[2*i]        <- V_1[i]
                }
                for (i in 2:d_0) {
                    V_bulk[2*(i-1)] <- V_0[i]
                    V[2*i-1]        <- V_0[i]
                }
                for (i in 1:N_bulk) {
                    K[i] ~ dcat(V_bulk)
                    X[i] ~ dchisqr(K[i])
                    Y[i] ~ dchisqr(d-K[i])
                }
            }"
        }
    } else if (dim==d && lin>0) {
        if ((d-lin)%%2==1) {
            model_string <- "model {
                V_0 ~ ddirch(Dir_prior_0)
                V_1 ~ ddirch(Dir_prior_1)
                V_extreme[1] <- sum(V_0)+sum(V_1)-V_1[d_1+1]
                V_extreme[2] <- V_1[d_1+1]
                count_bulk ~ dmulti(V_extreme, N)
                for (i in 1:lin) {
                    V[i] <- 0
                }
                V[d+1] <- V_extreme[2]
                for (i in 1:(d_0+1)) {
                    V_bulk[2*i-1] <- V_0[i]
                    V[lin+2*i-1]  <- V_0[i]
                }
                for (i in 1:d_1) {
                    V_bulk[2*i] <- V_1[i]
                    V[lin+2*i]  <- V_1[i]
                }
                for (i in 1:N_bulk) {
                    K[i] ~ dcat(V_bulk)
                    X[i] ~ dchisqr(lin+K[i])
                    Y[i] ~ dchisqr(d-K[i]-lin)
                }
            }"
        } else {
            model_string <- "model {
                V_0 ~ ddirch(Dir_prior_0)
                V_1 ~ ddirch(Dir_prior_1)
                V_extreme[1] <- sum(V_0)+sum(V_1)-V_0[d_0+1]
                V_extreme[2] <- V_0[d_0+1]
                count_bulk ~ dmulti(V_extreme, N)
                for (i in 1:lin) {
                    V[i] <- 0
                }
                V[d+1] <- V_extreme[2]
                for (i in 1:d_0) {
                    V_bulk[2*i-1] <- V_0[i]
                    V[lin+2*i-1]  <- V_0[i]
                }
                for (i in 1:(d_1+1)) {
                    V_bulk[2*i] <- V_1[i]
                    V[lin+2*i]  <- V_1[i]
                }
                for (i in 1:N_bulk) {
                    K[i] ~ dcat(V_bulk)
                    X[i] ~ dchisqr(lin+K[i])
                    Y[i] ~ dchisqr(d-K[i]-lin)
                }
            }"
        }
    } else if (dim<d && lin==0) {
        model_string <- "model {
            V_0 ~ ddirch(Dir_prior_0)
            V_1 ~ ddirch(Dir_prior_1)
            V_extreme[1] <- V_0[1]
            V_extreme[2] <- sum(V_0)+sum(V_1)-V_0[1]
            count_bulk ~ dmulti(V_extreme, N)
            V[1] <- V_extreme[1]
            for (i in (dim+2):(d+1)) {
                V[i] <- 0
            }
            for (i in 1:(d_1+1)) {
                V_bulk[2*i-1] <- V_1[i]
                V[2*i]        <- V_1[i]
            }
            for (i in 2:(d_0+1)) {
                V_bulk[2*(i-1)] <- V_0[i]
                V[2*i-1]        <- V_0[i]
            }
            for (i in 1:N_bulk) {
                K[i] ~ dcat(V_bulk)
                X[i] ~ dchisqr(K[i])
                Y[i] ~ dchisqr(d-K[i])
            }
        }"
    } else if (dim<d && lin>0) {
        model_string <- "model {
            V_0 ~ ddirch(Dir_prior_0)
            V_1 ~ ddirch(Dir_prior_1)
            for (i in 1:lin) {
                V[i] <- 0
            }
            for (i in (dim+2):(d+1)) {
                V[i] <- 0
            }
            for (i in 1:(d_0+1)) {
                V_bulk[2*i-1] <- V_0[i]
                V[lin+2*i-1]  <- V_0[i]
            }
            for (i in 1:(d_1+1)) {
                V_bulk[2*i] <- V_1[i]
                V[lin+2*i]  <- V_1[i]
            }
            for (i in 1:N_bulk) {
                K[i] ~ dcat(V_bulk)
                X[i] ~ dchisqr(lin+K[i])
                Y[i] ~ dchisqr(d-K[i]-lin)
            }
        }"
    }

    out                <- list()
    out$model          <- model_string
    out$data           <- data_list
    out$variable.names <- "V"
    return(out)
}



