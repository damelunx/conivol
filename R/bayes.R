#' JAGS model creation for Bayesian posterior given samples of bivariate chi-bar-squared distribution
#'
#' \code{find_ivols_jags} generates inputs for JAGS (data list and model string or external file)
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
#' @param filename filename for output (NA by default, in which case the return is a string)
#' @param overwrite logical; determines whether the output should overwrite an existing file
#'
#' @return If \code{filename==NA} then the output of \code{find_ivols_jags} is a list containing the following elements:
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
#' If \code{filename!=NA} then the model string will be written to the file with
#' the specified name and the output will only contain the elements \code{data}
#' and \code{variable.names}.
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
#' v <- circ_ivols(D, alpha, product = TRUE)
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
find_ivols_jags <- function(samples, d, dim=d, lin=0, prior="noninformative", v_prior=NA, filename=NA, overwrite=FALSE) {

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
    out$data           <- data_list
    out$variable.names <- "V"

    if ( !is.na(filename) && file.exists(filename) && overwrite==FALSE )
        stop("\n File with given filename exists and overwrite==FALSE.")
    else if ( !is.na(filename) ) {
        if (file.exists(filename))
            file.remove(filename)
        file_conn<-file(filename)
        writeLines(model_string, file_conn)
        close(file_conn)
    } else {
        out$model <- model_string
    }
    return(out)
}






#' Stan model creation for Bayesian posterior given samples of bivariate chi-bar-squared distribution
#'
#' \code{find_ivols_stan} generates inputs for Stan (data list and model string or external file)
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
#' @param filename filename for output (NA by default, in which case the return is a string)
#' @param overwrite logical; determines whether the output should overwrite an existing file
#'
#' @return If \code{filename==NA} then the output of \code{find_ivols_stan} is a list containing the following elements:
#' \itemize{
#'   \item \code{model}: a string that forms the description of the Stan model,
#'   \item \code{data}: a data list containing the prepared data to be used
#'                    for defining a Stan model object,
#'                    ## rest still has to be adapted ##
#'   \item \code{variable.names}: the single string "V" to be used as additional
#'                    parameter when creating samples from the JAGS object to
#'                    indicate that only this vector should be tracked.
#' }
#' If \code{filename!=NA} then the model string will be written to the file with
#' the specified name and the output will only contain the elements \code{data}
#' and \code{variable.names}.
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
#' v <- circ_ivols(D, alpha, product = TRUE)
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
find_ivols_stan <- function(samples, d, dim=d, lin=0, prior="noninformative", v_prior=NA, filename=NA, overwrite=FALSE) {

    I_pol  <- which(samples[ ,1]==0)
    I_prim <- which(samples[ ,2]==0)
    samples_bulk <- samples[ -c(I_pol,I_prim) , ]
    N_pol  <- length(I_pol)
    N_prim <- length(I_prim)
    N_bulk <- N-N_pol-N_prim

    I_0 <- 2*(0:floor((dim-lin)/2)) + 1       # final +1 is because of R indices start at 1
    I_1 <- 1+2*(0:floor((dim-lin-1)/2)) + 1   # final +1 is because of R indices start at 1

    if (is.na(v_prior))
        v_prior <- rep(1,dim-lin+1)/(dim-lin+1)

    if (prior=="informative")
        prior_weight <- 1
    else
        prior_weight <- 0

    if (dim==d && lin==0) {
        data_list <- list(
            d            = d ,
            N_0          = N_pol ,
            N_bulk       = N_bulk ,
            N_d          = N_prim ,
            X            = samples_bulk[ ,1] ,
            Y            = samples_bulk[ ,2] ,
            V_nonz_prior = v_prior ,
            prior_weight = prior_weight
        )
        model_string <- "
data {
    int <lower=1> d;                                         \\ ambient dimension
    int <lower=0> N_0;                                       \\ number of samples in polar cone
    int <lower=0> N_bulk;                                    \\ number of samples in neither cone
    int <lower=0> N_d;                                       \\ number of samples in primal cone
    vector <lower=0>[N_bulk] X;                              \\ squared length of projection on primal
    vector <lower=0>[N_bulk] Y;                              \\ squared length of projection on polar
    vector <lower=0>[d+1] V_nonz_prior ;                     \\ prior estimate of nonzero intrinsic volumes
    real <lower=0> prior_weight ;                            \\ added weight for Dirichlet priors (0=noninformative, 1=informative)
}
transformed data {
    int sample_multinom[3] ;                                 \\ collecting numbers of points in primal/polar cone
    sample_mutinom[1] = N_0 ;
    sample_mutinom[2] = N_bulk ;
    sample_mutinom[3] = N_d ;

    int d_0 ;                                                \\ dimension of even simplex
    int d_1 ;                                                \\ dimension of odd simplex
    d_0 =  ceil((d-1)/2) ;
    d_1 = floor((d-1)/2) ;

    int I_0[d_0+1] ;                                         \\ set of even indices of nonzero entries
    int I_1[d_1+1] ;                                         \\ set of odd indices of nonzero entries
    I_0 = 2*(0:floor(d/2)) + 1                               \\ final +1 is because indices start at 1
    I_1 = 1+2*(0:floor((d-1)/2)) + 1                         \\ final +1 is because indices start at 1

    vector [d_0+1] alpha ;                                   \\ weight for even Dirichlet prior
    vector [d_1+1] beta ;                                    \\ weight for odd Dirichlet prior
    alpha = V_nonz_prior[I_0] + prior_weight ;
    beta  = V_nonz_prior[I_1] + prior_weight ;

    row_vector[d-1] log_dens_XY[N_bulk] ;
    for (i in 1:N_bulk) {
        for (k in 1:(d-1)) {
            log_dens_XY[i][k] = chi_square_lpdf(X[i]|k) + chi_square_lpdf(Y[i]|(d-k)) ;
        }
    }
}
parameters {
    simplex[d_0+1] V_0 ;
    simplex[d_1+1] V_1 ;
}
transformed parameters {
    simplex[d+1] V ;
    V[I_0] = V_0/2 ;
    V[I_1] = V_1/2 ;

    simplex[3] V_extreme ;
    V_extreme[1] = V[1] ;
    V_extreme[2] = sum(V[2:d]) ;
    V_extreme[3] = V[d+1] ;

    row_vector[d-1] logV_bulk ;
    logV_bulk = to_row_vector( log(V[2:d]) ) ;
}
model {
    V_0 ~ dirichlet(alpha)
    V_1 ~ dirichlet(beta)

    sample_multinom ~ multinomial(V_extreme)

    \\ the following is the main part of the modelling where we marginalized the discrete latent variable K
    for (i in 1:N_bulk) {
        target += log_sum_exp( logV_bulk + log_dens_XY[i] ) ;
    }
}"
    } else if (dim==d && lin>0) {
        data_list <- list(
            d            = d ,
            lin          = lin ,
            N_bulk       = N_bulk ,
            N_d          = N_prim ,
            X            = samples_bulk[ ,1] ,
            Y            = samples_bulk[ ,2] ,
            V_nonz_prior = v_prior ,
            prior_weight = prior_weight
        )
        model_string <- "
data {
    int <lower=1> d;                                         \\ ambient dimension
    int <lower=1> lin;                                       \\ lineality
    int <lower=0> N_bulk;                                    \\ number of samples in neither cone
    int <lower=0> N_d;                                       \\ number of samples in primal cone
    vector <lower=0>[N_bulk] X;                              \\ squared length of projection on primal
    vector <lower=0>[N_bulk] Y;                              \\ squared length of projection on polar
    vector <lower=0>[d-lin+1] V_nonz_prior ;                 \\ prior estimate of nonzero intrinsic volumes
    real <lower=0> prior_weight ;                            \\ added weight for Dirichlet priors (0=noninformative, 1=informative)
}
transformed data {
    int sample_multinom[2] ;                                 \\ collecting numbers of points in primal/polar cone
    sample_mutinom[1] = N_bulk ;
    sample_mutinom[2] = N_d ;

    int d_0 ;                                                \\ dimension of even simplex
    int d_1 ;                                                \\ dimension of odd simplex
    d_0 =  ceil((d-lin-1)/2) ;
    d_1 = floor((d-lin-1)/2) ;

    int I_0[d_0+1] ;                                         \\ set of even indices of nonzero entries
    int I_1[d_1+1] ;                                         \\ set of odd indices of nonzero entries
    I_0 = 2*(0:floor((d-lin)/2)) + 1                         \\ final +1 is because indices start at 1
    I_1 = 1+2*(0:floor((d-lin-1)/2)) + 1                     \\ final +1 is because indices start at 1

    vector [d_0+1] alpha ;                                   \\ weight for even Dirichlet prior
    vector [d_1+1] beta ;                                    \\ weight for odd Dirichlet prior
    alpha = V_nonz_prior[I_0] + prior_weight ;
    beta  = V_nonz_prior[I_1] + prior_weight ;

    row_vector[d-lin] log_dens_XY[N_bulk] ;
    for (i in 1:N_bulk) {
        for (k in lin:(d-1)) {
            log_dens_XY[i][k-lin+1] = chi_square_lpdf(X[i]|k) + chi_square_lpdf(Y[i]|(d-k)) ;
        }
    }
}
parameters {
    simplex[d_0+1] V_0 ;
    simplex[d_1+1] V_1 ;
}
transformed parameters {
    simplex[d+1] V ;
    V[1:lin] = 0 ;
    V[lin+I_0] = V_0/2 ;
    V[lin+I_1] = V_1/2 ;

    simplex[2] V_extreme ;
    V_extreme[1] = sum(V[(lin+1):d]) ;
    V_extreme[2] = V[d+1] ;

    row_vector[d-lin] logV_bulk ;
    logV_bulk = to_row_vector( log(V[(lin+1):d]) ) ;
}
model {
    V_0 ~ dirichlet(alpha)
    V_1 ~ dirichlet(beta)

    sample_multinom ~ multinomial(V_extreme)

    \\ the following is the main part of the modelling where we marginalized the discrete latent variable K
    for (i in 1:N_bulk) {
        target += log_sum_exp( logV_bulk + log_dens_XY[i] ) ;
    }
}"
    } else if (dim<d && lin==0) {
        data_list <- list(
            d            = d ,
            dim          = dim ,
            N_0          = N_pol ,
            N_bulk       = N_bulk ,
            X            = samples_bulk[ ,1] ,
            Y            = samples_bulk[ ,2] ,
            V_nonz_prior = v_prior ,
            prior_weight = prior_weight
        )
        model_string <- "
data {
    int <lower=1> d;                                         \\ ambient dimension
    int <lower=1, upper=d-1> dim;                            \\ dimension of linear span
    int <lower=0> N_0;                                       \\ number of samples in polar cone
    int <lower=0> N_bulk;                                    \\ number of samples in neither cone
    vector <lower=0>[N_bulk] X;                              \\ squared length of projection on primal
    vector <lower=0>[N_bulk] Y;                              \\ squared length of projection on polar
    vector <lower=0>[dim+1] V_nonz_prior ;                   \\ prior estimate of nonzero intrinsic volumes
    real <lower=0> prior_weight ;                            \\ added weight for Dirichlet priors (0=noninformative, 1=informative)
}
transformed data {
    int sample_multinom[2] ;                                 \\ collecting numbers of points in primal/polar cone
    sample_mutinom[1] = N_0 ;
    sample_mutinom[2] = N_bulk ;

    int d_0 ;                                                \\ dimension of even simplex
    int d_1 ;                                                \\ dimension of odd simplex
    d_0 =  ceil((dim-1)/2) ;
    d_1 = floor((dim-1)/2) ;

    int I_0[d_0+1] ;                                         \\ set of even indices of nonzero entries
    int I_1[d_1+1] ;                                         \\ set of odd indices of nonzero entries
    I_0 = 2*(0:floor(dim/2)) + 1                             \\ final +1 is because indices start at 1
    I_1 = 1+2*(0:floor((dim-1)/2)) + 1                       \\ final +1 is because indices start at 1

    vector [d_0+1] alpha ;                                   \\ weight for even Dirichlet prior
    vector [d_1+1] beta ;                                    \\ weight for odd Dirichlet prior
    alpha = V_nonz_prior[I_0] + prior_weight ;
    beta  = V_nonz_prior[I_1] + prior_weight ;

    row_vector[dim] log_dens_XY[N_bulk] ;
    for (i in 1:N_bulk) {
        for (k in 1:dim) {
            log_dens_XY[i][k] = chi_square_lpdf(X[i]|k) + chi_square_lpdf(Y[i]|(d-k)) ;
        }
    }
}
parameters {
    simplex[d_0+1] V_0 ;
    simplex[d_1+1] V_1 ;
}
transformed parameters {
    simplex[d+1] V ;
    V[(dim+2):(dim+1)] = 0 ;
    V[I_0] = V_0/2 ;
    V[I_1] = V_1/2 ;

    simplex[2] V_extreme ;
    V_extreme[1] = V[1] ;
    V_extreme[2] = sum(V[2:(dim+1)]) ;

    row_vector[dim] logV_bulk ;
    logV_bulk = to_row_vector( log(V[2:(dim+1)]) ) ;
}
model {
    V_0 ~ dirichlet(alpha)
    V_1 ~ dirichlet(beta)

    sample_multinom ~ multinomial(V_extreme)

    \\ the following is the main part of the modelling where we marginalized the discrete latent variable K
    for (i in 1:N_bulk) {
        target += log_sum_exp( logV_bulk + log_dens_XY[i] ) ;
    }
}"
    } else if (dim<d && lin>0) {
        data_list <- list(
            d            = d ,
            dim          = dim ,
            lin          = lin ,
            N            = N_bulk ,
            X            = samples_bulk[ ,1] ,
            Y            = samples_bulk[ ,2] ,
            V_nonz_prior = v_prior ,
            prior_weight = prior_weight
        )
        model_string <- "
data {
    int <lower=1> d;                                         \\ ambient dimension
    int <lower=1, upper=d-1> dim;                            \\ dimension of linear span
    int <lower=1, upper=dim-1> lin;                          \\ lineality
    int <lower=0> N;                                         \\ number of samples
    vector <lower=0>[N] X;                                   \\ squared length of projection on primal
    vector <lower=0>[N] Y;                                   \\ squared length of projection on polar
    vector <lower=0>[dim-lin+1] V_nonz_prior ;               \\ prior estimate of nonzero intrinsic volumes
    real <lower=0> prior_weight ;                            \\ added weight for Dirichlet priors (0=noninformative, 1=informative)
}
transformed data {
    int d_0 ;                                                \\ dimension of even simplex
    int d_1 ;                                                \\ dimension of odd simplex
    d_0 =  ceil((dim-lin-1)/2) ;
    d_1 = floor((dim-lin-1)/2) ;

    int I_0[d_0+1] ;                                         \\ set of even indices of nonzero entries
    int I_1[d_1+1] ;                                         \\ set of odd indices of nonzero entries
    I_0 = 2*(0:floor((dim-lin)/2)) + 1                       \\ final +1 is because indices start at 1
    I_1 = 1+2*(0:floor((dim-lin-1)/2)) + 1                   \\ final +1 is because indices start at 1

    vector [d_0+1] alpha ;                                   \\ weight for even Dirichlet prior
    vector [d_1+1] beta ;                                    \\ weight for odd Dirichlet prior
    alpha = V_nonz_prior[I_0] + prior_weight ;
    beta  = V_nonz_prior[I_1] + prior_weight ;

    row_vector[dim-lin+1] log_dens_XY[N] ;
    for (i in 1:N) {
        for (k in lin:dim) {
            log_dens_XY[i][k-lin+1] = chi_square_lpdf(X[i]|k) + chi_square_lpdf(Y[i]|(d-k)) ;
        }
    }
}
parameters {
    simplex[d_0+1] V_0 ;
    simplex[d_1+1] V_1 ;
}
transformed parameters {
    simplex[dim-lin+1] V ;
    V[1:lin] = 0 ;
    V[(dim+2):(d+1)] = 0 ;
    V[lin+I_0] = V_0/2 ;
    V[lin+I_1] = V_1/2 ;

    row_vector[dim-lin+1] logV_bulk ;
    logV_bulk = to_row_vector( log(V[(lin+1):(dim+1)]) ) ;
}
model {
    V_0 ~ dirichlet(alpha)
    V_1 ~ dirichlet(beta)

    \\ the following is the main part of the modelling where we marginalized the discrete latent variable K
    for (i in 1:N) {
        target += log_sum_exp( logV_bulk + log_dens_XY[i] ) ;
    }
}"
    }

    out                <- list()
    out$data           <- data_list
    out$variable.names <- "V"

    if ( !is.na(filename) && file.exists(filename) && overwrite==FALSE )
        stop("\n File with given filename exists and overwrite==FALSE.")
    else if ( !is.na(filename) ) {
        if (file.exists(filename))
            file.remove(filename)
        file_conn<-file(filename)
        writeLines(model_string, file_conn)
        close(file_conn)
    } else {
        out$model <- model_string
    }
    return(out)
}


