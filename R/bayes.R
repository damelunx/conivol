

samp_post <- function(d, m_samp, N=20, v_init=NULL, init_mode=0, data=NULL, inform=TRUE, eps=1) {
    #################################
    # some general internal options:
    #
    #################################

    # find the values of the chi-squared densities at the sample points
    if (is.null(data))
        data <- conivol::prepare_data(d, m_samp)

    out_samples <- matrix(0,N+1,d+1)

    # set the starting point at initial guess for intrinsic volumes
    if (length(v_init)==d+1)
        v <- v_init
    else {
        est <- conivol::estimate_statdim_var(d, m_samp)
        v <- conivol::init_v(d,init_mode,delta=est$delta,var=est$var)
    }
    out_iterates[1, ] <- v

    for (i in 1:N) {
        # draw uniform samples for update
        u <- runif(d+1)
        # determine update order
        pi <- sample(0:d)
        for (k in pi) {
            e <- 1
            r_lowbd <- -inf
            r_upbd  <- inf
            while ( r_upbd-r_lowbd < eps ) {
                r_low <- asdf
                r_up  <- asdf
                integral <- asdf

                # https://www.r-bloggers.com/tailoring-univariate-probability-distributions/
                # uniroot

                e <- e/2
            }
            v[pi[k]]
        }
        out_iterates[i+1, ] <- v
    }

    return(out_samples)
}










