#' EM iterates for product of circular cones example.
#'
#' The iterates obtained from the EM algorithm to find the weights of the
#' bivariate chi-bar-squared distribution.
#'
#' @format A matrix of fourteen columns, first row gives the starting point.
#'
#' @source
#' \code{n <- 10^5}
#'
#' \code{D <- c(5,8)}
#'
#' \code{alpha <- c( asin(sqrt(0.9)) , asin(sqrt(0.8)))}
#'
#' \code{set.seed(1234)}
#'
#' \code{m_samp <- rbichibarsq_circ(n,D,alpha)}
#'
#' \code{EM_iterates_mode0 <- bichibarsq_find_weights( m_samp, d, N=1000, mode=0)}
#'
#' (takes about ten minutes)
#'
"EM_iterates_mode0"


#' EM iterates for product of circular cones example.
#'
#' The iterates obtained from the EM algorithm to find the weights of the
#' bivariate chi-bar-squared distribution.
#'
#' @format A matrix of fourteen columns, first row gives the starting point.
#'
#' @source
#' \code{n <- 10^5}
#'
#' \code{D <- c(5,8)}
#'
#' \code{alpha <- c( asin(sqrt(0.9)) , asin(sqrt(0.8)))}
#'
#' \code{set.seed(1234)}
#'
#' \code{m_samp <- rbichibarsq_circ(n,D,alpha)}
#'
#' \code{EM_iterates_mode0 <- bichibarsq_find_weights( m_samp, d, N=1000, mode=1)}
#'
#' (takes about ten minutes)
#'
"EM_iterates_mode1"


#' EM iterates for product of circular cones example.
#'
#' The iterates obtained from the EM algorithm to find the weights of the
#' bivariate chi-bar-squared distribution.
#'
#' @format A matrix of fourteen columns, first row gives the starting point.
#'
#' @source
#' \code{n <- 10^5}
#'
#' \code{D <- c(5,8)}
#'
#' \code{alpha <- c( asin(sqrt(0.9)) , asin(sqrt(0.8)))}
#'
#' \code{set.seed(1234)}
#'
#' \code{m_samp <- rbichibarsq_circ(n,D,alpha)}
#'
#' \code{EM_iterates_mode0 <- bichibarsq_find_weights( m_samp, d, N=1000, mode=2)}
#'
#' (takes about ten minutes)
#'
"EM_iterates_mode2"
