#' mcnuts: A package with some gradient-based MCMC algorithms
#'
#' This package implements Hamiltonion Monte Carlo and the No-U-Turn Sampler
#' for easy integration into a custom MCMC algorithm that might, for example,
#' include Gibbs steps.
#'
#' @importFrom stats rnorm runif
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom foreach foreach %dopar%
#' @docType package
#' @name mcnuts
NULL
