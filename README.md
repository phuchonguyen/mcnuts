
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mcnuts

<!-- badges: start -->

<!-- badges: end -->

The goal of mcnuts is to provide implementation of one iteration of
Hamiltonian Monte Carlo and the No-U-Turn Sampler for easy integration
into a custom MCMC algorithm that might, for example, include Gibbs
steps.

**References**:

  - Betancourt, M. (2018). A Conceptual Introduction to Hamiltonian
    Monte Carlo.

  - Hoffman, M.D., & Gelman, A. (2014). The No-U-Turn Sampler:
    Adaptively Setting Path Lengths in Hamiltonian Monte Carlo. Journal
    of Machine Learning Research.

  - Gelman, A., Carlin, J.B., Stern, H.S., Dunson, D.B., Vehtari, A., &
    Rubin, D.B. (2013). Bayesian Data Analysis (3rd ed.). Chapman and
    Hall/CRC. <https://doi.org/10.1201/b16018>

## Installation

You can install the working development version through GitHub

``` r
devtools::install_github("phuchonguyen/mcnuts")
```

## Example

This is a basic example of how to fit a hierarchical mean model using
NUTS. The example is taken from Gelman et. al. 2013.

``` r
library(mcnuts)
```

### Make synthetic data

``` r
J = 8
sigma = 2
mu = 50
tau = 15
alpha = rnorm(J, mu, tau)
y = MASS::mvrnorm(100, alpha, diag(sigma^2, J))
y = colMeans(y)
```

### Define the log posterior density and gradient for the model

``` r
# Log posterior density
log_p_th = function(th) {
  J = length(th)-2
  theta = th[1:J]
  mu = th[J+1]
  tau = th[J+2]
  if (is.na(tau) | tau<=0)  # returns zero if tau is negative
    return(-Inf)
  else{
    log_hyperprior = 1
    log_prior = sum(dnorm(theta, mu, tau, log=TRUE))
    log_likelihood = sum(dnorm(y, theta, sigma, log=TRUE))
    return(sum(log_hyperprior, log_prior, log_likelihood))
  }
}
# Gradient th
gradient_th = function(th){
  J = length(th) - 2
  theta = th[1:J]
  mu = th[J+1]
  tau = th[J+2]
  if (tau <= 0 | is.na(tau))
    return(rep(0, length(th)))
  else{
    d_theta = -(theta - y)/sigma^2 - (theta-mu)/tau^2
    d_mu = -sum(mu-theta)/tau^2
    d_tau = -J/tau + sum((mu-theta)^2)/tau^3
    return(c(d_theta, d_mu, d_tau))
  }
}
```

### Run NUTS

``` r
parameter_names = c(paste0("theta[",1:8,"]"), "mu", "tau")
d = length(parameter_names)
chains = 4
mass_vector = c(rep(1/15^2, d-2), 1/5, 1/10)
starts = array(NA, c(chains, d), dimnames = list(NULL, parameter_names))
for (j in 1:chains){
  starts[j,] = rnorm(d, 0, 15)
  starts[j,10] = runif(1, 0, 15)
}
out = nuts_run_parallel(log_p_th, gradient_th, starting_values = starts, iter = 2000, M = mass_vector)
#> 
#> Inference for the input samples (4 chains: each with iter = 2000; warmup = 1000):
#> 
#>            Q5  Q50  Q95 Mean  SD  Rhat Bulk_ESS Tail_ESS
#> theta[1] 66.8 69.9 73.0 69.9 1.9  1.00     4062     2485
#> theta[2] 38.0 41.2 44.4 41.2 2.0  1.00     5098     2325
#> theta[3] 39.3 42.5 45.7 42.5 2.0  1.00     7423     2821
#> theta[4] 39.5 42.6 45.9 42.6 1.9  1.00     5754     2606
#> theta[5] 47.6 50.7 53.9 50.7 1.9  1.00     7387     2844
#> theta[6] 57.2 60.5 63.6 60.4 1.9  1.01     6301     2482
#> theta[7] 50.7 53.9 57.2 53.9 2.0  1.00     5494     2617
#> theta[8] 56.1 59.4 62.6 59.4 2.0  1.00     5755     2523
#> mu       44.3 53.0 62.3 53.1 5.4  1.01      216      290
#> tau       7.6 12.4 22.7 13.7 5.7  1.04      152      121
#> 
#> For each parameter, Bulk_ESS and Tail_ESS are crude measures of 
#> effective sample size for bulk and tail quantities respectively (an ESS > 100 
#> per chain is considered good), and Rhat is the potential scale reduction 
#> factor on rank normalized split chains (at convergence, Rhat <= 1.05).
#> Avg trajectory length 0.71 0.78 0.67 0.67
```
