
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

You can install the released version of mcnuts from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("mcnuts")
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
#> theta[1] 46.1 49.4 52.7 49.4 2.0  1.00     6841     2237
#> theta[2] 40.6 43.9 47.2 43.9 2.0  1.00     6655     2333
#> theta[3] 50.7 53.9 57.2 53.9 2.0  1.00     6661     2586
#> theta[4] 54.4 57.7 61.1 57.7 2.0  1.00     6772     2236
#> theta[5] 36.1 39.1 42.3 39.2 1.9  1.00     6986     2291
#> theta[6] 51.7 55.1 58.3 55.0 2.0  1.00     7403     2031
#> theta[7] 75.5 78.7 82.1 78.7 2.0  1.00     6042     2149
#> theta[8] 68.8 72.0 75.1 72.0 1.9  1.00     5790     2482
#> mu       46.4 56.3 66.5 56.3 6.1  1.04      131      200
#> tau       9.9 15.5 27.5 16.8 5.9  1.01      214      220
#> 
#> For each parameter, Bulk_ESS and Tail_ESS are crude measures of 
#> effective sample size for bulk and tail quantities respectively (an ESS > 100 
#> per chain is considered good), and Rhat is the potential scale reduction 
#> factor on rank normalized split chains (at convergence, Rhat <= 1.05).
#> Avg trajectory length 0.73 0.79 0.77 0.79
```
