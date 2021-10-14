#' Run one Hamiltonian Monte Carlo iteration
#'
#' Code based on Gelman, A., Carlin, J.B., Stern, H.S., Dunson, D.B., Vehtari, A., &
#' Rubin, D.B. (2013). Bayesian Data Analysis (3rd ed.). Chapman and Hall/CRC
#'
#' @param th a vector of parameters of interests
#' @param log_p_th a function that takes a vector of parameters and returns a log likelihood
#' @param grad_th a function that takes a vector of parameters and returns the gradient
#' @param M covariance for the momentum
#' @param L an integer, trajectory length
#' @param epsilon step size
#' @return
#'      th: new sample of the parameters of interest
#'      p_jump: current iteration's acceptance probability
#' @export
hmc_iteration <- function(th, log_p_th, grad_th, epsilon, L, M) {
  M_inv <- 1 / M # M is vector
  d <- length(th)
  phi <- rnorm(d, 0, sqrt(M))
  th_old <- th
  log_p_old <- log_p_th(th_old) - 0.5 * sum(M_inv * phi^2)
  phi <- phi + 0.5 * epsilon * gradient_th(th_old)
  for (i in 1:L) {
    th <- th + epsilon * M_inv * phi
    phi <- phi + epsilon * gradient_th(th)
  }
  phi <- -phi
  log_p_star <- log_p_th(th) - 0.5 * sum(M_inv * phi^2)
  r <- exp(log_p_star - log_p_old)
  if (is.na(r)) r <- 0
  p_jump <- min(r, 1)
  th_new <- if (runif(1) < p_jump) th else th_old
  return(list(th = th_new, p_jump = p_jump))
}

#' Run many chains of a Hamiltonian Monte Carlo sampler
#'
#' @param log_p_th a function that takes a vector of parameters and returns a log likelihood
#' @param grad_th a function that takes a vector of parameters and returns the gradient
#' @param starting_values c x d matrix of starting values, c: number of chains, d: number of parameters
#' @param iter number of MCMC iterations
#' @param M covariance for the momentum
#' @param L_0 an integer, initial trajectory length
#' @param epsilon_0 initial step size
#' @return
#'      th: new sample of the parameters of interest
#'      p_jump: current iteration's acceptance probability
#' @export
hmc_run <- function(log_p_th, grad_th, starting_values, iter, epsilon_0, L_0, M) {
  chains <- nrow(starting_values)
  d <- ncol(starting_values)
  sims <- array(NA,
    dim = c(iter, chains, d),
    dimnames = list(NULL, NULL, colnames(starting_values))
  )
  warmup <- 0.5 * iter
  p_jump <- array(NA, c(iter, chains))
  for (j in 1:chains) {
    th <- starting_values[j, ]
    for (t in 1:iter) {
      epsilon <- runif(1, 0, 2 * epsilon_0)
      L <- ceiling(2 * L_0 * runif(1))
      temp <- hmc_iteration(th, log_p_th, grad_th, epsilon, L, M)
      p_jump[t, j] <- temp$p_jump
      sims[t, j, ] <- temp$th
      th <- temp$th
    }
  }
  rstan::monitor(sims, warmup)
  cat(
    "Avg acceptance probs:",
    round(colMeans(p_jump[(warmup + 1):iter, ]), 2), "\n"
  )
  return(list(sims = sims, p_jump = p_jump))
}
