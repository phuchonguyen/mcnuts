#' Run NUTS with many chains in parallel using the foreach and doParallel packages
#'
#' @param log_p_th a function that takes a vector of parameters and returns a log likelihood
#' @param grad_th a function that takes a vector of parameters and returns the gradient
#' @param starting_values c x d matrix of starting values, c: number of chains, d: number of parameters
#' @param iter number of MCMC iterations
#' @param M covariance for the momentum
#' @param p_adapt the fragtion of iterations where step size epsilon is adapted
#' @param epsilon_bar initial guess for epsilon
#' @param verbose boolean indicating to print progress or not
#' @return
#'      sims: iter x c x d array of MCMC sample
#'      L: vector of trajectory length
#'      eps: vector of step size
#' @export
nuts_run_parallel <- function(log_p_th, grad_th, starting_values, iter, M = NULL,
                              p_adapt = 0.5, epsilon_bar = 1, H = 0, gamma = 0.05, t0 = 10, kappa = 0.75, delta = 0.65,
                              verbose = FALSE) {
  # Set up parallel backend
  cores <- parallel::detectCores()
  cl <- parallel::makeCluster(cores[1] - 1) # not to overload your computer
  doParallel::registerDoParallel(cl)

  chains <- nrow(starting_values)
  d <- ncol(starting_values)
  if (is.null(M)) M <- rep(1, d)
  warmup <- round(p_adapt * iter, 0)
  pb_list <- list(rep(NA, chains))
  out_list <- foreach::foreach(
    c = 1:chains,
    .export = c(
      "find_reasonable_epsilon", "nuts_iteration",
      "build_tree", "leapfrog", "y", "sigma"
    )
  ) %dopar% {
    # Get initial values
    th <- starting_values[c, ]
    # Get a reasonable starting epsilon
    epsilon <- find_reasonable_epsilon(th, log_p_th, grad_th, M = M, verbose = verbose)
    mu_eps <- log(10 * epsilon)
    # store sims
    sims <- array(NA,
      dim = c(iter, d),
      dimnames = list(NULL, colnames(starting_values))
    )
    L <- rep(NA, iter)
    eps <- rep(NA, iter)
    for (t in 1:iter) {
      # Update with NUTS
      temp <- nuts_iteration(th, log_p_th, grad_th, epsilon,
        warmup = (t <= warmup), M = M, verbose = verbose,
        t = t, epsilon_bar = epsilon_bar, H = H, gamma = gamma, t0 = t0, kappa = kappa, delta = delta, mu_eps = mu_eps
      ) # for adapting epsilon
      L[t] <- temp$L
      sims[t, ] <- temp$th
      th <- temp$th
      eps[t] <- temp$epsilon
      epsilon <- temp$epsilon
      epsilon_bar <- temp$epsilon_bar
      H <- temp$H
    }
    list(sims = sims, L = L, eps = eps)
  }

  parallel::stopCluster(cl)

  sims_list <- lapply(out_list, function(i) i$sims)
  sims <- abind::abind(sims_list, along = -1) # has dimension (chains * iter * params)
  sims <- aperm(sims, perm = c(2, 1, 3)) # permutes to (iter * chains * dimension)
  L <- abind::abind(lapply(out_list, function(i) i$L), along = -1)
  eps <- abind::abind(lapply(out_list, function(i) i$eps), along = -1)
  cat("\n")
  rstan::monitor(sims, warmup)
  cat(
    "Avg trajectory length",
    round(rowMeans(L[, (warmup + 1):iter] * eps[, (warmup + 1):iter]), 2), "\n"
  )
  return(list(sims = sims, L = L, eps = eps))
}

#' Run NUTS with many chains
#'
#' @param log_p_th a function that takes a vector of parameters and returns a log likelihood
#' @param grad_th a function that takes a vector of parameters and returns the gradient
#' @param starting_values c x d matrix of starting values, c: number of chains, d: number of parameters
#' @param iter number of MCMC iterations
#' @param M covariance for the momentum
#' @param p_adapt the fragtion of iterations where step size epsilon is adapted
#' @param epsilon_bar initial guess for epsilon
#' @param verbose boolean indicating to print progress or not
#' @return
#'      sims: iter x c x d array of MCMC sample
#'      L: vector of trajectory length
#'      eps: vector of step size
#' @export
nuts_run <- function(log_p_th, grad_th, starting_values, iter, M = NULL,
                     p_adapt = 0.5, epsilon_bar = 1, H = 0, gamma = 0.05, t0 = 10, kappa = 0.75, delta = 0.65,
                     verbose = FALSE) {
  chains <- nrow(starting_values)
  d <- ncol(starting_values)
  if (is.null(M)) M <- rep(1, d)
  sims <- array(NA,
    dim = c(iter, chains, d),
    dimnames = list(NULL, NULL, colnames(starting_values))
  )
  warmup <- round(p_adapt * iter, 0)
  L <- array(NA, c(iter, chains))
  eps <- array(NA, c(iter, chains))
  pb_list <- list(rep(NA, chains))
  for (c in 1:chains) {
    # Set progress bar
    pb_list[[c]] <- txtProgressBar(style = 3)
    # Get initial values
    th <- starting_values[c, ]
    # Get a reasonable starting epsilon
    epsilon <- find_reasonable_epsilon(th, log_p_th, grad_th, M = M, verbose = verbose)
    mu_eps <- log(10 * epsilon)
    for (t in 1:iter) {
      # Update with NUTS
      temp <- nuts_iteration(th, log_p_th, grad_th, epsilon,
        warmup = (t <= warmup), M = M, verbose = verbose,
        t = t, epsilon_bar = epsilon_bar, H = H, gamma = gamma, t0 = t0, kappa = kappa, delta = delta, mu_eps = mu_eps
      ) # for adapting epsilon
      L[t, c] <- temp$L
      sims[t, c, ] <- temp$th
      th <- temp$th
      eps[t, c] <- temp$epsilon
      epsilon <- temp$epsilon
      epsilon_bar <- temp$epsilon_bar
      H <- temp$H

      # Update progress bar
      setTxtProgressBar(pb_list[[c]], t / iter)
    }
    cat("\033[1B") # draw new progress bar on a new line
  }
  cat("\n")
  rstan::monitor(sims, warmup)
  cat(
    "Avg trajectory length",
    round(colMeans(L[(warmup + 1):iter, ] * eps[(warmup + 1):iter, ]), 2), "\n"
  )
  return(list(sims = sims, L = L, eps = eps))
}

#' Run one NUTS iteration
#'
#' @param th a vector of parameters of interests
#' @param log_p_th a function that takes a vector of parameters and returns a log likelihood
#' @param grad_th a function that takes a vector of parameters and returns the gradient
#' @param epsilon step size
#' @param M covariance of the momentums
#' @param max_depth max depth for recursive leap-frog tree
#' @param warmup boolean of whether to update epsilon or not
#' @return List of outputs
#'       th: a new sample of the vector of parmeters
#'       epsilon: step size for the next iteration
#' @export
nuts_iteration <- function(th, log_p_th, grad_th, epsilon, M, max_depth = 5, warmup = FALSE, verbose = FALSE,
                           t = NULL, epsilon_bar = NULL, H = NULL, gamma = 0.05, t0 = 10, kappa = 0.75, delta = 0.65, mu_eps = NULL) {
  d <- length(th)
  # sample momentum
  r0 <- rnorm(d, 0, sqrt(M))
  # sample slice sample
  log_p_old <- log_p_th(th) - 0.5 * sum(r0^2 / M)
  log_u <- log(runif(1, 0, exp(log_p_old)))
  # initialize NUTs params
  th_new <- th
  th_minus <- th
  th_plus <- th
  r_minus <- r0
  r_plus <- r0
  n <- 1
  s <- 1
  j <- 0
  # get NUTS proposal by recursively building a binary tree
  while (s == 1) {
    # pick a direction at random
    v <- sample(c(-1, 1), 1, prob = c(0.5, 0.5))
    # build tree
    if (v == -1) {
      out <- build_tree(th_minus, r_minus, log_u, v, j, epsilon, log_p_old, log_p_th, grad_th, warmup = warmup, M = M)
      th_minus <- out$th_minus
      r_minus <- out$r_minus
    } else {
      out <- build_tree(th_plus, r_plus, log_u, v, j, epsilon, log_p_old, log_p_th, grad_th, warmup = warmup, M = M)
      th_plus <- out$th_plus
      r_plus <- out$r_plus
    }
    th_prime <- out$th_prime
    n_prime <- out$n_prime
    s_prime <- out$s_prime
    alpha <- out$alpha
    n_alpha <- out$n_alpha
    # get candidate proposal from C
    if ((s_prime == 1) & (runif(1) < min(1, n_prime / n))) {
      th_new <- th_prime
    }
    n <- n + n_prime
    s <- as.numeric(s_prime * (crossprod(th_plus - th_minus, r_minus) >= 0) * (crossprod(th_plus - th_minus, r_plus) >= 0))
    j <- j + 1
    if (j > max_depth & verbose) {
      warning("NUTS reached max tree depth, terminating early")
      break
    }
  }

  # Adapt epsilon
  if (warmup) {
    H <- (1 - 1 / (t + t0)) * H + 1 / (t + t0) * (delta - alpha / n_alpha)
    log_e <- mu_eps - (sqrt(t) / gamma) * H
    log_e_bar <- t^(-kappa) * log_e + (1 - t^(-kappa)) * log(epsilon_bar)
    epsilon <- exp(log_e)
    epsilon_bar <- exp(log_e_bar)
  } else {
    epsilon <- epsilon_bar
  }

  return(list(th = th_new, L = 2^(j - 1), epsilon = epsilon, epsilon_bar = epsilon_bar, H = H))
}


#' Proposes a candidate for th
build_tree <- function(th, r, log_u, v, j, epsilon, log_p_old, log_p_th, grad_th, warmup, M) {
  DELTA_MAX <- 1000
  # base case
  if (j == 0) {
    out <- leapfrog(th, r, v * epsilon, grad_th, M = M)
    th_prime <- out$th
    r_prime <- out$r
    log_p <- log_p_th(th_prime) - 0.5 * sum(r_prime^2 / M)
    n_prime <- 1 * (log_u <= log_p)
    s_prime <- 1 * (log_u < DELTA_MAX + log_p)
    if (s_prime == 0 & !warmup) {
      warning(paste("Divergent transition after warmup, current eps=", epsilon, "log_u=", log_u))
    }
    alpha <- min(1, exp(log_p - log_p_old))
    return(list(
      th_minus = th_prime, r_minus = r_prime,
      th_plus = th_prime, r_plus = r_prime,
      th_prime = th_prime, n_prime = n_prime, s_prime = s_prime,
      alpha = alpha, n_alpha = 1
    ))
  } else {
    # Recursion
    out0 <- build_tree(th, r, log_u, v, j - 1, epsilon, log_p_old, log_p_th, grad_th, warmup = warmup, M = M)
    th_minus <- out0$th_minus
    r_minus <- out0$r_minus
    th_plus <- out0$th_plus
    r_plus <- out0$r_plus
    th_prime <- out0$th_prime
    n_prime <- out0$n_prime
    s_prime <- out0$s_prime
    alpha_prime <- out0$alpha
    n_alpha_prime <- out0$n_alpha
    if (out0$s_prime == 1) {
      if (v == -1) {
        out1 <- build_tree(th_minus, r_minus, log_u, v, j - 1, epsilon, log_p_old, log_p_th, grad_th, warmup = warmup, M = M)
        th_minus <- out1$th_minus
        r_minus <- out1$r_minus
      } else {
        out1 <- build_tree(th_plus, r_plus, log_u, v, j - 1, epsilon, log_p_old, log_p_th, grad_th, warmup = warmup, M = M)
        th_plus <- out1$th_plus
        r_plus <- out1$r_plus
      }
      n_prime <- out0$n_prime + out1$n_prime
      if (n_prime > 0) {
        accept_p <- out1$n_prime / n_prime
        if (runif(1) < accept_p) {
          th_prime <- out1$th_prime
        }
      }
      s_prime <- as.numeric(out1$s_prime * (crossprod(th_plus - th_minus, r_minus) >= 0) * (crossprod(th_plus - th_minus, r_plus) >= 0))
      alpha_prime <- out0$alpha + out1$alpha
      n_alpha_prime <- out0$n_alpha + out1$n_alpha
    }
    return(list(
      th_minus = th_minus, r_minus = r_minus,
      th_plus = th_plus, r_plus = r_plus,
      th_prime = th_prime, n_prime = n_prime, s_prime = s_prime,
      alpha = alpha_prime, n_alpha = n_alpha_prime
    ))
  }
}


#' Leap frog
leapfrog <- function(th, r, epsilon, grad_th, M) {
  r <- r + 0.5 * epsilon * grad_th(th)
  th <- th + epsilon * r / M
  r <- r + 0.5 * epsilon * grad_th(th)
  return(list(th = th, r = r))
}

#' Calculate reasonable initial value for epsilon
#'
#' @param th a vector of parameters of interests
#' @param log_p_th a function that takes a vector of parameters and returns a log likelihood
#' @param grad_th a function that takes a vector of parameters and returns the gradient
#' @param M covariance of the momentums
#' @param max_iter max number of search iterationns
#' @return a value for epsilon
#' @export
find_reasonable_epsilon <- function(th, log_p_th, grad_th, M, max_iter = 100, verbose = FALSE) {
  epsilon <- 1
  d <- length(th)
  r <- rnorm(d)
  log_p <- log_p_th(th) - 0.5 * sum(r^2 / M)
  out <- leapfrog(th, r, epsilon, grad_th, M)
  log_p_new <- log_p_th(out$th) - 0.5 * sum(out$r^2 / M)
  log_r <- log_p_new - log_p
  a <- 2 * (log_r > log(0.5)) - 1
  i <- 1
  while (a * log_r > (-a) * log(2)) {
    epsilon <- (2^a) * epsilon
    out <- leapfrog(th, r, epsilon, grad_th, M)
    log_p_new <- log_p_th(out$th) - 0.5 * sum(out$r^2 / M)
    log_r <- log_p_new - log_p
    i <- i + 1
    if (i > max_iter & verbose) {
      warning("find_reasonable_epsilon reached max iter, terminating early")
      break
    }
  }
  if (verbose) message("Reasonable epsilon=", epsilon, " after ", i, " steps")
  return(epsilon)
}
