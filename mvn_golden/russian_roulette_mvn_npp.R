library(MASS)      # for mvrnorm
library(mvtnorm)   # for dmvnorm
library(cmdstanr)
library(posterior)
library(ggplot2)
library(future.apply)

set.seed(42)
# historical ("prior") sample
N0     <- 100;   p <- 3
sigma  <- 1.5;   sigma0 <- 10
Sigma  <- sigma^2  * diag(p)
Sigma0 <- sigma0^2 * diag(p)

true_theta <- rnorm(p, 0, sigma0)
D0 <- mvrnorm(N0, mu = true_theta, Sigma = Sigma)

# current data
N  <- 50
D  <- mvrnorm(N,  mu = true_theta, Sigma = Sigma)

# hyper‐parameters for eta ~ Beta(alpha0,beta0)
alpha0 <- 1;  beta0 <- 1

# pre-compute for power‐posterior draws
dbar  <- colMeans(D0)
invSigma  <- solve(Sigma)
invSigma0 <- solve(Sigma0)
S <- t(D0 - matrix(dbar, N0, p, byrow = TRUE)) %*% (D0 - matrix(dbar, N0, p, byrow = TRUE))

# gold‐standard via Stan (closed‐form c0)
stan_data <- list(
  p      = p,
  N0     = N0,
  N      = N,
  D0     = D0,
  D      = D,
  Sigma  = Sigma,
  Sigma0 = Sigma0,
  alpha0 = alpha0,
  beta0  = beta0
)
mod <- cmdstan_model("npp_mvn.stan")
fit_gold <- mod$sample(
  data             = stan_data,
  seed             = 123,
  chains           = 4,
  parallel_chains  = 4,
  iter_warmup      = 1000,
  iter_sampling    = 5000,
  refresh          = 0
)

draws_gold <- as_draws_df(fit_gold$draws(c("theta","eta")))

# draw exactly from the power‐posterior
sample_theta_power <- function(n, tau) {
  # posterior precision:  V*tau^{-1} = tau*N0*Sigma^{-1} + Sigma0^{-1}
  Q <- tau * N0 * invSigma + invSigma0
  V <- solve(Q)
  mu <- V %*% (tau * N0 * invSigma %*% dbar)
  return(mvrnorm(n, mu, V))
}

# full log‐likelihood of historical data (including constants)
log_lik0 <- function(D0, theta) {
  sum(dmvnorm(D0, mean = theta, sigma = Sigma, log = TRUE))
}

log_c0 <- function(eta) {
  p <- length(dbar)
  Sigma_inv <- ginv(Sigma)
  Sigma0_inv <- ginv(Sigma0)
  
  A <- eta * N0 * Sigma_inv + Sigma0_inv
  A_inv <- ginv(A)
  B <- eta * N0 * Sigma_inv %*% dbar
  
  term1 <- -eta * N0 * p / 2 * log(2 * pi)
  term2 <- -eta * N0 / 2 * log(det(Sigma))
  term3 <- -eta / 2 * sum(diag(S %*% Sigma_inv))
  term4 <- -eta * N0 / 2 * t(dbar) %*% Sigma_inv %*% dbar
  term5 <- t(B) %*% A_inv %*% B / 2
  term6 <- -0.5 * log(det(Sigma0)) - 0.5 * log(det(A))
  
  return(as.numeric(term1 + term2 + term3 + term4 + term5 + term6))
}

unbiased_ratio_exact <- function(eta0, eta1) {
  # Exact computation of c0(eta1) / c0(eta0)
  
  if (eta0 == eta1) return(1.0)
  
  log_c0_eta0 <- log_c0(eta0)
  log_c0_eta1 <- log_c0(eta1)
  
  
  return(exp(log_c0_eta1 - log_c0_eta0))
}

log_unbiased_ratio_RR <- function(eta0, eta1, K = 10, n=50) {
  # Numerical integration using K Monte Carlo samples
  
  if (eta0 == eta1) {
    return(1.0)
  } else if (eta0 > eta1) {
    etl <- eta1
    etu <- eta0
  } else {
    etl <- eta0
    etu <- eta1
  }
  
  tau_star <- runif(K, etl, etu)
  
  # Estimate E[log L0(theta) | τ_mid]
  interval_width <- (eta1 - eta0)
  
  estimate_integral <- function(tau) {
    # Sample from power posterior at tau_mid
    theta_star <- sample_theta_power(n, tau)
    log_lik_val <- mean(future_apply(theta_star, 1, function(th) {
      log_lik0(D0, th)
    }, future.seed=TRUE))
    return(log_lik_val)
  }
  # Russian roulette expansion of exp(integral_estimate)
  
  # integral_estimate <- 0
  # for (k in 1:K) {
  #   integral_estimate <- integral_estimate + estimate_integral(tau_star[k])
  # }
  
  integral_estimate <- sum(future_sapply(tau_star, FUN = estimate_integral, future.seed=TRUE))
  
  Delta <- integral_estimate * interval_width / K
  
  return(Delta)
}

unbiased_ratio_RR <- function(eta0, eta1, p_rr = 0.6, K = 200, n=200, M=100) {
  Delta <- log_unbiased_ratio_RR(eta0, eta1, K, n)
  
  # Russian‐roulette truncation for exp(Delta)
  if (abs(Delta) < 1e-30) {
    return(1.0)  # avoid numerical issues when eta0 approx. eta1
  }
  
  results <- numeric(M)
  for (m in 1:M) {
    # Geometric truncation
    N_rr <- rgeom(1, p_rr)
    
    # Compute truncated series: sum_{k=0}^{N_rr} Delta^k / k!
    results[m] <- 0
      for (k in 0:N_rr) {
        if (k == 0) {
          term <- 1
        } else {
          term <- Delta^k / factorial(k)
        }
        results[m] <- results[m] + term / (1 - p_rr)^k
      }
  }
  
  return(mean(results))
}

# random‐walk proposals (symmetric, cancels in ratio)
propose_theta <- function(theta, sigma_prop = 0.5) {
  as.numeric(mvrnorm(1, theta, sigma_prop^2 * diag(p)))
}

propose_eta <- function(eta, sigma_prop = 0.05) {
  x <- runif(1, eta - sigma_prop, eta + sigma_prop)
  if (x < 0 || x > 1) eta else x
}

# full MCMC step: updates theta then eta 
log_prior_theta <- function(theta) {
  sum(dmvnorm(theta, mean = rep(0,p), sigma = Sigma0, log = TRUE))
}
log_lik_data <- function(D, theta) {
  sum(dmvnorm(D, mean = theta, sigma = Sigma, log = TRUE))
}
log_prior_eta <- function(eta) {
  dbeta(eta, alpha0, beta0, log = TRUE)
}

run_joint_mcmc <- function(niter, init_theta, init_eta,
                           sigma_prop_theta = 0.5, sigma_prop_eta = 0.05,
                           p_rr = 0.5) {
  chain_th <- matrix(NA, nrow = niter, ncol = p)
  chain_et <- numeric(niter)
  chain_th[1,] <- init_theta
  chain_et[1] <- init_eta
  
  accept_theta <- 0
  accept_eta <- 0
  
  for (i in 2:niter) {
    # update theta via RW‐MH given eta
    th0 <- chain_th[i-1,]
    et0 <- chain_et[i-1]
    th1 <- propose_theta(th0, sigma_prop_theta)
    lp0 <- log_prior_theta(th0) + log_lik_data(D, th0) + et0 * log_lik0(D0, th0)
    lp1 <- log_prior_theta(th1) + log_lik_data(D, th1) + et0 * log_lik0(D0, th1)
    
    if (log(runif(1)) < lp1 - lp0) {
      chain_th[i,] <- th1
      accept_theta <- accept_theta + 1
    } else {
      chain_th[i,] <- th0
    }
    
    # update eta via pseudo-marginal MH 
    et0 <- chain_et[i-1]
    th_curr <- chain_th[i,]   # use updated theta
    et1 <- propose_eta(et0, sigma_prop_eta)
    
    # log‐prior ratio
    lpa0 <- log_prior_eta(et0)
    lpa1 <- log_prior_eta(et1)
    
    # tempered historical‐lik ratio
    hist_ratio <- (et1 - et0) * log_lik0(D0, th_curr)
    
    # unbiased estimator of c0(et1)/c0(et0) using Russian roulette
    # Rhat <- unbiased_ratio_exact(et0, et1)
    # Rhat <- unbiased_ratio_RR(et0, et1, p_rr)
    # logRhat <- log(Rhat)
    logRhat <- log_unbiased_ratio_RR(et0, et1)
    
    log_r <- (lpa1 - lpa0) + hist_ratio - logRhat
    if (log(runif(1)) < log_r) {
      chain_et[i] <- et1
      accept_eta <- accept_eta + 1
    } else {
      chain_et[i] <- et0
    }
    # cat("Iteration", i, ": eta = ", round(et0, 4), 
        # "->", round(et1, 4), ", log_r =", round(log_r, 6), "\n")
    
    # Progress reporting
    if (i %% 100 == 0) {
      cat("Iteration", i, "/ Accept rates: theta =", 
          round(accept_theta/i, 3), ", eta =", round(accept_eta/i, 3), "\n")
    }
  }
  
  cat("Final acceptance rates: theta =", round(accept_theta/niter, 3), 
      ", eta =", round(accept_eta/niter, 3), "\n")
  
  list(theta = chain_th, eta = chain_et)
}

# Test the Russian roulette estimator
test_rr <- FALSE
if (test_rr) {
  cat("Testing Russian roulette estimator...\n")
  set.seed(123)
  eta_test <- c(0.55, 0.6)
  rr_estimates <- replicate(10, log_unbiased_ratio_RR(eta_test[1], eta_test[2]))
  cat("RR estimates summary:\n")
  print(summary(rr_estimates))
  cat("Mean:", mean(rr_estimates), "SD:", sd(rr_estimates), "\n")
  cat("Exact ratio:", log(unbiased_ratio_exact(eta_test[1], eta_test[2])), "\n")
}


# Run MCMC
cat("\nRunning MCMC...\n")
set.seed(123)
mcmc_res <- run_joint_mcmc(niter = 10000,
                           init_theta = rep(0, p),
                           init_eta   = 0.5,
                           sigma_prop_theta = 0.25,
                           sigma_prop_eta   = 0.04,
                           p_rr = 0.7)

# Analysis
burnin <- 1000
eta_samples <- mcmc_res$eta[-(1:burnin)]

cat("\nPosterior summary for eta:\n")
cat("Mean:", mean(eta_samples), "\n")
cat("Median:", median(eta_samples), "\n")
cat("95% CI:", quantile(eta_samples, c(0.025, 0.975)), "\n")

# Extract theta samples from MCMC (after burnin)
theta_samples <- mcmc_res$theta[-(1:burnin), ]

# Create comparison plots
create_comparison_plots <- function(draws_gold, mcmc_res, burnin = 1000) {
  
  # Extract samples after burnin
  eta_rr <- mcmc_res$eta[-(1:burnin)]
  theta_rr <- mcmc_res$theta[-(1:burnin), ]
  # make theta_rr work as 2d matrix
  if (is.null(dim(theta_rr))) {
    theta_rr <- matrix(theta_rr, ncol = 1)
  }
  
  # Prepare data for eta comparison
  eta_comparison <- data.frame(
    value = c(draws_gold$eta, eta_rr),
    method = rep(c("Stan Gold", "Russian Roulette"), 
                 times = c(length(draws_gold$eta), length(eta_rr)))
  )
  
  # Plot eta density comparison
  p_eta <- ggplot(eta_comparison, aes(x = value, fill = method, color = method)) +
    geom_density(alpha = 0.5, size = 1) +
    scale_fill_manual(values = c("Stan Gold" = "blue", "Russian Roulette" = "red")) +
    scale_color_manual(values = c("Stan Gold" = "blue", "Russian Roulette" = "red")) +
    labs(title = "Posterior Density Comparison: eta",
         x = "eta", y = "Density",
         fill = "Method", color = "Method") +
    theme_minimal() +
    theme(legend.position = "top")
  
  # Create theta comparison plots for each dimension
  theta_plots <- list()
  J <- ncol(theta_rr)
  for (j in 1:J) {
    # Extract theta[j] from Stan results
    theta_gold_j <- draws_gold[[paste0("theta[", j, "]")]]
    
    theta_j_comparison <- data.frame(
      value = c(theta_gold_j, theta_rr[, j]),
      method = rep(c("Stan Gold", "Russian Roulette"), 
                   times = c(length(theta_gold_j), nrow(theta_rr)))
    )
    
    theta_plots[[j]] <- ggplot(theta_j_comparison, aes(x = value, fill = method, color = method)) +
      geom_density(alpha = 0.5, size = 1) +
      scale_fill_manual(values = c("Stan Gold" = "blue", "Russian Roulette" = "red")) +
      scale_color_manual(values = c("Stan Gold" = "blue", "Russian Roulette" = "red")) +
      labs(title = paste("Posterior Density: theta[", j, "]", sep = ""),
           x = paste("theta[", j, "]", sep = ""), y = "Density",
           fill = "Method", color = "Method") +
      theme_minimal() +
      theme(legend.position = "top")
  }
  
  # Print summary statistics
  cat("\n=== COMPARISON SUMMARY ===\n")
  cat("ETA:\n")
  cat("  Stan Gold    - Mean:", round(mean(draws_gold$eta), 4), 
      "SD:", round(sd(draws_gold$eta), 4), "\n")
  cat("  Russian Roul - Mean:", round(mean(eta_rr), 4), 
      "SD:", round(sd(eta_rr), 4), "\n")
  cat("  Difference in means:", round(abs(mean(draws_gold$eta) - mean(eta_rr)), 4), "\n\n")
  
  for (j in 1:J) {
    theta_gold_j <- draws_gold[[paste0("theta[", j, "]")]]
    cat("THETA[", j, "]:\n")
    cat("  Stan Gold    - Mean:", round(mean(theta_gold_j), 4), 
        "SD:", round(sd(theta_gold_j), 4), "\n")
    cat("  Russian Roul - Mean:", round(mean(theta_rr[, j]), 4), 
        "SD:", round(sd(theta_rr[, j]), 4), "\n")
    cat("  Difference in means:", round(abs(mean(theta_gold_j) - mean(theta_rr[, j])), 4), "\n\n")
  }
  
  return(list(eta_plot = p_eta, theta_plots = theta_plots))
}

# Simple diagnostic plots (always available)
plot_df <- data.frame(
  iteration = (burnin+1):length(mcmc_res$eta),
  eta = eta_samples
)

# Trace plot for eta
p_trace_eta <- ggplot(plot_df, aes(x = iteration, y = eta)) +
  geom_line(alpha = 0.7, color = "darkblue") +
  labs(title = "Trace plot of eta (Russian Roulette MCMC)", 
       x = "Iteration", y = "eta") +
  theme_minimal()

# Trace plots for theta
theta_trace_data <- data.frame(
  iteration = rep((burnin+1):nrow(mcmc_res$theta), p),
  theta_value = as.vector(theta_samples),
  parameter = rep(paste0("theta[", 1:p, "]"), each = nrow(theta_samples))
)

p_trace_theta <- ggplot(theta_trace_data, aes(x = iteration, y = theta_value)) +
  geom_line(alpha = 0.7, color = "darkred") +
  facet_wrap(~ parameter, scales = "free_y", ncol = 1) +
  labs(title = "Trace plots of theta (Russian Roulette MCMC)",
       x = "Iteration", y = "theta value") +
  theme_minimal()

# Density plot for eta (RR only)
p_density_eta <- ggplot(plot_df, aes(x = eta)) +
  geom_density(fill = "skyblue", alpha = 0.7, color = "darkblue") +
  geom_vline(xintercept = mean(eta_samples), color = "red", linetype = "dashed") +
  labs(title = "Posterior density of eta (Russian Roulette)", 
       x = "eta", y = "Density") +
  theme_minimal()

p_density_golden_eta <- ggplot(draws_gold, aes(x = eta)) +
  geom_density(fill = "skyblue", alpha = 0.7, color = "darkblue") +
  geom_vline(xintercept = mean(draws_gold$eta), color = "red", linetype = "dashed") +
  labs(title = "Posterior density of eta (Stan Golden)", 
       x = "eta", y = "Density") +
  theme_minimal()

print(p_trace_eta)
print(p_trace_theta)
print(p_density_eta)
print(p_density_golden_eta)

comparison_plots <- create_comparison_plots(draws_gold, mcmc_res, burnin = 1000)
print(comparison_plots$eta_plot)
if (length(comparison_plots$theta_plots) == 0) {
  cat("No theta plots to display.\n")
} else {
  cat("Displaying theta plots...\n")
}
for(i in 1:length(comparison_plots$theta_plots)) {
  print(comparison_plots$theta_plots[[i]])
}


# Test the unbiased ratio estimator with various eta pairs
unbiased_ratio_test <- FALSE
if (unbiased_ratio_test) {
  test_eta_pairs <- list(
    c(0.5, 0.6),
    c(0.3, 0.7),
    c(0.1, 0.9),
    c(0.8, 0.2)  # This should give inverse ratio
  )
  
  
  for (pair in test_eta_pairs) {
    eta0 <- pair[1]
    eta1 <- pair[2]
    
    # Multiple estimates
    rr_ests <- replicate(20, unbiased_ratio_exact(eta0, eta1, p_rr = 0.7))
    
    cat("eta0 =", eta0, ", eta1 =", eta1, "\n")
    cat("  RR estimates: mean =", round(mean(rr_ests), 4), 
        ", sd =", round(sd(rr_ests), 4), "\n")
    cat("  Range: [", round(min(rr_ests), 4), ",", round(max(rr_ests), 4), "]\n")
    
    # Check if any are extreme
    if (any(rr_ests > 100) || any(rr_ests < 0.01)) {
      cat("  WARNING: Extreme values detected!\n")
    }
    cat("\n")
  }
}