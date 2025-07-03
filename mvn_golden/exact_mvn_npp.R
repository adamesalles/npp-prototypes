library(MASS)
library(cmdstanr)
library(posterior)
library(ggplot2)    

set.seed(42)
# historical ("prior") sample
N0    <- 100
p     <- 4
sigma <- 1.5
sigma0<- 10
Sigma  <- sigma^2  * diag(p)
Sigma0 <- sigma0^2 * diag(p)

true_theta0 <- rnorm(p, 0, sigma0)
true_theta <- rnorm(p, 0, sigma0)
D0 <- mvrnorm(N0, mu = true_theta0, Sigma = Sigma)

# current data
N <- 50
D  <- mvrnorm(N,  mu = true_theta, Sigma = Sigma)

# hyper‐parameters for eta ~ Beta(alpha0,beta0)
alpha0 <- 1
beta0  <- 1


# assemble data list
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


# compile the Stan model
mod <- cmdstan_model("npp_mvn.stan")


# run the sampler
fit <- mod$sample(
  data             = stan_data,
  seed             = 123,
  chains           = 4,
  parallel_chains  = 4,
  iter_warmup      = 1000,
  iter_sampling    = 5000,
  refresh          = 500
)


# print a quick summary 
cat("True theta: ", true_theta, "\n")
fit$print(digits = 2)


# extract draws and make density plots
draws_df <- as_draws_df(fit$draws(c("theta", "eta")))

p1 <- ggplot(draws_df, aes(x = `theta[1]`)) +
  geom_density() +
  labs(title = expression(paste("Posterior of ", theta[1])))

p2 <- ggplot(draws_df, aes(x = `theta[2]`)) +
  geom_density() +
  labs(title = expression(paste("Posterior of ", theta[2])))

p3 <- ggplot(draws_df, aes(x = `theta[3]`)) +
  geom_density() +
  labs(title = expression(paste("Posterior of ", theta[3])))

p4 <- ggplot(draws_df, aes(x = eta)) +
  geom_density() +
  labs(title = expression(paste("Posterior of ", eta)))

# combine
(p1 | p2) / (p3 | p4)

# summary statistics
summarise_draws(draws_df)

# plot mcmc trace
t1 <- ggplot(draws_df, aes(x = `.iteration`, y = `theta[1]`)) +
  geom_line() +
  labs(title = expression(paste("Trace plot of ", theta[1])),
       x = "Iteration",
       y = expression(paste(theta[1]))) +
  theme_minimal()

t2 <- ggplot(draws_df, aes(x = `.iteration`, y = `theta[2]`)) +
  geom_line() +
  labs(title = expression(paste("Trace plot of ", theta[2])),
       x = "Iteration",
       y = expression(paste(theta[2]))) +
  theme_minimal()

t3 <- ggplot(draws_df, aes(x = `.iteration`, y = `theta[3]`)) +
  geom_line() +
  labs(title = expression(paste("Trace plot of ", theta[3])),
       x = "Iteration",
       y = expression(paste(theta[3]))) +
  theme_minimal()

t4 <- ggplot(draws_df, aes(x = `.iteration`, y = `eta`)) +
  geom_line() +
  labs(title = expression(paste("Trace plot of ", eta)),
       x = "Iteration",
       y = expression(paste(eta))) +
  theme_minimal()

(t1 | t2) / (t3 | t4)

