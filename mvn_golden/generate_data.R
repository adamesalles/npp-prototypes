library(MASS)

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

# save data to a file
saveRDS(list(
  D0 = D0,
  D = D,
  Sigma = Sigma,
  Sigma0 = Sigma0,
  alpha0 = alpha0,
  beta0 = beta0,
  true_theta = true_theta
), file = "../data/npp_mvn_data.rds")

# save D and D0 to csv
write.csv(D, file = "../data/D.csv", row.names = FALSE)
write.csv(D0, file = "../data/D0.csv", row.names = FALSE)