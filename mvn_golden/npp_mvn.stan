functions {
  // log normalizing constant for the power prior on D0
  real log_c0(real eta,
              int N0,
              matrix Sigma,
              matrix Sigma0,
              vector dbar,
              matrix S) {
    int p = num_elements(dbar);
    // inverses
    matrix[p,p] Sigma_inv  = inverse(Sigma);
    matrix[p,p] Sigma0_inv = inverse(Sigma0);
    matrix[p,p] A     = eta * N0 * Sigma_inv + Sigma0_inv;
    matrix[p,p] A_inv = inverse(A);
    vector[p] B = eta * N0 * (Sigma_inv * dbar);

    // terms from the closed-form integral
    real term1 = - eta * N0 * p / 2 * log(2 * pi());
    real term2 = - eta * N0 / 2 * log_determinant(Sigma);
    real term3 = - eta / 2 * trace(S * Sigma_inv);
    real term4 = - eta * N0 / 2 * quad_form(Sigma_inv, dbar);
    real term5 = quad_form(A_inv, B) / 2;
    real term6 = - 0.5 * log_determinant(Sigma0) 
                 - 0.5 * log_determinant(A);

    return term1 + term2 + term3 + term4 + term5 + term6;
  }
}
data {
  int<lower=1> p;
  int<lower=0> N0;
  int<lower=0> N;
  array[N0] vector[p] D0;
  array[N] vector[p] D;
  matrix[p, p] Sigma;
  matrix[p, p] Sigma0;
  real<lower=0> alpha0;
  real<lower=0> beta0;
}
transformed data {
  vector[p] dbar;
  matrix[p,p] S;
  {
    vector[p] sum_d0 = rep_vector(0.0, p);
    for (i in 1:N0)
      sum_d0 += D0[i];
    dbar = sum_d0 / N0;
  }
  {
    S = rep_matrix(0.0, p, p);
    for (i in 1:N0) {
      vector[p] diff = D0[i] - dbar;
      S += diff * diff';
    }
  }
}
parameters {
  vector[p] theta;
  real<lower=0,upper=1> eta;
}
model {
  // data likelihood
  target += multi_normal_lpdf(D | theta, Sigma);
  // tempered historical likelihood
  target += eta * multi_normal_lpdf(D0 | theta, Sigma);
  // prior on theta
  target += multi_normal_lpdf(theta | rep_vector(0.0, p), Sigma0);
  // prior on eta
  target += beta_lpdf(eta | alpha0, beta0);
  // subtract log c0(eta)
  target += - log_c0(eta, N0, Sigma, Sigma0, dbar, S);
}
