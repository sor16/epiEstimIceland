data {
  int<lower = 0> N_days;
  int<lower = 0> local[N_days]; // Local new cases
  vector<lower = 0>[N_days] total; // Total new cases
  vector<lower=0>[N_days] prop_quarantine;
  real<lower = 0> SI_shape; // 1.54 calculated from Icelandic data (Flaxman et al. use 4.5)
  real<lower = 0> SI_rate; // 0.28 calculated from Icelandic data (Flaxman et al. use 0.6)
}

transformed data {
  vector[N_days] SI_rev;
  vector[N_days - 1] lambda;
  vector[N_days - 1] log_lambda;
  
  // Calculate SI but put results in reverse vector for easy dot product calculations
  for (t in 1:N_days) {
  SI_rev[N_days - t + 1] = gamma_cdf(1.0 * t + 0.5, SI_shape, SI_rate) - gamma_cdf(1.0 * t - 0.5, SI_shape, SI_rate);
  }
  SI_rev = SI_rev / sum(SI_rev);
  
  // Calculate expectd number of new cases if R_t = 1 by convolving SI with past total new cases
  // The division is used such that the truncated SI at the start of the data sums to 1
  for (t in 2:N_days) {
    lambda[t - 1] = dot_product(head((1-prop_quarantine).*total, t - 1), tail(SI_rev, t - 1) / sum(tail(SI_rev, t - 1)));
  }
  
  log_lambda = log(lambda);
}

parameters {
  // Want to put normal prior on sqrt(phi) where Var(x) = mu + phi * mu^2
  // Stan's NB family is parametrized as Var(x) = mu + 1/phi * mu^2 so we have to use inverse of phi.
  real<lower = 0> phi_inv_sqrt;
  
  // Standard deviation of random jumps in second order random walk for R_t
  real<lower = 0> sigma_R;
  
  // Need to estimate first two R_t for second order random walk
  real log_R1;
  real log_R2;
  // Auxiliary variables for non-centered parametrization 
  // x = mu + sigma * z, where z ~ normal(0, 1)
  // instead of x ~ normal(mu, sigma))
  vector[N_days - 3] error;
} 

transformed parameters {
  // Model log(R_t) with random walk as log(R_t) \in (-inf, inf)
  vector[N_days - 1] log_R;
  
  
  real<lower = 0> phi_inv = square(phi_inv_sqrt);
  real<lower = 0> phi = inv(phi_inv);
  
  log_R[1] = log_R1;
  log_R[2] = log_R2;
  
  // Second order random walk assumes normally distributed second derivative
  // dfdt = f(t) - f(t - 1)
  // df2dt2 = (f(t) - f(t - 1)) - (f(t - 1) - f(t - 2)) = f(t) - 2f(t - 1) + f(t - 2) ~ normal(0, sigma)
  // Instead write f(t) ~ normal(2f(t - 1) - f(t - 2), sigma)
  for (t in 3:(N_days - 1)) {
    log_R[t] = 2 * log_R[t - 1] - log_R[t - 2] + sigma_R * error[t - 2];
  }
  
}

model {
  sigma_R ~ std_normal();
  error ~ std_normal();
  phi_inv_sqrt ~ std_normal();
  
  log_R1 ~ std_normal();
  log_R2 ~ std_normal();
  
  local[2:N_days] ~ neg_binomial_2_log(log_R + log_lambda, phi);
}

generated quantities {
  // Calculate R_t and y_hat for model checking
  vector[N_days - 1] R = exp(log_R);
  int y_hat[N_days - 1] = neg_binomial_2_log_rng(log_R + log_lambda, phi);
}




