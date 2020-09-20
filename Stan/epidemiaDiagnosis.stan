data {
    int<lower = 0> N_days;
    int<lower = 0> local[N_days]; // Local new cases
    vector<lower = 0>[N_days] total; // Total new cases
    real<lower = 0> SI_shape; // 1.54 calculated from Icelandic data (Flaxman et al. use 4.5)
    real<lower = 0> SI_rate; // 0.28 calculated from Icelandic data (Flaxman et al. use 0.6)
    real<lower = 0> pi_shape; // 1.54 calculated from Icelandic data (Flaxman et al. use 4.5)
    real<lower = 0> pi_rate; // 0.28 calculated from Icelandic data (Flaxman et al. use 0.6)
}

transformed data {
    vector[N_days] SI_rev;
    vector[N_days] pi_rev;
    
    // Calculate SI and pi but put results in reverse vector for easy dot product calculations
    for (t in 1:N_days) {
        SI_rev[N_days - t + 1] = gamma_cdf(1.0 * t + 0.5, SI_shape, SI_rate) - gamma_cdf(1.0 * t - 0.5, SI_shape, SI_rate);
        pi_rev[N_days - t + 1] = gamma_cdf(1.0 * t + 0.5, pi_shape, pi_rate) - gamma_cdf(1.0 * t - 0.5, pi_shape, pi_rate);
    }
    SI_rev = SI_rev / sum(SI_rev);
    pi_rev = pi_rev / sum(pi_rev);
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
  vector<lower=0>[N_days] log_mu_cases;
  vector<lower=0>[N_days] infections;
  // Model log(R_t) with random walk as log(R_t) \in (-inf, inf)
  vector[N_days - 1] log_R;
  real<lower = 0> phi_inv = square(phi_inv_sqrt);
  real<lower = 0> phi = inv(phi_inv);
  
  
  log_mu_cases[1]=log(total[1]);
  infections[1]=total[1];
  log_R[1] = log_R1;
  log_R[2] = log_R2;
  
  // Second order random walk assumes normally distributed second derivative
  // dfdt = f(t) - f(t - 1)
  // df2dt2 = (f(t) - f(t - 1)) - (f(t - 1) - f(t - 2)) = f(t) - 2f(t - 1) + f(t - 2) ~ normal(0, sigma)
  // Instead write f(t) ~ normal(2f(t - 1) - f(t - 2), sigma)
  for (t in 3:(N_days - 1)) {
    log_R[t] = 2 * log_R[t - 1] - log_R[t - 2] + sigma_R * error[t - 2];
  }
  real lambda_g;
  for (t in 2:N_days) {
      lambda_g = dot_product(head(infections, t - 1), tail(SI_rev, t - 1)) / sum(tail(SI_rev, t - 1));
      infections[t]=exp(log_R[t]) * lambda_g;
      log_mu_cases[t]=log(dot_product(head(infections, t - 1), tail(pi_rev, t - 1)) / sum(tail(pi_rev, t - 1)));
  }
}

model {
  sigma_R ~ std_normal();
  error ~ std_normal();
  phi_inv_sqrt ~ std_normal();
  
  log_R1 ~ std_normal();
  log_R2 ~ std_normal();
  
  local[2:N_days] ~ neg_binomial_2_log(log_mu_cases, phi);
}

generated quantities {
  // Calculate R_t and y_hat for model checking
  vector[N_days - 1] R = exp(log_R);
  int y_hat[N_days - 1] = neg_binomial_2_log_rng(log_mu_cases, phi);
  
}




