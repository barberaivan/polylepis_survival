data {
  int<lower=0> N;               // trees
  int<lower=0> K;               // number of fixed-effects parameters
  int<lower=0> Ks;
  int<lower=0> P;               // number of plots  (random effect)

  vector[N] y;                  // growth
  matrix[N, K] X;               // fixed-effects design matrix
  matrix[N, Ks] Xs;             // fixed-effects design matrix for sigma
  matrix[N, P] Zp;              // random-effects design matrix

  // priors
  real prior_intercept_sd;
  real prior_b_sd;
  real prior_sigma_sd;

  real prior_alphasig_sd;
  real prior_alphasig_mean;
  real prior_betasig_sd;
}

parameters {
  vector[K] b_raw;                   // fixed parameters
  vector[Ks] bs_raw;
  vector[P] e_plot_raw;              // plot random effects
  real<lower=0> sigma_plot_raw;      // raneff sd
}

transformed parameters {
  vector[K] b;                   // fixed parameters
  vector[Ks] bs;
  vector[P] e_plot;              // plot random effects
  real<lower=0> sigma_plot;      // raneff sd

  b[1] = b_raw[1] * prior_intercept_sd;
  b[2:K] = b_raw[2:K] * prior_b_sd;
  sigma_plot = sigma_plot_raw * prior_sigma_sd;
  e_plot = e_plot_raw * sigma_plot;

  bs[1] = bs_raw[1] * prior_alphasig_sd + prior_alphasig_mean;
  bs[2] = bs_raw[2] * prior_betasig_sd;
}

model {
  // priors
  b_raw ~ std_normal();
  sigma_plot_raw ~ std_normal();
  e_plot_raw ~ std_normal();
  bs_raw ~ std_normal();

  // likelihood
  {
    vector[N] mu = X * b + Zp * e_plot;
    vector[N] sigma = exp(Xs * bs);
    y ~ normal(mu, sigma);
  }
}