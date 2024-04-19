data {
  int<lower=0> N;               // trees
  int<lower=0> K;               // number of fixed-effects parameters
  int<lower=0> Ks;              // number of fixed-effects parameters for scale
  int<lower=0> P;               // number of plots  (random effect)

  vector[N] y;                  // dieback at logit scale
  matrix[N, K] X;               // fixed-effects design matrix
  matrix[N, Ks] Xs;             // fixed-effects design matrix for scale
  matrix[N, P] Zp;              // random-effects design matrix

  // priors
  real prior_intercept_sd;
  real prior_b_sd;
  real prior_sigma_sd;
  real prior_scale_sd;
  real prior_alpha_sd;
}

parameters {
  vector[K] b_raw;                   // fixed parameters
  vector[P] e_plot_raw;              // plot random effects
  real<lower=0> sigma_plot_raw;      // raneff sd
  vector<lower=0>[Ks] scale_raw;     // dispersion
  real alpha_raw;
}

transformed parameters {
  vector[K] b;                   // fixed parameters
  vector[P] e_plot;              // plot random effects
  real<lower=0> sigma_plot;      // raneff sd
  vector<lower=0>[Ks] scale;     // dispersion
  real alpha;

  b[1] = b_raw[1] * prior_intercept_sd;
  b[2:K] = b_raw[2:K] * prior_b_sd;
  sigma_plot = sigma_plot_raw * prior_sigma_sd;
  e_plot = e_plot_raw * sigma_plot;

  scale = scale_raw * prior_scale_sd;
  alpha = alpha_raw * prior_alpha_sd;
}

model {
  // priors
  b_raw ~ std_normal();
  sigma_plot_raw ~ std_normal();
  e_plot_raw ~ std_normal();
  scale_raw ~ std_normal();
  alpha_raw ~ std_normal();

  // likelihood
  {
    vector[N] xi = X * b + Zp * e_plot;
    vector[N] scale_vec = Xs * scale;
    y ~ skew_normal(xi, scale_vec, alpha);
  }
}