data {
  int<lower=0> N;
  int y[N]; // burned plots
  int S[N]; // total plots
}

parameters {
  vector[N] eta; // burn probabilities
}

transformed parameters {
  vector[N] p = inv_logit(eta);
}

model {
  eta ~ normal(0, 10);
  y ~ binomial(S, p);
}