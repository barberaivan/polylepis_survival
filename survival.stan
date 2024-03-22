data {
  int<lower=0> N;               // identified trees
  int<lower=0> N_unid_rows;     // unidentified rows (non-recorded in parcels with unidentified)
  int<lower=0> N_unid_trees;    // unidentified trees (recorded, not identified)

  int<lower=0> K;               // number of fixed-effects parameters
  int<lower=0> P;               // number of plots  (random effect)
  int<lower=0> P_unid;          // number of parcels with unidentified trees

  int y[N];                     // survival {0, 1}
  int y_unid[N_unid_trees];     // survival of unidentified trees {0, 1}

  // design matrices
  matrix[N, K] X;                // fixed effects
  matrix[N_unid_rows, K] X_unid;

  matrix[N, P] Zp;               // random effects
  matrix[N_unid_rows, P] Zp_unid;

  // data for marginalization over unidentified trees
  int y_unid_start[P_unid];      // beginning and end of each plot in y_unid
  int y_unid_end[P_unid];
  int y_unid_length[P_unid];

  int max_combs;                 // largest number of possible id distributions

  int combs_rows_matrix[N_unid_trees, max_combs];
  // matrix indicating the row_ids from X_unid that correpond to each combination
  // of possible ids.
  int combs_n[P_unid];           // number of possible combinations by plot

  int unid;                      // use unidentified trees in likelihood?

  // priors
  real prior_intercept_sd;
  real prior_b_sd;
  real prior_sigma_sd;
}

parameters {
  vector[K] b_raw;                   // fixed parameters
  vector[P] e_plot_raw;              // plot random effects
  real<lower=0> sigma_plot_raw;      // raneff sd
}

transformed parameters {
  vector[K] b;                   // fixed parameters
  vector[P] e_plot;              // plot random effects
  real<lower=0> sigma_plot;      // raneff sd

  b[1] = b_raw[1] * prior_intercept_sd;
  b[2:K] = b_raw[2:K] * prior_b_sd;
  sigma_plot = sigma_plot_raw * prior_sigma_sd;
  e_plot = e_plot_raw * sigma_plot;
}

model {
  // priors
  b_raw ~ std_normal();
  sigma_plot_raw ~ std_normal();
  e_plot_raw ~ std_normal();

  // likelihood for identified trees
  y ~ bernoulli_logit(X * b + Zp * e_plot);

  // likelihood for unidentified trees
  if(unid) {
    vector[P_unid] marginal_ll; // marginal loglik at each plot
    vector[N_unid_rows] p_logit = X_unid * b +
                                  Zp_unid * e_plot;
                                // logit-prob at each possible tree

    // loop over plots with unidentified trees
    for(p in 1:P_unid) {
      int nc = combs_n[p];
      int ylen = y_unid_length[p];
      int rows_response[ylen]; //= y_unid_start[p] : y_unid_end[p];
      vector[nc] loglik_combs;

      for(i in 1:ylen) {
        rows_response[i] = y_unid_start[p] + i - 1;
      }

      // loglik for every possible combination at plot p
      for(c in 1:nc) {
        int rows_pred[y_unid_length[p]] = combs_rows_matrix[rows_response, c];
        loglik_combs[c] = bernoulli_logit_lpmf(y_unid[rows_response] |
                                               p_logit[rows_pred]);
      }

      // add marginal likelihood to the log-posterior
      target += log_sum_exp(loglik_combs);
    }
  }
}