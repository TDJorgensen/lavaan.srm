
  matrix[Ng, Kd2] GG;               // group-level random effects
  vector<lower=0>[Kd2] S_g;         // SDs of group effects (+ covariates)
  cholesky_factor_corr[Kd2] chol_g; // Cholesky factor of group-level correlations
