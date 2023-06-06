
  matrix[Ng, allKg] GG;                 // group-level random effects
  vector<lower=0>[allKg] S_g;           // SDs of group effects (+ covariates)
  cholesky_factor_corr[allKg] chol_r_g; // Cholesky factor of group-level correlations
