
  matrix[Ng, Kd2] GG;           // group-level random effects
  vector<lower=0>[allKg] S_g;   // SDs of group effects (+ covariates)
  // correlations among group-level random effects (+ covariates)
  vector<lower=0,upper=1>[ (allKg*(allKg-1) ) %/% 2] Rg_vec;
