
  // priors for means and SDs, based on empirical ranges
  for (k in 1:Kd2) {
    S_g[k] ~ student_t(rr_group_t[k,1], rr_group_t[k,2], rr_group_t[k,3]) T[0, ];
  }

  // priors for correlations
  chol_g ~ lkj_corr_cholesky(group_lkj);

  // priors for random effects
  // FIXME: this assumes no group-level covariates
  for (n in 1:Ng) {
    GG[n,] ~ multi_normal_cholesky(rep_row_vector(0, Kd2), chol_g);
  }
