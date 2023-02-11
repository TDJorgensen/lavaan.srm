
// priors for means and SDs, based on empirical ranges
for (k in 1:Kd2) {
  S_g[k] ~ student_t(4, 0, halfRangeKd2[k] / 5) T[0, ];
}

// priors for correlations
chol_g ~ lkj_corr_cholesky(2);

// priors for random effects
for (n in 1:Ng) {
  GG[n,] ~ multi_normal_cholesky(rep_row_vector(0, Kd2), chol_g);
}
