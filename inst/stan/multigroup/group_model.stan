
  // priors for group-level SDs
  for (k in 1:Kd2) {
    S_g[k] ~ student_t(rr_group_t[k,1], rr_group_t[k,2], rr_group_t[k,3]) T[0, ];
  }

  // priors for group-level correlations
  {
    int vec_idx = 1;    // index Rg_vec (vector of correlations)

    for (k in 1:allKg) {
      // column-wise along lower triangle
      // == row-wise along upper triangle
      if (k < allKg) { for (kk in (k+1):allKg) {
        // hyperparameters:       alpha below diag,   beta above diag
        Rg_vec[vec_idx] ~ beta(  group_beta[kk, k],  group_beta[k, kk] );
        vec_idx += 1;
      }}
    }

    // end block for priors of group-level correlations
  }

  // priors for group-level random effects
  // FIXME: this assumes no group-level covariates
  for (n in 1:Ng) {
    GG[n,] ~ multi_normal_cholesky(rep_row_vector(0, Kd2), chol_g);
  }
