  // scale cholesky for case-level covariates
  {
    vector[allKp]         sub_S_p = rep_vector(1, allKp);
    matrix[allKp, allKp] chol_r_p = cholesky_decompose(Rp);

    for (k in (2*Kd2+1):allKp) sub_S_p[k] = S_p[k];
    chol_p = diag_pre_multiply(sub_S_p, chol_r_p);
  }
