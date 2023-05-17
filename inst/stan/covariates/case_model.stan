  // SDs of case-level covariates (start indexing after AP effects)
  for (k in 1:Kp) {
    S_p[2*Kd2 + k] ~ normal(halfRangeKp[k] / 5, halfRangeKp[k] / 5) T[0, ];
  }
