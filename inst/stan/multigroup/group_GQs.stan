
  // % variance at group, case, and dyad levels
  matrix[Kd2, 4] Rsq;

  // group-level correlation matrix
  matrix[Kd2, Kd2] Rg;

  // calculate % of each round-robin variance due to each random effect
  {
    matrix[Kd2, 4] vars; // group, actor, partner, and relationship variances
    vector[Kd2] totals;  // sum of variance components

    for (k in 1:Kd2) {
      vars[k,1] = square(S_g[k]);       // group
      vars[k,2] = square(S_p[2*k - 1]); // actor
      vars[k,3] = square(S_p[2*k]);     // partner
      vars[k,4] = square(s_rr[k]);      // relationship

      totals[k] = sum(vars[k, ]);
      Rsq[k,] = vars[k,] ./ totals[k];
    }
    // end R-squared block
  }

  // calculate group-level correlation matrix
  Rg = multiply_lower_tri_self_transpose(chol_g);
