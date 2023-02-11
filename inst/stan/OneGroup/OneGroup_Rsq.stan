
// % variance at case and dyad levels
matrix[Kd2, 3] Rsq;

// calculate % of each round-robin variance due to each random effect
{
  matrix[Kd2, 3] vars; // actor, partner, and relationship variances
  vector[Kd2] totals;  // sum of variance components

  for (k in 1:Kd2) {
    vars[k,1] = square(S_p[2*k - 1]); // actor
    vars[k,2] = square(S_p[2*k]);     // partner
    vars[k,3] = square(s_rr[k]);      // relationship

    totals[k] = sum(vars[k, ]);
    Rsq[k,] = vars[k,] ./ totals[k];
  }
  // end R-squared block
}

