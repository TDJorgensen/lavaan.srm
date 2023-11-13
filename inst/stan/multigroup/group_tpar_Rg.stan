  // assemble correlations among GROUP-level components of round-robin variables
  {
    int vec_idx = 1;    // index Rg_vec (vector of correlations)

    for (k in 1:allKg) {
      Rg[k,k] = 1;  // set diagonal = 1

      // column-wise along lower triangle
      // == row-wise along upper triangle
      if (k < allKg) { for (kk in (k+1):allKg) {
        Rg[kk,  k] = -1 + 2*Rg_vec[vec_idx];
        Rg[ k, kk] = -1 + 2*Rg_vec[vec_idx];
        vec_idx += 1;
      }}
    }

    // end block combining group-level correlation matrix
  }
