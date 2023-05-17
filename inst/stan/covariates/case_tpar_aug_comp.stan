  // augmented case-level data include random effects and covariates
  for (n in 1:Np) {
    // first 2*Kd2 columns == AP
    augYp[n, 1:(2*Kd2)] = AP[n, ]; // AP already sorted by case-ID
    for (k in 1:Kp) {
      // Yp might not be sorted, so use its IDs to assign the correct row
      augYp[n, 2*Kd2 + k] = Yp[ IDpp[n], k];
    }
  }
