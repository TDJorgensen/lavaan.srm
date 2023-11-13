
  int<lower=0> Ng;  // number of groups  (Level 3)
  int IDg[Nd];      // group-level IDs for dyad-level observations

  matrix[Kd2, 3] rr_group_t;  // for prior of group-level RR-component SDs
  // beta hyperparameters for group-level correlations
  // alpha below diagonal, beta above diagonal (ignore diagonal)
  matrix<lower=0>[allKg, allKg] group_beta;

