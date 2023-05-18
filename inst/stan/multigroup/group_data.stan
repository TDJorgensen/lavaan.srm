
  int<lower=0> Ng;  // number of groups  (Level 3)
  int IDg[Nd];      // group-level IDs

  matrix[Kd2, 3] rr_group_t;  // for prior of group-level RR component
  real group_lkj;             // LKJ parameter for case-level correlations

