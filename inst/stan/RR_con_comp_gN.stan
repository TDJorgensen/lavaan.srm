// Terrence D. Jorgensen
// Last updated: 16 May 2023

// Program to estimate covariance matrices for multivariate SRM components
// - standard round-robin design (indistinguishable subjects)
//    - multiple RR groups have random effects
// - continuous data only
// - no missing data


functions {}
data {
  // sample sizes
  int<lower=0> Nd;        // number of dyads   (Level 1)
  int<lower=0> Np;        // number of persons (Level 2, cross-classified)
  // number of observed measures (half the number of columns)
  int<lower=0> Kd2;       // number of dyad-level measures   (vary within dyad)
  // observed data
  matrix[Nd, 2*Kd2] Yd2;  // observed round-robin variables
  // ID variables
  int IDp[Nd, 2];         // person-level IDs (cross-classified)

#include /multigroup/group_data.stan

  // TODO: pass hyperparameters for priors as data?

}
transformed data {
  // number of pairs of round-robin variables
  int<lower=0> nPairs = Kd2*(Kd2 - 1) / 2;

#include /vanilla/tdata.stan

  // define counts of missing observations
  // int nMiss_d2 = 0;  // round-robin variables

  // save limits for default priors
  matrix[Kd2, 2] limKd2;
  vector[Kd2] halfRangeKd2;

  // calculate observed limits
  for (k in 1:Kd2) {
    limKd2[k,1] = min(Yd2[ : , (2*k - 1):(2*k)] );
    limKd2[k,2] = max(Yd2[ : , (2*k - 1):(2*k)] );
    halfRangeKd2[k] = (limKd2[k, 2] - limKd2[k, 1]) / 2;
  }

}
parameters {
  // means
#include /means/means_parameters.stan

  // SDs
  vector<lower=0>[Kd2]  s_rr;  // round-robin residuals
  vector<lower=0>[allKp] S_p;  // person-level covariates + AP effects

  // correlations among observed variables and random effects
  cholesky_factor_corr[allKp] chol_r_p; // person-level random-effect correlations
  // correlations among round-robin residuals
  real<lower=0,upper=1> r_d2[Kd2];          // dyadic reciprocity (within variable)
  real<lower=0,upper=1> r_intra[nPairs];    // intrapersonal residuals (between variable)
  real<lower=0,upper=1> r_inter[nPairs];    // interpersonal residuals (between variable)

  // random effects to sample on unit scale
  matrix[Np, 2*Kd2] AP; // matrix of actor-partner effects for all RR variables

#include /multigroup/group_parameters.stan

}
transformed parameters {
  // expected values, given random effects
  matrix[Nd, 2*Kd2] Yd2hat;   // dyad-level \hat{y}s
  // combined dyad-level SDs and correlations
  vector[2*Kd2] S_d;
  matrix[2*Kd2, 2*Kd2] Rd2;
  // cholesky decompositions of correlation matrices
  matrix[2*Kd2, 2*Kd2] chol_d;  // dyad-level
  matrix[allKp, allKp] chol_p;  // case-level


  // combine correlations among round-robin variables
  {
    int idx1;  // arbitrary iterators
    int idx2;
    int idp1;
    int idp2;
    int pair;

    pair = 1;
    for (k in 1:Kd2) {
      idx1 = k*2 - 1;
      idx2 = k*2;

      // within round-robin variable
      S_d[idx1] = s_rr[k];        // equal relationship SDs
      S_d[idx2] = s_rr[k];
      Rd2[idx1, idx1] = 1;        // diagonal = 1
      Rd2[idx2, idx2] = 1;
      Rd2[idx1, idx2] = -1 + 2*r_d2[k];  // equal dyadic reciprocity
      Rd2[idx2, idx1] = -1 + 2*r_d2[k];

      // bewteen round-robin variables
      if (k < Kd2) for (kk in (k+1):Kd2) {
        idp1 = kk*2 - 1;
        idp2 = kk*2;

        Rd2[idx1, idp1] = -1 + 2*r_intra[pair]; // within person
        Rd2[idx2, idp2] = -1 + 2*r_intra[pair];
        Rd2[idp1, idx1] = -1 + 2*r_intra[pair];
        Rd2[idp2, idx2] = -1 + 2*r_intra[pair];
        Rd2[idx1, idp2] = -1 + 2*r_inter[pair]; // between person
        Rd2[idp2, idx1] = -1 + 2*r_inter[pair];
        Rd2[idx2, idp1] = -1 + 2*r_inter[pair];
        Rd2[idp1, idx2] = -1 + 2*r_inter[pair];
        pair += 1;
      }

    }

    // end block combining dyad-level correlation matrix
  }
  // cholesky decompositions for model{} block
  chol_d = diag_pre_multiply(S_d, cholesky_decompose(Rd2));

  // scale cholesky for case-level covariates?
  // TODO: Use #include /vanilla/tpar_chol_p.stan
  chol_p = chol_r_p;

  // calculate and/or combine expected values
  {
    int idx1;  // arbitrary iterators
    int idx2;

    for (k in 1:Kd2) {
      idx1 = k*2 - 1; // actor effect for k^th measure
      idx2 = k*2;   // partner effect for k^th measure

      for (d in 1:Nd) {
        // expected values of round-robin variables, given random effects
        Yd2hat[d, idx1] = S_p[idx1]*AP[ IDp[d,1], idx1] + S_p[idx2]*AP[ IDp[d,2], idx2];
        Yd2hat[d, idx2] = S_p[idx1]*AP[ IDp[d,2], idx1] + S_p[idx2]*AP[ IDp[d,1], idx2];

#include /means/means_tparameters.stan
#include /multigroup/group_tparameters.stan
      }

    }
  }

}
model {
  // priors for means and SDs, based on empirical ranges
  for (k in 1:Kd2) {
    // means
#include /means/means_model.stan
    // residual/dyadic SDs
    s_rr[k]      ~ student_t(4, halfRangeKd2[k] / 5, halfRangeKd2[k] / 5) T[0, ];
    // actor effect SDs
    S_p[2*k - 1] ~ student_t(4, halfRangeKd2[k] / 5, halfRangeKd2[k] / 5) T[0, ];
    // partner effect SDs
    S_p[2*k]     ~ student_t(4, halfRangeKd2[k] / 5, halfRangeKd2[k] / 5) T[0, ];
  }

  // priors for correlations
  chol_r_p ~ lkj_corr_cholesky(2);
  for (k in 1:Kd2) r_d2[k] ~ beta(1.5, 1.5);
  if (Kd2 > 1) for (kk in 1:nPairs) {
    r_intra[kk] ~ beta(1.5, 1.5);
    r_inter[kk] ~ beta(1.5, 1.5);
  }

  // likelihoods for observed data (== priors for imputed data & random effects)
  for (n in 1:Nd) Yd2[n,] ~ multi_normal_cholesky(Yd2hat[n,], chol_d);
  for (n in 1:Np) AP[n,] ~ multi_normal_cholesky(rep_row_vector(0, 2*Kd2), chol_p);

  // priors for group effects + their SDs & correlations
#include /multigroup/group_model.stan
}
generated quantities{
  matrix[Nd, 2*Kd2] Yd2e;   // residuals (relationship effects + error)
  matrix[allKp, allKp] Rp;  // person-level correlation matrix

#include /multigroup/group_GQs.stan

  // calculate residuals to return as relationship effects
  Yd2e = Yd2 - Yd2hat;

  // calculate person-level correlation matrix
  Rp = multiply_lower_tri_self_transpose(chol_r_p);
}

#include /include/license.stan
