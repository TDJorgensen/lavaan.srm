// Terrence D. Jorgensen
// Last updated: 6 June 2023

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

  // (hyper)priors for model{}    (group-level priors #included above)
#include /means/means_data.stan
  // df,m,sd for SDs of round-robin variable components
  matrix[Kd2, 3] rr_rel_t;  // dyad (relationship)
  matrix[Kd2, 3] rr_out_t;  // out-going (actor / perceiver)
  matrix[Kd2, 3] rr_in_t ;  // in-coming (partner / target)
  // LKJ parameter for case-level correlations
  real case_lkj;
  // beta parameters (a and b) for dyad-level correlations
  matrix[Kd2, Kd2] rr_beta_a;   // (dyadic reciprocity on diagonal,
  matrix[Kd2, Kd2] rr_beta_b;   // intra/inter correlations above/below diagonal)
}
transformed data {
#include /vanilla/tdata_allKd.stan
#include /vanilla/tdata_allKp.stan
#include /vanilla/tdata_allKg.stan
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
  //  - dyadic reciprocity (within variable) on diagonal
  //  - intrapersonal correlations (between variable, within  case) below diagonal
  //  - interpersonal correlations (between variable, between case) above diagonal
  matrix<lower=0,upper=1>[Kd2, Kd2] r_d2;

  // random effects to sample on unit scale
  matrix[Np, 2*Kd2] AP; // matrix of actor-partner effects for all RR variables

#include /multigroup/group_parameters.stan

}
transformed parameters {
  // expected values, given random effects
  matrix[Nd, allKd] Yd2hat;   // dyad-level \hat{y}s
  // combined dyad-level SDs and correlations
  vector[allKd] S_d;
  matrix[allKd, allKd] Rd2;
  // cholesky decompositions of correlation matrices
  matrix[allKd, allKd] chol_d;  // dyad-level
  matrix[allKp, allKp] chol_p;  // case-level
#include /multigroup/group_tpar_declare.stan

  // augmented level-specific data sets
  // (can include observed covariates, random effects, imputed missing values)
  matrix[Nd, allKd] augYd;  // all dyad-level variables
  matrix[Np, allKp] augYp;  // all case-level variables

  // combine correlations among round-robin variables
  {
    int idx1;  // arbitrary iterators
    int idx2;
    int idp1;
    int idp2;

    for (k in 1:Kd2) {
      idx1 = k*2 - 1;
      idx2 = k*2;

      // within round-robin variable
      S_d[idx1] = s_rr[k];        // equal relationship SDs
      S_d[idx2] = s_rr[k];
      Rd2[idx1, idx1] = 1;        // diagonal = 1
      Rd2[idx2, idx2] = 1;
      Rd2[idx1, idx2] = -1 + 2*r_d2[k,k];  // equal dyadic reciprocity
      Rd2[idx2, idx1] = -1 + 2*r_d2[k,k];

      // bewteen round-robin variables
      if (k < Kd2) { for (kk in (k+1):Kd2) {
        idp1 = kk*2 - 1;
        idp2 = kk*2;

        Rd2[idx1, idp1] = -1 + 2*r_d2[kk, k ]; // within person (intra = BELOW)
        Rd2[idx2, idp2] = -1 + 2*r_d2[kk, k ];
        Rd2[idp1, idx1] = -1 + 2*r_d2[kk, k ];
        Rd2[idp2, idx2] = -1 + 2*r_d2[kk, k ];
        Rd2[idx1, idp2] = -1 + 2*r_d2[k , kk]; // between person (inter = ABOVE)
        Rd2[idp2, idx1] = -1 + 2*r_d2[k , kk];
        Rd2[idx2, idp1] = -1 + 2*r_d2[k , kk];
        Rd2[idp1, idx2] = -1 + 2*r_d2[k , kk];
      }}

    }

    // end block combining dyad-level correlation matrix
  }
  // cholesky decompositions for model{} block
  chol_d = diag_pre_multiply(S_d, cholesky_decompose(Rd2));

#include /vanilla/tpar_chol_p.stan
#include /vanilla/tpar_chol_g.stan

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

  // Augment observed values with ...
  // - level-specific covariates?
  // - estimates of missing values?
#include /vanilla/tpar_augYd.stan
#include /vanilla/tpar_augYp.stan

}
model {
  // priors for means and SDs, based on empirical ranges
  for (k in 1:Kd2) {
    // means
#include /means/means_model.stan
    // residual/dyadic SDs
    s_rr[k]      ~ student_t(rr_rel_t[k,1], rr_rel_t[k,2], rr_rel_t[k,3]) T[0, ];
    // actor effect SDs
    S_p[2*k - 1] ~ student_t(rr_out_t[k,1], rr_out_t[k,2], rr_out_t[k,3]) T[0, ];
    // partner effect SDs
    S_p[2*k]     ~ student_t(rr_in_t[k,1 ], rr_in_t[k, 2], rr_in_t[k, 3]) T[0, ];
  }

  // priors for correlations
  chol_r_p ~ lkj_corr_cholesky(case_lkj);
  for (k in 1:Kd2) {
    // dyadic correlations (priors on diagonal)
    r_d2[k,k] ~ beta(rr_beta_a[k,k], rr_beta_b[k,k]);
    // between-variable correlations
    if (k < Kd2) { for (kk in (k+1):Kd2) {
      r_d2[kk, k ] ~ beta(rr_beta_a[kk, k ], rr_beta_b[kk, k ]); // intra = BELOW
      r_d2[k , kk] ~ beta(rr_beta_a[k , kk], rr_beta_b[k , kk]); // inter = ABOVE
    }}
  }

  // likelihoods for observed data (== priors for imputed data & random effects)
  for (n in 1:Nd) augYd[n,] ~ multi_normal_cholesky(Yd2hat[n,], chol_d);
  for (n in 1:Np) augYp[n,] ~ multi_normal_cholesky(rep_row_vector(0, allKp), chol_p);

  // priors for group effects + their SDs & correlations
#include /multigroup/group_model.stan
}
generated quantities{
  matrix[Nd, allKd] Yd2e;      // residuals (relationship effects + error)
  matrix[allKp, allKp] Rp;     // person-level correlation matrix
  matrix[allKp, allKp] pSigma; // person-level covariance matrix
  matrix[allKp, allKp] dSigma; // dyad-level covariance matrix

#include /multigroup/group_GQs.stan

  // calculate residuals to return as relationship effects
  Yd2e = augYd - Yd2hat;

  // calculate person-level correlation matrix
  Rp = multiply_lower_tri_self_transpose(chol_r_p);
  // calculate group-level covariance matrix
  pSigma = multiply_lower_tri_self_transpose(chol_p);
  // calculate dyad-level covariance matrix
  dSigma = multiply_lower_tri_self_transpose(chol_d);
}

#include /include/license.stan
