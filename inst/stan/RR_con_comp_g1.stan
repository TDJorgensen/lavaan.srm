// Terrence D. Jorgensen
// Last updated: 13 November 2023

// Program to estimate covariance matrices for multivariate SRM components
// - standard round-robin design (indistinguishable subjects)
//    - only 1 RR group (means still estimated)
// - continuous data only
// - no missing data


functions {}
data {
  // sample sizes
  int<lower=0> Nd;    // number of dyads   (Level 1)
  int<lower=0> Np;    // number of persons (Level 2, cross-classified)

  // number of observed measures (half the number of columns)
  int<lower=0> Kd2;   // number of dyad-level    measures (    vary within dyad)
  int<lower=0> Kd1;   // number of dyad-level  covariates (constant within dyad)
  int<lower=0> Kp ;   // number of case-level  covariates
  int<lower=0> Kg ;   // number of group-level covariates
  // number of modeled components at each level
  int<lower=0> allKd; // number of  dyad-level components (ij, ji, and d1)
  int<lower=0> allKp; // number of  case-level components (of d2, d1, and c)
  int<lower=0> allKg; // number of group-level components (of d2, d1, c and g)

  // observed data
  matrix[Nd, 2*Kd2] Yd2;  // observed round-robin variables
  // ID variables
  int IDp[Nd, 2];         // person-level IDs (cross-classified) for dyad-level observations

  // (hyper)priors for model{}

#include /means/means_data.stan

  // t(df,m,sd) for SDs of round-robin variable components
  matrix[Kd2, 3] rr_rel_t;  // dyad (relationship)
  matrix[Kd2, 3] rr_out_t;  // out-going (actor / perceiver)
  matrix[Kd2, 3] rr_in_t ;  // in-coming (partner / target)

  // beta parameters (a and b) for dyad-level correlations
  matrix<lower=0>[Kd2, Kd2] rr_beta_a;   // (dyadic reciprocity on diagonal,
  matrix<lower=0>[Kd2, Kd2] rr_beta_b;   // intra/inter correlations below/above diagonal)

  // beta hyperparameters for case-level correlations
  // alpha below diagonal, beta above diagonal (ignore diagonal)
  matrix<lower=0>[allKp, allKp] case_beta;
}
transformed data {}
parameters {
  // means
#include /means/means_parameters.stan

  // SDs
  vector<lower=0>[Kd2]  s_rr;  // round-robin residuals
  vector<lower=0>[allKp] S_p;  // person-level AP effects (+ covariates)

  // correlations among person-level random effects (+ covariates)
  vector<lower=0,upper=1>[ ( allKp*(allKp-1) ) %/% 2] Rp_vec;

  // correlations among round-robin residuals
  //  - dyadic reciprocity (within variable) on diagonal
  //  - intrapersonal correlations (between variable, within  case) below diagonal
  //  - interpersonal correlations (between variable, between case) above diagonal
  matrix<lower=0,upper=1>[Kd2, Kd2] r_d2;

  // random effects to sample on unit scale
  matrix[Np, 2*Kd2] AP; // matrix of actor-partner effects for all RR variables
}
transformed parameters {
  // expected values, given random effects
  matrix[Nd, allKd] Yd2hat;   // dyad-level \hat{y}s
  // assembled dyad-level SDs and correlations
  vector[allKd] S_d;
  matrix[allKd, allKd] Rd2;
  // assembled case-level correlations
  matrix[allKp, allKp] Rp;
  // cholesky decompositions of correlation matrices
  matrix[allKd, allKd] chol_d;  // dyad-level
  matrix[allKp, allKp] chol_p;  // case-level

  // augmented level-specific data sets
  // (can include observed covariates, random effects, imputed missing values)
  matrix[Nd, allKd] augYd;  // all dyad-level variables
  matrix[Np, allKp] augYp;  // all case-level variables

  // assemble correlations among DYAD-level components of round-robin variables
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

  // assemble correlations among CASE-level components of round-robin variables
  {
    int vec_idx = 1;    // index Rp_vec (vector of correlations)

    for (k in 1:allKp) {
      Rp[k,k] = 1;  // set diagonal = 1

      // column-wise along lower triangle
      // == row-wise along upper triangle
      if (k < allKp) { for (kk in (k+1):allKp) {
        Rp[kk,  k] = -1 + 2*Rp_vec[vec_idx];
        Rp[ k, kk] = -1 + 2*Rp_vec[vec_idx];
        vec_idx += 1;
      }}
    }

    // end block combining case-level correlation matrix
  }

#include /vanilla/tpar_chol_p.stan

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

  // priors for dyad-level correlations
  for (k in 1:Kd2) {
    // within each round-robin variable (dyadic reciprocity)
    r_d2[k,k] ~ beta(rr_beta_a[k,k], rr_beta_b[k,k]);  // priors on diagonal
    // between each round-robin variable (intra/inter-personal)
    if (k < Kd2) { for (kk in (k+1):Kd2) {
      // dyad level
      r_d2[kk, k ] ~ beta(rr_beta_a[kk, k ], rr_beta_b[kk, k ]); // intra = BELOW
      r_d2[k , kk] ~ beta(rr_beta_a[k , kk], rr_beta_b[k , kk]); // inter = ABOVE
    }}
  }

  // priors for case-level correlations
  {
    int vec_idx = 1;    // index Rp_vec (vector of correlations)

    for (k in 1:allKp) {
      // column-wise along lower triangle
      // == row-wise along upper triangle
      if (k < allKp) { for (kk in (k+1):allKp) {
        // hyperparameters:       alpha below diag,   beta above diag
        Rp_vec[vec_idx] ~ beta(   case_beta[kk, k],   case_beta[k, kk] );
        vec_idx += 1;
      }}
    }

    // end block for priors of case-level correlations
  }

  // likelihoods for observed data (== priors for imputed data & random effects)
  for (n in 1:Nd) augYd[n,] ~ multi_normal_cholesky(Yd2hat[n,], chol_d);
  for (n in 1:Np) augYp[n,] ~ multi_normal_cholesky(rep_row_vector(0, allKp), chol_p);
}
generated quantities{
  matrix[Nd, allKd] Yd2e;      // residuals (relationship effects + error)
  matrix[allKp, allKp] pSigma; // person-level covariance matrix
  matrix[allKd, allKd] dSigma; // dyad-level covariance matrix

#include /OneGroup/OneGroup_Rsq.stan

  // calculate residuals to return as relationship effects
  Yd2e = augYd - Yd2hat;

  // scale correlation matrices to covariance matrices
  pSigma = quad_form_diag(Rp , S_p);  // case level
  dSigma = quad_form_diag(Rd2, S_d);  // dyad level
}

#include /include/license.stan
