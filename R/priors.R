### Terrence D. Jorgensen
### Last updated: 7 June 2023
### function to set default priors for mvsrm()

## - t_df, t_m, t_sd: matrix[Kd2, 3] for 3 RR-components, vector[K] for covariate SDs
## - lkj_p(d) prior parameter for cor_p (cor_g when relevant)
## - beta_a, beta_b: matrix[Kd2,Kd2] for dyadic reciprocity (diagonal),
##                   intra/inter correlations above/below diagonal

## rr.data must ALWAYS be long-uni format

srm_priors <- function(rr.data, cov_d, cov_p, cov_g,
                       modelG = FALSE, modelM = FALSE) {
  priors <- list()

  ## SDs for each RR component (out, in, rel)
  ## by default, half-range / 5 for each (so divide range by 10)
  halfRrr <- sapply(rr.data, function(v) diff(range(v, na.rm = TRUE)) / 2)
  priors$rr_rel_t <- data.frame(df = 4, m = halfRrr/5, sd = halfRrr/5)
  priors$rr_out_t <- data.frame(df = 4, m = halfRrr/5, sd = halfRrr/5)
  priors$rr_in_t  <- data.frame(df = 4, m = halfRrr/5, sd = halfRrr/5)
  if (modelG) priors$rr_group_t <- data.frame(df = 4, m = halfRrr/5, sd = halfRrr/5)

  if (!missing(cov_d)) {
    halfRd <- sapply(cov_d, function(v) diff(range(v, na.rm = TRUE)) / 2)
    priors$d1_dyad_t <- data.frame(df = 4, m = halfRd/5, sd = halfRd/2)
    priors$d1_case_t <- data.frame(df = 4, m = halfRd/5, sd = halfRd/2)
    if (modelG) priors$d1_group_t <- data.frame(df = 4, m = halfRd/5, sd = halfRd/5)
  }

  if (!missing(cov_p)) {
    halfRc <- sapply(cov_p, function(v) diff(range(v, na.rm = TRUE)) / 2)
    priors$case_cov_t <- data.frame(df = 4, m = halfRc/3, sd = halfRc/3)
    if (modelG) priors$case_group_t <- data.frame(df = 4, m = halfRc/5, sd = halfRc/5)
  }

  if (modelG & !missing(cov_g)) {
    halfRg <- sapply(cov_g, function(v) diff(range(v, na.rm = TRUE)) / 2)
    priors$group_cov_t <- data.frame(df = 4, m = halfRg/3, sd = halfRg/3)
  }

  ## LKJ priors for intact correlation matrices
  priors$case_lkj <- 2
  if (modelG) priors$group_lkj <- 2
  if (!missing(cov_d)) {
    if (ncol(cov_d) > 1L) priors$d1_lkj <- 2 # correlations among cov_d
  }

  ## beta priors for constrained dyad-level correlations
  priors$rr_beta_a <- matrix(1.5, nrow = length(rr.data), ncol = length(rr.data),
                             dimnames = list(names(rr.data), names(rr.data)))
  priors$rr_beta_b <- priors$rr_beta_a
  if (!missing(cov_d)) {
    ## correlations between cov_d and rr.data
    priors$d12_beta_a <- matrix(1.5, nrow = length(cov_d), ncol = length(rr.data),
                                dimnames = list(names(cov_d), names(rr.data)))
    priors$d12_beta_b <- priors$d12_beta_a
  }

  if (modelM) {
    priors$rr_Mvec_m  <- sapply(rr.data, median, na.rm = TRUE)
    priors$rr_Mvec_sd <- halfRrr
    #TODO: separate Mvec (RR) from group-means of each level's covariate(s)?
  }

  priors
}
