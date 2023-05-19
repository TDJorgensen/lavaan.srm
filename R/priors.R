### Terrence D. Jorgensen
### Last updated: 20 May 2023
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
  ## by default, half-range / 5 for each
  halfR5rr <- sapply(rr.data, function(v) diff(range(v, na.rm = TRUE)) / 5)
  priors$rr_rel_t <- data.frame(df = 4, m = halfR5rr, sd = halfR5rr)
  priors$rr_out_t <- data.frame(df = 4, m = halfR5rr, sd = halfR5rr)
  priors$rr_in_t  <- data.frame(df = 4, m = halfR5rr, sd = halfR5rr)
  if (modelG) priors$rr_group_t <- data.frame(df = 4, m = halfR5rr, sd = halfR5rr)

  if (!missing(cov_d)) {
    halfR5d <- sapply(cov_d, function(v) diff(range(v, na.rm = TRUE)) / 5)
    priors$d1_dyad_t <- data.frame(df = 4, m = halfR5d, sd = halfR5d)
    priors$d1_case_t <- data.frame(df = 4, m = halfR5d, sd = halfR5d)
    if (modelG) priors$d1_group_t <- data.frame(df = 4, m = halfR5d, sd = halfR5d)
  }

  if (!missing(cov_p)) {
    halfR5c <- sapply(cov_p, function(v) diff(range(v, na.rm = TRUE)) / 5)
    priors$case_cov_t <- data.frame(df = 4, m = halfR5c, sd = halfR5c)
    if (modelG) priors$case_group_t <- data.frame(df = 4, m = halfR5c, sd = halfR5c)
  }

  if (modelG & !missing(cov_g)) {
    halfR5g <- sapply(cov_g, function(v) diff(range(v, na.rm = TRUE)) / 5)
    priors$group_cov_t <- data.frame(df = 4, m = halfR5g, sd = halfR5g)
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
    priors$rr_Mvec_sd <- halfR5rr * 5 # use half-range, not divided by 5
    #TODO: separate Mvec (RR) from group-means of each level's covariate(s)?
  }

  priors
}
