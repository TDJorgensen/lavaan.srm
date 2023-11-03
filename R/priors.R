### Terrence D. Jorgensen
### Last updated: 3 November 2023
### function to set default priors for mvsrm()


##' Default priors for multivariate SRM
##'
##' Priors used during MCMC estimation of a multivariate social relations model
##' (mvSRM; Nestler, [2018](https://doi.org/10.3102/1076998617741106)) using the
##' Stan (Carpenter et al., [2017](https://doi.org/10.18637/jss.v076.i01)) software.
##'
##' @param data,group_data,case_data See the [mvsrm()] function for details.
##'        However, the data sets should contain **only** modeled variables (no
##'        ID variables), and `data=` must contain round-robin variables in long
##'        format (one column per variable, 2 rows per dyad).
  #TODO: @param cov_d dyad-constant covariates
  #      separate data.frame or names in `data`?  Could infer from rr.vars=
##' @param modelG `logical` indicating whether to estimate group-level (co)variances
##' @param modelM `logical` indicating whether to estimate means. Setting `TRUE`
##'        presupposes either that there is only 1 round-robin group
##'        (`modelG=FALSE`) or that group-level effects will not be partialed
##'        out as fixed effects (`fixed.groups=FALSE` when [mvsrm()] is called).
##' @param SDby `character` indicating how to choose priors for *SD*s of each
##'        round-robin variable component (i.e., dyad- and case-level, and
##'        potentially group-level if `modelG=TRUE`). Options include
##'        `"sd"` to base it on the (total) sample *SD*, or `"range"` to
##'        base it on the empirical minimum and maximum observation. The latter
##'        is likely to overestimate the actual *SD*, and it might only be a
##'        sensible heuristic when using Likert-type data, due to their
##'        limited response options that determine the maximum possible *SD*.
##' @param decomp `numeric` vector that sums to 1 (or will be rescaled to do so)
##'        indicating the hypothesized proportion of total variance attributable
##'        to each round-robin component. The vector should therefore have the
##'        following [base::names()]: `c("dyad","out","in")`, as well as
##'        `"group"` when `modelG=TRUE`. Note that the name `"in"` must be in
##'        backticks when assigning names within the [base::c()] function (see
##'        the default value in **Usage**), in order to avoid confusion with the
##'        [base::Reserved()] word in R.  The `"out"` and `"in"` components
##'        correspond to what are often called "perceivers" and "targets" in
##'        SRMs of interpersonal perception, or "actors" and "partners" in SRMs
##'        of group behavior. The same hypothesized decomposition is applied to
##'        each round-robin variable, unless the user provides a list of such
##'        vectors (with length equal to the number of variables in `data=`) to
##'        apply to each round-robin variable. A list of vectors can have names
##'        that correspond to those in `data=`, but will otherwise be assumed to
##'        follow the same order.
##' @param decomp_c `numeric` vector (or list of vectors) similar to `decomp`
##'        argument, but for the decomposition of `case_data` into case- and
##'        group-level variance components.  Ignored unless `modelG=TRUE`.
##'
##' @references
##'   Carpenter, B., Gelman, A., Hoffman, M. D., Lee, D., Goodrich, B.,
##'   Betancourt, M., Brubaker, M., Guo, J., Li, P., & Riddell, A. (2017).
##'   Stan: A probabilistic programming language.
##'   *Journal of Statistical Software, 76*(1).
##'   <https://doi.org/10.18637/jss.v076.i01>
##'
##'  Nestler, S. (2018). Likelihood estimation of the multivariate social
##'  relations model. *Journal of Educational and Behavioral Statistics, 43*(4),
##'  387--406. <https://doi.org/10.3102/1076998617741106>
##'
##' @examples
##' ## use example data from the ?srm::srm help page
##' data(data.srm01, package="srm")
##'
##' ## ignoring group-level variance, not estimating means
##' ## (i.e., mvsrm() called with fixed.groups = TRUE)
##' srm_priors(data.srm01[5:7])
##'
##' ## include group-level variance and means (fixed.groups = FALSE)
##' srm_priors(data.srm01[5:7], modelG = TRUE, modelM = TRUE)
##'
##' ## base SD hyperparameters on observed range instead of observed SD
##' srm_priors(data.srm01[5:7], SDby = "range")
##'
##' @importFrom stats sd
##' @export
srm_priors <- function(data, group_data, case_data, # cov_d or rr.vars = NULL,
                       modelG = FALSE, modelM = FALSE, SDby = "sd",
                       decomp = c(dyad = .5, out = .2, `in` = .2, group = .1),
                       decomp_c = c(case = .9, group = .1)) {
  ## Priors for Means:
  ## - normal priors, use M(edian) and SD of each RR variable
  ## TODO: also for covariate means

  ## Priors for SDs:
  ## - t_df, t_m, t_sd: matrix[Kd2, 3] for 3 RR-components
  ## TODO: also for covariate SDs

  ## Priors for group/case-level correlations:
  ## - lkj_p(d) prior parameter for cor_p (cor_g when relevant)
  ## - beta_a, beta_b: matrix[Kd2,Kd2] for dyadic reciprocity (diagonal),
  ##                   intra/inter correlations above/below diagonal
  ## TODO: Enable Wishart prior: df and a target matrix.
  ##       Use cor(data) for target?  Need to expand to out/in components.
  ##       Add argument to specify expected generalized reciprocities, set
  ##       other out--in correlations (between variables) to 0?

  #TODO: When Stan model is ready, separate cov_d from data[rr.vars]
  rr.data <- data#[rr.vars]
  # cov_d <- data[ setdiff(names(data), rr.vars) ]
  cov_d <- NULL

  if (inherits(decomp, "list")) {
    ## as many decompositions as variables?
    stopifnot(length(decomp) == ncol(rr.data))
    ## all decompositions have the same length?
    stopifnot(length(unique(sapply(decomp, length)) > 1L))

    ## should all have the same names
    if (is.null(names(decomp))) names(decomp) <- names(rr.data)
    ## convert to a matrix
    decomp <- sapply(decomp, function(v) v[c("dyad","out","in","group")])

  } else {

    stopifnot(inherits(decomp, "numeric"))
    stopifnot(length(decomp) >= 3L + modelG)
    if (is.null(names(decomp))) {
      names(decomp)[1:3] <- c("dyad","out","in")
      if (modelG) names(decomp)[4] <- "group"
    } else {
      stopifnot(all(names(decomp) %in% c("dyad","out","in","group")))
    }
    ## repeat values (in a matrix) per RR variable
    decomp <- sapply(names(rr.data), function(v) decomp[c("dyad","out","in","group")])
  }

  ## will case-level covariates be decomposed?
  if (modelG & !missing(case_data)) {

    if (inherits(decomp_c, "list")) {
      ## as many decompositions as variables?
      stopifnot(length(decomp_c) == ncol(case_data))
      ## all decompositions have the same length?
      stopifnot(length(unique(sapply(decomp_c, length)) > 1L))

      ## should all have the same names
      if (is.null(names(decomp_c))) names(decomp_c) <- names(case_data)
      ## convert to a matrix
      decomp_c <- sapply(decomp_c, function(v) v[c("case","group")])

    } else {

      stopifnot(inherits(decomp_c, "numeric"))
      stopifnot(length(decomp_c) == 2L)
      if (is.null(names(decomp_c))) {
        names(decomp_c) <- c("case","group")
      } else {
        stopifnot(all(names(decomp_c) %in% c("case","group")))
      }
      ## repeat values (in a matrix) per RR variable
      decomp_c <- sapply(names(case_data), function(v) decomp_c[c("case","group")])
    }
  }
  #TODO: also for dyadic and group-level covariates, when implemented


  priors <- list()

  if (SDby == "sd") {
    #FIXME: only applicable for continuous variables
    #TODO:  theta parameterization (fix dyad SD=1). How to set case/group SDs?
    rrSD <- sapply(rr.data, sd, na.rm = TRUE)

    if (         !is.null(cov_d)     ) dSD <- sapply(cov_d     , sd, na.rm = TRUE)
    if (         !missing(case_data) ) cSD <- sapply(case_data , sd, na.rm = TRUE)
    if (modelG & !missing(group_data)) gSD <- sapply(group_data, sd, na.rm = TRUE)

  } else if (SDby == "range") {
    # each "propRange" is heuristic.  Make it an argument?

    maxPossibleSD <- sapply(rr.data, function(v) diff(range(v, na.rm = TRUE)) / 2)
    propRange <- 1 / 2
    rrSD <- maxPossibleSD * propRange

    if (!is.null(cov_d)) {
      maxPossibleSD_d <- sapply(cov_d, function(v) diff(range(v, na.rm = TRUE)) / 2)
      propRange_d <- 1 / 2
      dSD <- maxPossibleSD_d * maxPossibleSD_d
    }
    if (!missing(case_data)) {
      maxPossibleSD_c <- sapply(case_data, function(v) diff(range(v, na.rm = TRUE)) / 2)
      propRange_c <- 1 / 2
      cSD <- maxPossibleSD_c * maxPossibleSD_c
    }
    if (modelG & !missing(group_data)) {
      maxPossibleSD_g <- sapply(group_data, function(v) diff(range(v, na.rm = TRUE)) / 2)
      propRange_g <- 1 / 2
      gSD <- maxPossibleSD_g * maxPossibleSD_g
    }

  } else stop("Invalid choice for 'SDby=' argument")

  ## SDs for each RR component (out, in, rel)
  priors$rr_rel_t <- priors$rr_out_t <- priors$rr_in_t <- data.frame(df = rep(4, ncol(rr.data)))

  priors$rr_rel_t$sd <- priors$rr_rel_t$m <- sqrt(rrSD^2 * decomp["dyad",])
  priors$rr_out_t$sd <- priors$rr_out_t$m <- sqrt(rrSD^2 * decomp["out" ,])
  priors$rr_in_t$sd  <- priors$rr_in_t$m  <- sqrt(rrSD^2 * decomp["in"  ,])
  if (modelG) {
    priors$rr_group_t <- data.frame(df = rep(4, length(rrSD)))
    priors$rr_group_t$sd <- priors$rr_group_t$m <- sqrt(rrSD^2 * decomp["group",])
  }

  if (!is.null(cov_d)) {
    ## heuristic decomposition below assumes negligible case/group-level variance
    priors$d1_dyad_t <- data.frame(df = 4,
                                   m  = sqrt(dSD^2 * .8),
                                   sd = sqrt(dSD^2 * .8))
    priors$d1_case_t <- data.frame(df = 4,
                                   m  = sqrt(dSD^2 * .1),
                                   sd = sqrt(dSD^2 * .1))
    if (modelG) priors$d1_group_t <- data.frame(df = 4,
                                                m  = sqrt(dSD^2 * .1),
                                                sd = sqrt(dSD^2 * .1))
  }

  if (!missing(case_data)) {
    if (modelG) {
      priors$case_cov_t <- data.frame(df = 4,
                                      m  = sqrt(cSD^2 * decomp_c["case"]),
                                      sd = sqrt(cSD^2 * decomp_c["case"]))
      priors$case_group_t <- data.frame(df = 4,
                                        m  = sqrt(cSD^2 * decomp_c["group"]),
                                        sd = sqrt(cSD^2 * decomp_c["group"]))
    } else priors$case_cov_t <- data.frame(df = 4, m  = cSD, sd = cSD)
  }

  if (modelG & !missing(group_data)) {
    priors$group_cov_t <- data.frame(df = 4, m = gSD, sd = gSD)
  }

  ## LKJ priors for intact correlation matrices
  priors$case_lkj <- 2
  if (modelG) priors$group_lkj <- 2
  if (!is.null(cov_d)) {
    if (ncol(cov_d) > 1L) priors$d1_lkj <- 2 # correlations among cov_d
  }

  ## beta priors for constrained dyad-level correlations
  priors$rr_beta_a <- matrix(1.5, nrow = length(rr.data), ncol = length(rr.data),
                             dimnames = list(names(rr.data), names(rr.data)))
  priors$rr_beta_b <- priors$rr_beta_a
  if (!is.null(cov_d)) {
    ## correlations between cov_d and rr.data
    priors$d12_beta_a <- matrix(1.5, nrow = length(cov_d), ncol = length(rr.data),
                                dimnames = list(names(cov_d), names(rr.data)))
    priors$d12_beta_b <- priors$d12_beta_a
  }

  if (modelM) {
    priors$rr_Mvec_m  <- sapply(rr.data, median, na.rm = TRUE)
    priors$rr_Mvec_sd <- rrSD
    #TODO: separate Mvec (RR) from group-means of each level's covariate(s)?
  }

  priors
}

# Set uninformative Wishart prior as identity with df = dim + 2?
# hist(as.numeric(apply(rWishart(n = 1000, df = 5, Sigma = diag(3)), MARGIN = 3,
#                       function(x) cov2cor(x)[lower.tri(diag(3), diag = FALSE)])))


