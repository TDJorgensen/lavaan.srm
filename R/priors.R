### Terrence D. Jorgensen
### Last updated: 8 November 2023
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
##' @param decomp_rr `numeric` vector of weights (internally rescaled so sum ==
##'        1) indicating hypothesized proportions of total variance attributable
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
##' @param decomp_c `numeric` vector (or list of vectors) similar to `decomp_rr`
##'        argument, but for the decomposition of `case_data` into case- and
##'        group-level variance components.  Ignored unless `modelG=TRUE`.
##'
##' @return A `list` of hyperparameters, with an `attr(,"Nvars")` indicating
##'         the number of (components of) variables modeled at each level (e.g.,
##'         the `case` level includes 2 components per round-robing variable
##'         plus the number of variables in `case_data=`).
##'
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
##' data(data.srm01, package = "srm")
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
                       decomp_rr = c(dyad = .5, out = .2, `in` = .2, group = .1),
                       decomp_c = c(case = .9, group = .1)) {
  ## Priors for Means:
  ## - normal priors, use M(edian) and SD of each RR variable
  ## TODO: also for covariate means

  ## Priors for SDs:
  ## - t_df, t_m, t_sd: matrix[Kd2, 3] for 3 RR-components
  ## TODO: also for covariate SDs

  ## Priors for dyad-level correlations:
  ## - rrD_beta_a, rrD_beta_b: matrix[Kd2,Kd2] for dyadic reciprocity (on diag),
  ##                           intra/inter correlations (above/below the diag)
  ## - d1_beta: matrix[Kd1,Kd1] with shapes a (above) and b (below) the diag
  ## - d21_beta_a, d21_beta_b: matrix[Kd2,Kd1] equality-constrained correlations
  ##                           between each RR and each symmetric (d1) variable
  ## TODO: Switch to beta priors for case and group levels
  ## -  case_beta: matrix[Kc,Kc] with shapes a (above) and b (below) the diag
  ## - group_beta: matrix[Kg,Kg] with shapes a (above) and b (below) the diag

  #TODO: When Stan model is ready, separate cov_d from data[rr.vars]
  rr.data <- data#[rr.vars]
  # cov_d <- data[ setdiff(names(data), rr.vars) ]
  cov_d <- NULL

  ## Count number of variables/components at each level
  ## (used to define dimensions for correlation priors)
  Nvars <- list(rr = ncol(rr.data), d1 = max(ncol(cov_d), 0),
                case = 2*ncol(rr.data) + max(ncol(cov_d), 0))
  if (!missing(case_data)) Nvars$case  <- Nvars$case + ncol(case_data)
  if (modelG) {
    Nvars$group <- ncol(rr.data) + max(ncol(cov_d), 0)
    if (!missing( case_data)) Nvars$group <- Nvars$group + ncol(case_data)
    if (!missing(group_data)) Nvars$group <- Nvars$group + ncol(group_data)
  }


  ## decomposition of round-robin variables
  if (inherits(decomp_rr, "list")) {
    ## as many decompositions as variables?
    stopifnot(length(decomp_rr) == ncol(rr.data))
    ## all decompositions have the same length?
    stopifnot(length(unique(sapply(decomp_rr, length)) > 1L))

    ## should all have the same names
    if (is.null(names(decomp_rr))) names(decomp_rr) <- names(rr.data)
    ## convert to a matrix
    decomp_rr <- sapply(decomp_rr, function(v) {
      components <- c("dyad","out","in")
      if (modelG) components <- c(components, "group")
      V <- v[components]
      V / sum(V) # rescale so sum(V) == 1
    })

  } else {

    stopifnot(inherits(decomp_rr, "numeric"))
    stopifnot(length(decomp_rr) >= 3L + modelG)
    if (is.null(names(decomp_rr))) {
      names(decomp_rr)[1:3] <- c("dyad","out","in")
      if (modelG) names(decomp_rr)[4] <- "group"
    } else {
      stopifnot(all(names(decomp_rr) %in% c("dyad","out","in","group")))
    }
    ## repeat values (in a matrix) per RR variable
    decomp_rr <- sapply(names(rr.data), function(v) {
      components <- c("dyad","out","in")
      if (modelG) components <- c(components, "group")
      V <- decomp_rr[components]
      V / sum(V) # rescale so sum(V) == 1
    })
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
      decomp_c <- sapply(decomp_c, function(v) {
        V <- v[c("case","group")]
        V / sum(V) # rescale so sum(V) == 1
      })

    } else {

      stopifnot(inherits(decomp_c, "numeric"))
      stopifnot(length(decomp_c) == 2L)
      if (is.null(names(decomp_c))) {
        names(decomp_c) <- c("case","group")
      } else {
        stopifnot(all(names(decomp_c) %in% c("case","group")))
      }
      ## repeat values (in a matrix) per RR variable
      decomp_c <- sapply(names(case_data), function(v) {
        V <- decomp_c[c("case","group")]
        V / sum(V) # rescale so sum(V) == 1
      })
    }
  }
  #TODO: also for dyadic and group-level covariates, when implemented


  priors <- list()
  attr(priors, "Nvars") <- Nvars # for mvsrm() to pass as data= to stan()

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

  ## STANDARD DEVIATIONS

  ## SDs for each RR component (out, in, rel)
  priors$rr_rel_t <- priors$rr_out_t <- priors$rr_in_t <- data.frame(df = rep(4, ncol(rr.data)))

  priors$rr_rel_t$sd <- priors$rr_rel_t$m <- sqrt(rrSD^2 * decomp_rr["dyad",])
  priors$rr_out_t$sd <- priors$rr_out_t$m <- sqrt(rrSD^2 * decomp_rr["out" ,])
  priors$rr_in_t$sd  <- priors$rr_in_t$m  <- sqrt(rrSD^2 * decomp_rr["in"  ,])
  if (modelG) {
    priors$rr_group_t <- data.frame(df = rep(4, length(rrSD)))
    priors$rr_group_t$sd <- priors$rr_group_t$m <- sqrt(rrSD^2 * decomp_rr["group",])
  }

  ## SDs for covariates at each level

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


  ## CORRELATIONS using beta priors

  priors$case_lkj <- 2
  if (modelG) priors$group_lkj <- 2
  if (!is.null(cov_d)) {
    if (ncol(cov_d) > 1L) priors$d1_lkj <- 2 # correlations among cov_d
  }

  ## equality-constrained dyad-level correlations
  priors$rrD_beta_a <- matrix(1.5, nrow = length(rr.data), ncol = length(rr.data),
                              dimnames = list(names(rr.data), names(rr.data)))
  priors$rrD_beta_b <- priors$rrD_beta_a
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


