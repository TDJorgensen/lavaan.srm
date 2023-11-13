### Terrence D. Jorgensen
### Last updated: 10 November 2023
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
##'         the `case` level includes 2 components per round-robin variable
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
  ## - rr_beta_a, rr_beta_b: matrix[Kd2,Kd2] for dyadic reciprocity (on diag),
  ##                         intra/inter correlations (above/below the diag)
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
  } else Nvars$group <- 0L


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


  ## begin STANDARD DEVIATIONS

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
    } else priors$case_cov_t <- data.frame(df = 4, m  = cSD, sd = cSD) #TODO: update /covariates/case_tdata & _model
  }

  if (modelG & !missing(group_data)) {
    priors$group_cov_t <- data.frame(df = 4, m = gSD, sd = gSD)
  }

  ## end priors for STANDARD DEVIATIONS


  ## begin CORRELATIONS using beta priors

  # priors$case_lkj <- 2
  # if (modelG) priors$group_lkj <- 2
  # if (!is.null(cov_d)) {
  #   if (ncol(cov_d) > 1L) priors$d1_lkj <- 2 # correlations among cov_d
  # }

  ## DYAD-level correlations
  priors$rr_beta_a <- matrix(1.5, nrow = Nvars$rr, ncol = Nvars$rr,
                             dimnames = list(names(rr.data), names(rr.data)))
  class(priors$rr_beta_a) <- append(class(priors$rr_beta_a),
                                    values = "lavaan.matrix", after = 0)
  priors$rr_beta_b <- priors$rr_beta_a
  ## add description to guide users
  rr_header <- paste('The matrix below contains shape hyperparameters "alpha" for the',
                     'beta-distribution priors of each round-robin variable\'s dyadic',
                     'reciprocity (on the diagonal), as well as for intra-personal',
                     '(lower triangle) and inter-personal (upper triangle) correlations.')
  attr(priors$rr_beta_a, "header") <- rr_header
  attr(priors$rr_beta_b, "header") <- gsub(x = rr_header, pattern = "alpha",
                                           replacement = "beta")
  if (!is.null(cov_d)) {
    ## correlations between cov_d and rr.data
    priors$d12_beta_a <- matrix(1.5, nrow = Nvars$d1, ncol = Nvars$rr,
                                dimnames = list(names(cov_d), names(rr.data)))
    class(priors$d12_beta_a) <- c("lavaan.matrix", class(priors$d12_beta_a))
    priors$d12_beta_b <- priors$d12_beta_a
    d12_header <- paste('The matrix below contains shape hyperparameters "alpha"',
                        'for the beta-distribution priors of correlations between',
                        'each dyad-level covariate (see rownames) and each',
                        'round-robin variable (see colnames).')
    attr(priors$d12_beta_a, "header") <- d12_header
    attr(priors$d12_beta_b, "header") <- gsub(x = d12_header, pattern = "alpha",
                                              replacement = "beta")

    ## correlations among multiple cov_d?
    if (Nvars$d1 > 1L) {
      priors$d1_beta <- matrix(1.5, nrow = Nvars$d1, ncol = Nvars$d1,
                               dimnames = list(names(cov_d), names(rr.data)))
      diag(priors$d1_beta) <- 1 # ignored
      class(priors$d1_beta) <- c("lavaan.matrix", class(priors$d1_beta))
      d1_header <- paste('The matrix below contains shape hyperparameters',
                         '"alpha" (lower triangle) and "beta" (upper triangle)',
                         'for the beta-distribution priors of correlations among',
                         'dyad-level covariates. Diagonal values are ignored.')
      attr(priors$d1_beta, "header") <- d1_header
    }
  }

  ## CASE-level correlations
  cNames <- paste0(rep(names(rr.data), each = 2), c("_out", "_in"))
  if (!is.null(cov_d)    ) cNames <- c(cNames, names(cov_d))
  if (!missing(case_data)) cNames <- c(cNames, names(case_data))

  priors$case_beta <- matrix(1.5, nrow = Nvars$case, ncol = Nvars$case,
                             dimnames = list(cNames, cNames))
  diag(priors$case_beta) <- 1 # ignored
  class(priors$case_beta) <- c("lavaan.matrix", class(priors$case_beta))

  ## conditionally list variables/components in case-level matrix
  c_header <- paste('The matrix below contains shape hyperparameters',
                    '"alpha" (lower triangle) and "beta" (upper triangle)',
                    'for the beta-distribution priors of correlations among',
                    'case-level components of each round-robin variable')
  if (length(Nvars$case) > 2*Nvars$rr) {
    c_header <- paste0(c_header, ", followed by")
    if (!is.null(cov_d)) {
      c_plus <- "case-level components of dyad-level covariates"
    } else c_plus <- character(0)
    if (!missing(case_data)) {
      c_plus <- c(c_plus,
                  paste0(ifelse(modelG, yes = "case-level components of ", no = ""),
                         "case-level covariates"))
    }
    c_header <- paste(c_header, paste(c_plus, collapse = " and "))
  }
  attr(priors$case_beta, "header") <- paste0(c_header,
                                             '. Diagonal values are ignored.')

  ## GROUP-level correlations
  if (modelG) {
    gNames <- names(rr.data)
    if (!is.null(cov_d)     ) gNames <- c(gNames, names(cov_d))
    if (!missing(case_data) ) gNames <- c(gNames, names(case_data))
    if (!missing(group_data)) gNames <- c(gNames, names(group_data))

    priors$group_beta <- matrix(1.5, nrow = Nvars$group, ncol = Nvars$group,
                                dimnames = list(gNames, gNames))
    diag(priors$group_beta) <- 1 # ignored
    class(priors$group_beta) <- c("lavaan.matrix", class(priors$group_beta))

    ## conditionally list variables/components in group-level matrix
    g_header <- paste('The matrix below contains shape hyperparameters',
                      '"alpha" (lower triangle) and "beta" (upper triangle)',
                      'for the beta-distribution priors of correlations among',
                      'group-level components of each round-robin variable')
    if (length(Nvars$group) > Nvars$rr) {
      g_header <- paste0(g_header, ", followed by")
      ## check each level
      if (!is.null(cov_d)) {
        g_plus <- "group-level components of dyad-level covariates"
      } else g_plus <- character(0)
      if (!missing(case_data)) {
        g_plus <- c(g_plus, "group-level components of case-level covariates")
      }
      if (!missing(group_data)) {
        g_plus <- c(g_plus, "group-level covariates")
      }
      ## add "and"?
      if (length(g_plus) > 1L) {
        g_plus[length(g_plus)] <- paste("and", g_plus[length(g_plus)])
      }
      ## add commas?
      if (length(g_plus) > 2L) {
        g_plus <- paste(g_plus, collapse = ", ")
      }
      g_header <- paste(g_header, g_plus)
    }
    attr(priors$group_beta, "header") <- paste0(g_header,
                                                '. Diagonal values are ignored.')
  }


  ## end CORRELATIONS


  ## priors for MEANS
  if (modelM) {
    priors$rr_Mvec_m  <- sapply(rr.data, median, na.rm = TRUE)
    priors$rr_Mvec_sd <- rrSD
    #TODO: separate Mvec (RR) from group-means of each level's covariate(s)?
  }

  priors
}




