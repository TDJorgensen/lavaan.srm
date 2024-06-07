### Terrence D. Jorgensen
### Last updated: 7 June 2024
### function to implement Stage-1 of 2-stage SR-SEM estimator


##' MCMC estimation of multivariate SRM
##'
##' Estimate a multivariate social relations model (mvSRM; Nestler,
##' [2018](https://doi.org/10.3102/1076998617741106)) using the Stan
##' (Carpenter et al., [2017](https://doi.org/10.18637/jss.v076.i01)) software.
##' Stan uses MCMC estimation, enabling Bayesian inferential methods for SRM
##' parameters.  In addition to extending the mvSRM with level-specific
##' covariates (including the round-robin-group level of analysis), this result
##' provides Stage-1 estimates to use as input for a 2-stage estimator of the
##' social relations structural equation model (SR-SEM; Nestler et al.,
##' [2020](https://doi.org/10.1007/s11336-020-09728-z)).
##'
##'
##' @param data Dyad-level `data.frame` containing round-robin variables, as
##'   well as symmetric dyad-level variables (i.e., constant within dyads).
##'   Must be in long format (one column per variable, 2 rows per dyad) and
##'   contain case-level ID variables. Cases are often individual subjects,
##'   but could be households or countries.  Can also contain group-level ID
##'   (named by the `IDgroup=` argument), but cannot contain any other types of
##'   variable.  Case-level and group-level covariates are passed to other
##'   arguments.
##'   *Optionally*, `data` can be a 3-dimensional
##'   `array` (for 1 round-robin group) or a `list` of such arrays.  The first
##'   2 dimensions of the `array`(s) must be of equal `length` (representing the
##'   same people in rows and columns of a sociomatrix), and every "slice" of
##'   the third dimension is a different round-robin variable (so the `dimnames`
##'   attribute must have names for the third dimension).
##' @param rr.vars Optional `character` indicating a subset of columns in `data`
##'   to fit the mvSRM to.  By default (`NULL`), all columns are analyzed.
##' @param IDout `character` indicating the name of the case-level ID variable
##'   in `data` that represents the out-going effect (ego, actor, or perceiver)
##' @param IDin `character` indicating the name of the case-level ID variable
##'   in `data` that represents the in-coming effect (alter, partner, or target)
##' @param IDgroup Optional `character` indicating the name of the ID variable
##'    in `data` that distinguishes between multiple round-robin groups
##'
##' @param fixed.groups `logical` indicating whether means of multiple
##'   round-robin groups should be treated as fixed effects.  When `TRUE`,
##'   means will be "partialed out" by group-mean centering each round-robin
##'   variable, as well as any case-level covariates in `case_data`.
##' @param group_data Optional `data.frame` containing group-level covariates.
##'   Must contain the same `IDgroup` variable as `data`.
##'   Ignored when `fixed.groups=TRUE`
##'   Not used yet.
##' @param case_data Optional `data.frame` containing case-level covariates.
##'   Must contain a variable named `ID`, with **all the same values** that
##'   appear in `data[[IDout]]` and `data[[IDin]]` when `is.data.frame(data)`;
##'   or all the same `dimnames(data)[1:2]` when `data=` is a matrix/array (or
##'   all `dimnames()` when `data=` is a list of arrays).
##' @param block Optional `character` indicating the name of a case-level
##'   grouping variable in `data` that distinguishes between groups of people
##'   within the same round-robin group.
##'   Not used yet.
##' @param priors A named `list` of hyperparameters for prior distributions
##'   The [srm_priors()] function provides default priors, and users can edit
##'   its output to substitute more informative hyperparameters.
##' @param return_stan_data `logical`. Set `TRUE` to return the list passed to
##'   `rstan::sampling(data=)`. Helpful for creating reprex when Stan fails.
##' @param return_stanfit `logical`. Set `TRUE` to return the
##'   \code{\linkS4class{stanfit}}-class object, without making it a
##'   \code{\linkS4class{mvSRM}}-class object.
##' @param saveComp `logical` indicating whether to save the posterior samples
##'   of group- (if relevant), case-, and dyad-level components of each variable
##'   in `data=` (and decomposed `case_data=`, if relevant)
##' @param ... Arguments passed to [rstan::sampling()] (e.g. iter, chains).
##'
##' @return An object of \code{\linkS4class{mvSRM}}, which inherits from a
##'   \code{\linkS4class{stanfit}} object returned by [rstan::sampling()].
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
##'  Nestler, S., Lüdtke, O., & Robitzsch, A. (2020). Maximum likelihood
##'  estimation of a social relations structural equation model.
##'  *Psychometrika, 85*(4), 870--889.
##'  <https://doi.org/10.1007/s11336-020-09728-z>
##'
##' @seealso \code{\linkS4class{mvSRM}}
##' @examples
##'
##'
##' \dontrun{
##' ## use example data from the ?srm::srm help page
##' data(data.srm01, package="srm")
##' ## do NOT use the same case-level IDs across groups
##' data.srm01$egoID   <- paste(data.srm01$Group, data.srm01$Actor  , sep = "_")
##' data.srm01$alterID <- paste(data.srm01$Group, data.srm01$Partner, sep = "_")
##'
##' srmOut <- mvsrm(data = data.srm01, rr.vars = paste0("Wert", 1:3),
##'                 IDout = "egoID", IDin = "alterID", IDgroup = "Group",
##'                 ## rstan::sampling() arguments to run this quickly
##'                 chains = 2, iter = 20, seed = 12345,
##'                 cores = ifelse(parallel::detectCores() < 3L, 1L, 2L))
##'
##'
##' ## default: posterior means of round-robin variable SDs and correlations
##' summary(srmOut)
##'
##' ## round-robin variable means are always a group-level statistic
##' summary(srmOut, srm.param = "mean")
##'
##' ## other point (and interval) estimates from posterior
##' summary(srmOut, srm.param = "mean", posterior.est = "median",
##'         ## choose 99% highest-density uncertainty intervals, overwrite
##'         interval = "hdi", credMass = .99)
##'
##' ## for the mode, other arguments passed to modeest::mlv()
##' summary(srmOut, srm.param = "mean",
##'         posterior.est = "mode", method = "shorth")
##'
##' ## SDs & correlations of group effects
##' summary(srmOut, component = "group")
##'
##'
##' ## group-mean centered results (no group-level effects)
##' gmcOut <- mvsrm(data = data.srm01, rr.vars = paste0("Wert", 1:3),
##'                 IDout = "egoID", IDin = "alterID", IDgroup = "Group",
##'                 fixed.groups = TRUE,
##'                 chains = 2, iter = 20, seed = 12345,
##'                 cores = ifelse(parallel::detectCores() < 3L, 1L, 2L))
##' summary(gmcOut)
##'
##'
##' ## single round-robin group: no group-level random effects or (co)variances
##' g1out <- mvsrm(data = data.srm01[data.srm01$Group == 1,],
##'                rr.vars = paste0("Wert", 1:3),      # still an option:
##'                IDout = "egoID", IDin = "alterID",  # fixed.groups = TRUE,
##'                chains = 2, iter = 20, seed = 12345,
##'                cores = ifelse(parallel::detectCores() < 3L, 1L, 2L))
##' summary(g1out)
##' summary(g1out, srm.param = "mean") # available unless fixed.groups = TRUE
##' }
##'
##' @importFrom rstan sampling
##' @importFrom methods as
##' @export
mvsrm <- function(data, rr.vars = NULL, IDout, IDin, IDgroup = NULL,
                  fixed.groups = FALSE, group_data = NULL,  case_data = NULL,
                  block = NULL, priors, return_stan_data = FALSE,
                  return_stanfit = FALSE, saveComp = FALSE, ...) {
  MC <- match.call(expand.dots = TRUE) # to store in mvSRM-class slot

  if (is.data.frame(data)) {
    ## assume it is univariate long format
    if (!is.null(rr.vars)) {
      ## extract selected columns + IDs
      longUni <- data[ , c(IDgroup, IDout, IDin, rr.vars), drop = FALSE]
    } else {
      longUni <- data
      ## save round-robin variable names
      rr.vars <- setdiff(names(data), c(IDgroup, IDout, IDin))
    }
    ## convert to bivariate long format
    Yd2 <- longUni2longBi(y = rr.vars, data = longUni, groupMC = fixed.groups,
                          IDout = IDout, IDin = IDin, group = IDgroup)

  } else if (is.array(data) || is.list(data)) {
    ## extract selected round-robin variables
    if (!is.null(rr.vars)) {
      if (is.array(data)) {
        SQ <- data[ , , rr.vars]
      } else {
        SQ <- lapply(data, function(x) x[ , , rr.vars])
      }
    } else SQ <- data

    ## save round-robin variable names
    if (is.null(rr.vars)) {
      rr.vars <- if (is.list(SQ)) {
        dimnames(SQ[[1]])[[3]] # use first RR group's dimnames
      } else dimnames(SQ)[[3]] # use only  RR group's dimnames
    }

    ## convert to long (bivariate) format
    Yd2 <- square2longBi(mat = SQ, groupMC = fixed.groups,
                         dropNA = is.array(data), #FIXME: handle NAs in Stan
                         label = ifelse(is.null(rr.vars), "y", rr.vars),
                         group = ifelse(is.null(IDgroup), "group", IDgroup))
    if (is.null(IDgroup) && !is.null(Yd2$group)) IDgroup <- "group"

  } else stop('Unrecognized format for data= argument')

  #FIXME: Remove this after missing data are allowed?
  #       Or add missing = c("augment","listwise") option?
  Yd2 <- Yd2[complete.cases(Yd2),]


  ## check for enough groups
  if (!is.null(IDgroup) && !fixed.groups) {
    Ng <- length(unique(Yd2[ , IDgroup]))
    msg <- 'fixed.groups=TRUE is recommended.'
    if (Ng == 1L) stop('Only 1 value of IDgroup= was found.  Group-level ',
                       'variance cannot be calculated.\nLeave IDgroup=NULL.')
    if (Ng == 2L) stop('Only 2 IDgroup= values found. Group-level correlations',
                       ' can only be \u22121, 0, or +1.\n', msg)
    if (Ng >= 3L && Ng < 10) warning('Group-level sample size might be too ',
                                     'small to model group-level (co)variance.',
                                     '\n', msg)
  } else Ng <- 0L

  ## check that variables are numeric
  #TODO: allow ordered() factors?
  rrNum <- sapply(Yd2[-1:ifelse(is.null(IDgroup), -2, -3)], is.numeric)
  if (!all(rrNum)) stop('Only numeric round-robin data are currently allowed')
  ## repeat below for group/case-level data

  ## store dyad-constant covariates separately
  Yd1 <- Yd2[c(IDgroup, "ID_i", "ID_j")]
  dyad_constant_vars <- character(0)
  ## store round-robin variable names for mvSRM-class slot
  rr.names <- list()

  for (rr in rr.vars) {
    sym <- Yd2[, paste0(rr, "_ij")] == Yd2[, paste0(rr, "_ji")]

    #FIXME: enable data-augmentation in Stan script; set flag here?
    dyadNAs <- which(is.na(sym))
    if (length(dyadNAs)) {
      NAlist <- paste('Cases', Yd2$ID_i[dyadNAs], 'and', Yd2$ID_j[dyadNAs])
      if (is.null(IDgroup)) {
        NAlist <- paste(NAlist, "\n\t")
      } else NAlist <- paste(NAlist, 'in group', Yd2[dyadNAs, IDgroup], "\n\t")

      stop('Currently, each observed dyad must have complete data on all ',
           'round-robin variables.  Missing data on variable ', rr,
           ' found in the following dyad(s):\n\t', NAlist)
    }

    if (all(sym)) {
      ## copy 1 column to dyad-constant data
      Yd1[, rr] <- Yd2[, paste0(rr, "_ij")]
      ## remove from round-robin data
      Yd2[, paste0(rr, "_ij")] <- NULL
      Yd2[, paste0(rr, "_ji")] <- NULL
      ## remove from list of round-robin variables
      rr.vars <- setdiff(rr.vars, rr)
      ## save to list of dyad-constant variables
      dyad_constant_vars <- c(dyad_constant_vars, rr)
    } else {
      rr.names[[rr]] <- paste0(rr, c("_ij", "_ji"))
    }
  }

  ## synchronize names of ID variables across dyad/case/group-level data sets

  if (!is.null(group_data)) {
    if (fixed.groups) stop('Cannot model group_data when group-level variance ',
                           'is removed via fixed.groups=TRUE')
    ## need group ID if there are group covariates
    stopifnot(!is.null(IDgroup))
    if (!IDgroup %in% colnames(group_data)) {
      stop('The "IDgroup" variable is missing from "group_data"')
    }
    if (ncol(group_data) < 2L) {
      stop('If "group_data" is provided, it must contain the "IDgroup" ',
           'variable and at least one covariate')
    }
    ## ensure group IDs appear first
    group_data <- group_data[, c(IDgroup, setdiff(colnames(group_data), IDgroup))]
    ## check that variables are numeric
    #TODO: allow ordered() factors?
    groupNum <- sapply(group_data[-1], is.numeric)
    if (!all(groupNum)) stop('Only numeric group-level data are currently allowed')
  }

  if (!is.null(case_data)) {
    if (!("ID" %in% colnames(case_data))) {
      stop('The "case_data=" argument must contain a casewise ID variable ',
           'named "ID".')
    }

    if (!is.null(IDgroup)) {
      ## verify there is a group ID variable
      #FIXME: unnecessary? can be assigned within Stan program
      if (!IDgroup %in% colnames(case_data))
        stop('The group-ID variable is missing from "case_data=". It must ',
             'contain group IDs in a variable with the same name as in "data="',
             ' (i.e., the "IDgroup=" argument).')
    }

    if (ncol(case_data) < (2L + !is.null(IDgroup)))
      stop('If "case_data" is provided, it must contain a casewise ID variable',
           if (!is.null(IDgroup)) ', a group-ID variable,',
           ' and at least one case-level covariate')

    ## ensure group and case IDs appear first
    case_data <- case_data[ , c(IDgroup, "ID",
                                setdiff(colnames(case_data), c(IDgroup, "ID")))]
    ## check that variables are numeric
    #TODO: allow ordered() factors?
    caseNum <- sapply(case_data[-1:ifelse(is.null(IDgroup), -1, -2)], is.numeric)
    if (!all(caseNum)) stop('Only numeric case-level data are currently allowed')
  }


  ## save ID values from dyad-level data to check levels in covariate data
  dyad_gIDs <- if (is.null(IDgroup)) NULL else Yd2[ , IDgroup]
  dyad_pIDs <- Yd2[ , c("ID_i", "ID_j")]
  dyad_pIDs_vec <- c(dyad_pIDs, recursive = TRUE, use.names = FALSE)

  ## check for duplicate case-IDs across groups
  if (!is.null(IDgroup)) {
    pgID_tab <- table(dyad_pIDs_vec, c(dyad_gIDs, dyad_gIDs))
    pgID_dup <- apply(pgID_tab, MARGIN = 1, function(x) sum(x > 0L) > 1L)
    if (any(pgID_dup)) {
      dupList <- apply(pgID_tab[pgID_dup, , drop = FALSE], 1,
                       function(x) which(x != 0), simplify = FALSE)
      dupMessage <- character(length(dupList))
      for (i in names(dupList)) {
        dupMessage[i] <- paste('\tcase ID', i, 'appears in the following groups:',
                               paste(dupList[[i]], collapse = ', '))
      }
      stop('The following case-IDs in data[IDout] or data[IDin] were found in ',
           'multiple groups:', paste(dupMessage, collapse = '\n'))
    }
  }


  ## verify group/case-level IDs match across level-specific data sets

  ## GROUP-level covariates
  if (!is.null(group_data)) {
    ## any group IDs in data missing from group_data?
    not_in_gcov <- setdiff(dyad_gIDs, group_data[ , IDgroup])
    if (length(not_in_gcov)) {
      stop('Found the following "IDgroup" level(s) in "data" but not in ',
           '"group_data":\n', paste(not_in_gcov, collapse = ", "))
    }
    ## any group IDs from group_data missing in data?
    g_not_in_data <- setdiff(group_data[ , IDgroup], dyad_gIDs)
    if (length(g_not_in_data)) {
      stop('Found the following "IDgroup" level(s) in "group_data" but not in ',
           '"data":\n', paste(g_not_in_data, collapse = ", "))
    }

    ## use group_data$IDgroup values to set levels in case/dyad data
    gID_levels <- group_data[ , IDgroup]
  } else gID_levels <- unique(dyad_gIDs) # can overwrite below

  ## CASE-level covariates
  if (!is.null(case_data)) {

    ## any person IDs in data missing from case_data?
    p_not_in_pcov <- setdiff(dyad_pIDs_vec, case_data$ID)
    if (length(p_not_in_pcov)) {
      stop('Found the following "IDout/IDin" level(s) in "data" but not in ',
           '"case_data":\n', paste(p_not_in_pcov, collapse = ", "))
    }
    ## any person IDs in case_data missing from data?
    p_not_in_data <- setdiff(case_data$ID, dyad_pIDs_vec)
    if (length(p_not_in_data)) {
      stop('Found the following "IDout/IDin" level(s) in "case_data" but not',
           ' in "data":\n', paste(p_not_in_data, collapse = ", "))
    }


    if (!is.null(IDgroup)) {
      ## any group IDs in data missing from case_data?
      pg_not_in_pcov <- setdiff(dyad_gIDs, case_data[ , IDgroup])
      if (length(pg_not_in_pcov)) {
        stop('Found the following "IDgroup" level(s) in "data" but not in ',
             '"case_data":\n', paste(pg_not_in_pcov, collapse = ", "))
      }
      ## any group IDs in case_data missing from data?
      pg_not_in_data <- setdiff(case_data[ , IDgroup], dyad_gIDs)
      if (length(pg_not_in_data)) {
        stop('Found the following "IDgroup" level(s) in "case_data" but not in',
             ' "data":\n', paste(pg_not_in_data, collapse = ", "))
      }

      if (is.null(group_data)) {
        ## if no group_data, use case_data$IDgroup values to set dyad-data levels
        gID_levels <- unique(case_data[ , IDgroup])
      }
    }

    ## use case_data$ID values to set levels in case/dyad data
    pID_levels <- case_data$ID
  } else pID_levels <- unique(dyad_pIDs_vec) # can overwrite below


  ## ensure IDs contain integers from 1:N(p/g)

  ## GROUP-level IDs
  ## ~~~~~~~~~~~~~~~
  if (!is.null(IDgroup)) {
    ## DYAD-level data
    Yd2[, IDgroup] <- unclass(factor(Yd2[, IDgroup], levels = gID_levels))
    ## GROUP-level data
    if (!is.null(group_data)) {
      group_data[, IDgroup] <- unclass(factor(group_data[, IDgroup],
                                              levels = gID_levels))
    }
    ## CASE-level data
    if (!is.null(case_data)) {
      case_data[, IDgroup] <- unclass(factor(case_data[, IDgroup],
                                             levels = gID_levels))
    }
  }

  ## CASE-level IDs
  ## ~~~~~~~~~~~~~~~
  ## DYAD-level data
  Yd2$ID_i <- unclass(factor(Yd2$ID_i, levels = pID_levels))
  Yd2$ID_j <- unclass(factor(Yd2$ID_j, levels = pID_levels))
  ## CASE-level data
  if (!is.null(case_data)) {
    case_data$ID <- unclass(factor(case_data$ID, levels = pID_levels))

    ## (group-)mean center variables?
    if (fixed.groups) {
      #FIXME: don't allow with incomplete data
      if (is.null(IDgroup)) {
        ## grand mean center
        case_z <- data.frame(scale(case_data[,-1, drop = FALSE], scale = FALSE))
        case_data <- cbind(case_data[,1, drop = FALSE], data.frame(case_z))

      } else {
        ## group mean center
        case_gM <- aggregate(x = case_data[, -1:-2, drop = FALSE],
                             by = case_data[IDgroup], FUN = mean)
        for (g in case_gM[,IDgroup]) {
          case.idx  <- case_data[,IDgroup] == g
          group.idx <- case_gM[ , IDgroup] == g
          for (v in names(case_data)[-1:-2]) {
            case_data[case.idx, v] <- case_data[case.idx, v] - case_gM[group.idx, v]
          }
        }

      }
    }

  }

  ## store IDs in list for mvSRM-class slot
  ID_slot <- list(dyad = Yd2[ , c(IDgroup, "ID_i", "ID_j")],
                  case = case_data[, c(IDgroup, "ID")])
  if (!is.null(IDgroup) && !fixed.groups) {
    ID_slot <- c(ID_slot, list(group = group_data[, IDgroup, drop = FALSE]))
  }


  ## format data for rstan::sampling(),
  ## while assembling the correct model name
  SRMname <- c(design = paste0("RR",ifelse(length(dyad_constant_vars), "d", ""),
                               ifelse(length(dyad_constant_vars) == 1, "1", ""),
                               ifelse(!is.null(case_data), "c", ""))) # add groups below
  SRMname["scale"] <- "con" #TODO: check "ord" or "mix" when thresholds available
  SRMname["NAs"] <- "comp" #TODO: enable "miss" with covariate suffixes "D", "C", or "G"

  ## sample sizes: Nd, Np, Ng
  knowns <- list(Nd = nrow(Yd2),
                 Np = length(pID_levels),
                 ## number of observed measures
                 Kd1 = length(dyad_constant_vars),
                 Kd2 = length(rr.vars))

  ## observed-data: Yd2, Yd1, Yp, Yg
  knowns$Yd2 <- as.matrix(Yd2[, -1:ifelse(is.null(IDgroup), -2, -3)])

  ## begin assembling a call for default hyperparameters (in case not provided)
  priorCall <- list(quote(srm_priors))
  ## in order to call srm_priors() below, make sure "data" is long-uni format
  dat_ij <- Yd2[paste0(rr.vars, "_ij")]
  dat_ji <- Yd2[paste0(rr.vars, "_ji")]
  ## set common names for round-robin variables
  names(dat_ij) <- rr.vars
  names(dat_ji) <- rr.vars
  priorCall$data <- rbind(dat_ij, dat_ji)

  if (length(dyad_constant_vars)) {
    #TODO: knowns$Yd1 <- as.matrix(Yd1[, -1:ifelse(is.null(IDgroup), -2, -3)])
    #TODO: priorCall$cov_d <- Yd1[, -1:ifelse(is.null(IDgroup), -2, -3), drop = FALSE]
    #FIXME: or is cov_d inferred from argument rr.vars= argument?
  }

  ## case IDs and data
  knowns$IDp <- as.matrix(Yd2[, c("ID_i", "ID_j")])
  if (!is.null(case_data)) {
    knowns$Kp <-      ncol(case_data)   -  ifelse(is.null(IDgroup), 1L, 2L) # minus 1 or 2 ID variables
    knowns$Yp <- as.matrix(case_data[ , -1:ifelse(is.null(IDgroup), -1, -2)])
    knowns$IDpp <- case_data$ID
    priorCall$case_data <- case_data[ , -1:ifelse(is.null(IDgroup), -1, -2), drop = FALSE]
    #TODO: use aggregate() or anova() to decompose case vs. group variance,
    #      then assign to priorCall$decomp_c
  } else knowns$Kp <- 0L

  ## group IDs and data
  if (!is.null(IDgroup) && !fixed.groups) {
    knowns$Ng <- length(gID_levels)
    knowns$IDg <- as.array(Yd2[, IDgroup]) # 1 dimension (not a matrix)
    if (!is.null(group_data)) {
      knowns$Kg <-      ncol(group_data) - 1L
      #TODO: knowns$Yg <- as.matrix(group_data[, -1])
      #TODO: priorCall$group_data <- group_data[, -1, drop = FALSE]
      SRMname["design"] <- paste0(SRMname["design"], "g")
    } else knowns$Kg <- 0L

    ## group ID per case
    if (!is.null(case_data)) {
      #TODO: knowns$IDgp <- case_data[[IDgroup]]
    }

    SRMname["groups"] <- "gN"

    ## model both group-level effects and means
    priorCall$modelG <- TRUE
    priorCall$modelM <- TRUE

  } else if (fixed.groups) {
    knowns$Kg <- 0L # group-mean centered, so no group-level covariates
    SRMname["groups"] <- "g0"
    ## model neither group-level effects nor means
    priorCall$modelG <- FALSE
    priorCall$modelM <- FALSE

  } else {
    knowns$Kg <- 0L # single group, so no group-level covariates
    SRMname["groups"] <- "g1" # must be is.null(IDgroup)
    ## model means, but not group-level effects
    priorCall$modelG <- FALSE
    priorCall$modelM <- TRUE
  }

  default_priors <- eval(as.call(priorCall))
  ## any user-specified hyperparameters?
  if (!missing(priors)) {
    stopifnot(is.list(priors))

    ## check that each replacement matches the expected format
    names2replace <- intersect(names(priors), names(default_priors))
    for (n in names2replace) {

      ## verify classes match
      if (!isTRUE(all.equal(class(priors[[n]]), class(default_priors[[n]])))) {
        stop('User-specified priors= element ', dQuote(n), ' has class ',
             paste(dQuote(class(priors[[n]])), collapse = ","), 'rather than ',
             'expected class: ', dQuote(class(priors[[n]])))
      }

      ## verify dimensions match
      if (inherits(priors[[n]], c("data.frame","matrix"))) {
        user_dim <- dim(        priors[[n]])
        exp_dim  <- dim(default_priors[[n]])
        if (!isTRUE(all.equal(user_dim, exp_dim))) {
          stop('User-specified priors= element ', dQuote(n), ' has ',
               paste(user_dim, c("rows","columns"), collapse = " and "),
               'rather than expected dimensions: [',
               paste(exp_dim, collapse = ","), ']')
        }

      } else if (inherits(priors[[n]], c("integer","numeric"))) {
        user_len <- length(        priors[[n]])
        exp_len  <- length(default_priors[[n]])
        if (user_len != exp_len) {
          stop('User-specified priors= vector ', dQuote(n), ' has ', user_len,
               ' values, but requires ', exp_len, 'values')
        }

      } else stop('Invalid class for priors= element ', dQuote(n))
    }

    ## any required hyperparameters not specified by user?
    names2add <- setdiff(names(default_priors), names(priors))
    if (length(names2add)) {
      ## Use default values
      for (n in names2add)  priors[[n]] <- default_priors[[n]]
    }

    # else, priors is missing, so use defaults:
  } else priors <- default_priors

  knowns <- c(knowns, priors)

  ## add number of components (allKd, allKp, and maybe allKg)
  Nvars <- attr(default_priors, "Nvars")
  knowns$allKd <- knowns$Kd2*2 + knowns$Kd1
  knowns$allKp <- Nvars$case
  knowns$allKg <- Nvars$group

  #TODO: after priors are chosen, then replace NAs in knowns$Y* with missCode=


  if (return_stan_data) return(knowns)

  ## save names of unknown parameters for Stan to sample
  mu <- if (fixed.groups) NULL else "Mvec"
  sigma <- c("s_rr","S_p")
  corr <- c("r_d2","Rd2","Rp")
  derived <- c("Rsq", "dSigma","pSigma")
  if (!is.null(IDgroup) && !fixed.groups) {
    sigma <- c(sigma, "S_g")
    corr  <- c(corr, "Rg")
    derived <- c(derived, "gSigma")
  }
  if (saveComp) {
    ## add the group-, case-, and dyad-level components
    ranFX <- c("Yd2e","AP")
    if (!is.null(IDgroup) && !fixed.groups) ranFX <- c(ranFX, "GG")
  } else ranFX <- NULL
  #TODO: option to save imputed observations (with random effects)
  ## imps <- "augYd","augYp",  if (!fixed.groups) "augYg"
  unknowns <- c(mu, sigma, corr, derived, ranFX)

  ## check for ...$init_r, set default to 0.5
  dots <- list(...)
  if (is.null(dots$init_r)) dots$init_r <- 0.5
  ## assemble args in a call
  stanCall <- c(list(quote(rstan::sampling),
                     object = stanmodels[[paste(SRMname, collapse = "_")]],
                     data = knowns, pars = unknowns), dots)
  fit <- try(eval(as.call(stanCall)), silent = TRUE)
  ## check whether it worked
  if (inherits(fit, "try-error")) {
    message('rstan::sampling() call was unsuccessful. Returned the error ',
            'message as a "try-error" object (see ?try documentation)')
    return(fit)
  }

  if (return_stanfit) return(fit)

  ## add slots to construct mvSRM object
  fit      <- as(fit, "mvSRM")
  fit@call <- MC
  fit@nobs <- c(group = knowns$Ng, case = knowns$Np, dyad = knowns$Nd)
  fit@IDs  <- ID_slot

  fit@varNames <- list(RR    = rr.names,
                       dyad  = dyad_constant_vars,
                       case  = c(paste0(rep(names(rr.names), each = 2),
                                        c("_out", "_in")),
                                 #TODO: add dyad-level covariates (case components)
                                 colnames(case_data[, -1:ifelse(is.null(IDgroup),
                                                                -1, -2),
                                                    drop = FALSE])),
                       #TODO: add group components of dyad/case-level covariates
                       group = c(colnames(group_data[, -1, drop = FALSE])))
  if ("S_g" %in% sigma) fit@varNames$group <- c(names(rr.names),
                                                #TODO: add dyad/case covariates
                                                fit@varNames$group)

  fit@parNames <- list(mu      = mu,
                       sigma   = sigma,
                       corr    = corr,
                       derived = derived)
  if (saveComp) fit@parNames$components <- ranFX
  #TODO: if (saveImps) fit@parNames$imputations <- imps

  fit
}

#TODO:
#   - model A(=P) random effects for dyad-constant covariates
# Stan models:
#   - univariate SRM (with level-specific covariates)
#   - integrate covariates at each level
#   - allow incomplete data at each level
#   - allow ordered/mixed data
#   - enable block designs (fSRM as a special case: blocks of n=1)

