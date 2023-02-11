### Terrence D. Jorgensen
### Last updated: 11 February 2023
### function to implement Stage-1 of 2-stage SR-SEM estimator


##' MCMC estimation of multivariate SRM
##'
##' Estimate a multivariate social relations model (mvSRM; Nestler,
##' [2018](https://doi.org/10.3102/1076998617741106)) using the Stan
##' (Carpenter et al., [2017](https://doi.org/10.18637/jss.v076.i01)) software.
##' Stan uses MCMC estimation, enabling Bayesian inferential methods for SRM
##' parameters.  In addition to extending the mvSRM with level-specific
##' covariates (including the round-robin group-level of analysis), this result
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
##'   Not used yet.
##' @param block Optional `character` indicating the name of a case-level
##'   grouping variable in `data` that distinguishes between multiple
##'   round-robin groups.
##'   Not used yet.
##' @param return_stan_data `logical`. Set `TRUE` to return the list passed to
##'   `rstan::sampling(data=)`. Helpful for creating reprex when Stan fails.
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
##'   Betancourt, M., ... & Riddell, A. (2017).
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
##' ## posterior means (default) of round-robin variable means
##' as.matrix(srmOut, stat = "mean")
##' ## other point estimates from posterior
##' as.matrix(srmOut, stat = "mean", point = "median")
##' ## other arguments passed to modeest::mlv()
##' as.matrix(srmOut, stat = "mean", point = "mode", method = "shorth")
##'
##' ## level-specific summary statistics
##'
##' as.matrix(srmOut, stat = "sd", component = "group") # SDs of group effects,
##' as.matrix(srmOut, stat = "sd", component = "case")  # case (person) effects,
##' as.matrix(srmOut, stat = "sd", component = "dyad")  # relationship effects
##'
##' as.matrix(srmOut, stat = "cor", component = "group") # group effects,
##' as.matrix(srmOut, stat = "cor", component = "case")  # case (person) effects
##' as.matrix(srmOut, stat = "cor", component = "dyad")  # relationship effects
##'
##'
#TODO: replace with summary() method for mvSRM class
##'
##'  ## group-mean centered results (no group-level effects)
##'  gmcOut <- mvsrm(data = data.srm01, rr.vars = paste0("Wert", 1:3),
##'                  IDout = "egoID", IDin = "alterID", IDgroup = "Group",
##'                  fixed.groups = TRUE,
##'                  chains = 2, iter = 20, seed = 12345,
##'                  cores = ifelse(parallel::detectCores() < 3L, 1L, 2L))
##'  }
##'
##' @importFrom rstan sampling
##' @importFrom methods as
##' @export
mvsrm <- function(data, rr.vars = NULL, IDout, IDin, #TODO: na.code = -9999L,
                  IDgroup = NULL, fixed.groups = FALSE, group_data = NULL,
                  case_data = NULL, block = NULL, return_stan_data = FALSE,
                  saveComp = FALSE, ...) {
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

  ## check for enough groups
  if (!is.null(IDgroup) && !fixed.groups) {
    Ng <- length(unique(Yd2[ , IDgroup]))
    msg <- 'fixed.groups=TRUE is recommended.'
    if (Ng == 1L) stop('Only 1 value of IDgroup= was found.  Group-level ',
                       'variance cannot be calculated.\nLeave IDgroup=NULL.')
    if (Ng == 2L) stop('Only 2 IDgroup= values found. Group-level correlations',
                       ' can only be \u22121, 0, or +1.\n', msg)
    if (Ng >= 3L && Ng < 10) warning('Insufficient group-level sample size to model ',
                                     'group-level (co)variance.\n', msg)
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
      stop('Currently, each observed dyad must have complete data on all ',
           'round-robin variables.  Missing data on variable', rr,
           'found in the following dyad(s):\n\t',
           paste('Cases', Yd2$ID_i[dyadNAs], 'and', Yd2$IDj[dyadNAs],
                 ifelse(is.null(IDgroup), NULL,
                        paste('in group', Yd2[dyadNAs, IDgroup])),
                 sep = "\n\t"))
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
    not_in_gcov <- setdiff(group_data[ , IDgroup], dyad_gIDs)
    if (length(not_in_gcov)) {
      stop('Found the following "IDgroup" level(s) in "data" but not in ',
           '"group_data":\n', paste(not_in_gcov, collapse = ", "))
    }
    ## any group IDs from group_data missing in data?
    g_not_in_data <- setdiff(dyad_gIDs, group_data[ , IDgroup])
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
    p_not_in_pcov <- setdiff(case_data$ID, dyad_pIDs_vec)
    if (length(p_not_in_pcov)) {
      stop('Found the following "IDout/IDin" level(s) in "data" but not in ',
           '"case_data":\n', paste(p_not_in_pcov, collapse = ", "))
    }
    ## any person IDs in case_data missing from data?
    p_not_in_data <- setdiff(dyad_pIDs_vec, case_data$ID)
    if (length(p_not_in_data)) {
      stop('Found the following "IDout/IDin" level(s) in "case_data" but not',
           ' in "data":\n', paste(p_not_in_data, collapse = ", "))
    }


    if (!is.null(IDgroup)) {
      ## any group IDs in data missing from case_data?
      pg_not_in_pcov <- setdiff(case_data[ , IDgroup], dyad_gIDs)
      if (length(pg_not_in_pcov)) {
        stop('Found the following "IDgroup" level(s) in "data" but not in ',
             '"case_data":\n', paste(pg_not_in_pcov, collapse = ", "))
      }
      ## any group IDs in case_data missing from data?
      pg_not_in_data <- setdiff(dyad_gIDs, case_data[ , IDgroup])
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
      #TODO: save row.idx, sort by group, groupMC, re-sort
    }
  }

  ## store IDs in list for mvSRM-class slot
  ID_slot <- list(dyad = Yd2[ , c(IDgroup, "ID_i", "ID_j")],
                  case = case_data[, c(IDgroup, "ID")])
  if (!is.null(IDgroup) && !fixed.groups) {
    ID_slot <- c(ID_slot, list(group = group_data[, IDgroup, drop = FALSE]))
  }


  ## format data for rstan::sampling()

  ## sample sizes: Nd, Np, Ng
  knowns <- list(Nd = nrow(Yd2),
                 Np = length(pID_levels),
                 ## number of observed measures
                 #TODO: Kd1 = length(dyad_constant_vars),
                 Kd2 = length(rr.vars))

  ## observed-data: Yd2, Yd1, Yp, Yg
  knowns$Yd2 <- as.matrix(Yd2[, -1:ifelse(is.null(IDgroup), -2, -3)])
  #TODO: knowns$Yd1 <- as.matrix(Yd1[, -1:ifelse(is.null(IDgroup), -2, -3)])

  ## case IDs and data
  knowns$IDp <- as.matrix(Yd2[, c("ID_i", "ID_j")])
  if (!is.null(case_data)) {
    #TODO: knowns$Kp <-      ncol(case_data)   -  ifelse(is.null(IDgroup), 1L, 2L)
    #TODO: knowns$Yp <- as.matrix( case_data[, -1:ifelse(is.null(IDgroup), -1, -2)])
  } # else knowns$Yp <- matrix[1:knowns$Np]

  ## group IDs and data
  if (!is.null(IDgroup) && !fixed.groups) {
    knowns$Ng <- length(gID_levels)
    knowns$IDg <- as.array(Yd2[, IDgroup]) # 1 dimension (not a matrix)
    if (!is.null(group_data)) {
      #TODO: knowns$Kg <-      ncol(group_data) - 1L
      #TODO: knowns$Yg <- as.matrix(group_data[, -1])
    }
  }
  if (return_stan_data) return(knowns)

  ## save names of unknown parameters for Stan to sample
  mu <- if (fixed.groups) NULL else "Mvec"
  sigma <- c("s_rr","S_p")
  corr <- c("Rd2","Rp")
  derived <- c("Rsq")
  if (!is.null(IDgroup) && !fixed.groups) {
    sigma <- c(sigma, "S_g")
    corr  <- c(corr , "Rg")
  }
  if (saveComp) {
    ## add the group-, case-, and dyad-level components
    ranFX <- c("Yd2e","AP")
    if (!is.null(IDgroup) && !fixed.groups) ranFX <- c(ranFX, "GG")
  } else ranFX <- NULL
  #TODO: option to save imputed observations
  ## imps <- "impYd","impYp",  if (!fixed.groups) "impYg"
  unknowns <- c(mu, sigma, corr, derived, ranFX)

  if (fixed.groups) {
    fit <- try(sampling(stanmodels$RR_con_comp_g0, data = knowns,
                        pars = unknowns, ...), silent = TRUE)
  } else if (is.null(IDgroup)) {
    fit <- try(sampling(stanmodels$RR_con_comp_g1, data = knowns,
                        pars = unknowns, ...), silent = TRUE)
  } else {
    fit <- try(sampling(stanmodels$RR_con_comp_gN, data = knowns,
                        pars = unknowns, ...), silent = TRUE)
  }

  if (inherits(fit, "try-error")) {
    message('rstan::sampling() call was unsuccessful. Returned the error ',
            'message as a "try-error" object (see ?try documentation)')
    return(fit)
  }

  ## add slots to construct mvSRM object
  fit      <- as(fit, "mvSRM")
  fit@call <- MC
  fit@IDs  <- ID_slot

  fit@varNames <- list(RR    = rr.names,
                       dyad  = dyad_constant_vars,
                       case  = colnames(case_data[, -1:ifelse(is.null(IDgroup),
                                                              -1, -2)]),
                       group = colnames(group_data[, -1]))

  fit@parNames <- list(mu      = mu,
                       sigma   = sigma,
                       corr    = corr,
                       derived = derived)
  if (saveComp) fit@parNames$components <- ranFX

  fit
}

#TODO:
#   - group-mean center case_data when fixed.groups=TRUE
# Stan models:
#   - integrate covariates at each level
#   - allow incomplete data at each level
#   - allow ordered/mixed data



## -----------------------------
## data-transformation functions
## -----------------------------


#TODO: allow passing covariates at different levels

## transform from long-long (univariate outcome) format to square format,
## like an adjacency matrix / sociometric data
## - y = vector of round-robin variable names/indices in "data"
## - IDout, IDin = actor & partner (perceiver & target) ID names/indices
## - group = optional group ID name/index
## - returnList indicates whether multiple groups should be treated as 1 group (FALSE)
longUni2square <- function(y, data, IDout, IDin, group = NULL,
                           returnList = TRUE) {
  ## all variable names in the data?
  if (is.numeric(y)) {
    y <- colnames(data)[y]
  } else if (is.character(y)) {
    y.miss <- setdiff(y, colnames(data))
    if (length(y.miss)) stop('The following variable names in "y" are missing ',
                             'from the data set:\n', paste(y.miss, collapse = ","))
  } else stop('"y" must be a character vector indicating variable name(s) or',
              ' a numeric vector of indices in "data"')

  if (is.numeric(IDout)) IDout <- colnames(data)[IDout]
  if (is.numeric(IDin)) IDin <- colnames(data)[IDin]
  if (is.numeric(group)) group <- colnames(data)[group]
  id.miss <- setdiff(c(IDout, IDin, group), colnames(data))
  if (length(id.miss)) stop('The following ID variable names are missing from ',
                            'the data set:\n', paste(id.miss, collapse = ","))
  data <- as.data.frame(data[ , c(IDout, IDin, group, y)], stringsAsFactors = FALSE)

  ## unique IDs across groups
  #FIXME: necessary?  mvsrm() checks for duplicate case IDs.
  #       This overwrites user's IDs before mvsrm() saves them in @ID slot
  if (!is.null(group)) {
    gFac <- factor(data[ , group])

    if (!returnList) {
      ## unique IDs across groups
      data[ , group] <- as.integer(gFac)
      data[ , IDout] <- paste0("g", data[ , group], "_", data[ , IDout])
      data[ , IDin ] <- paste0("g", data[ , group], "_", data[ , IDin ])
    }
  }

  ## list of data frames (one per group)
  dataList <- if (is.null(group) || !returnList) list(data) else {
    lapply(levels(gFac), function(g) {
      subset(data, subset = eval(parse(text = paste(group, "==", g))))
    })
  }

  ## loop over groups to store y in square matrix
  outList <- lapply(dataList, function(gdat) {
    ## all person/case-level IDs
    gdat[ , IDout] <- as.character(gdat[ , IDout])
    gdat[ , IDin] <- as.character(gdat[ , IDin])
    allIDs <- unique(c(gdat[ , IDout], gdat[ , IDin]))
    ## empty array of matrices
    gmat <- array(NA, dim = c(length(allIDs), length(allIDs), length(y)),
                  dimnames = list(ego = allIDs, alter = allIDs, variable = y))
    attr(gmat, "warn") <- 0L
    ## loop over rows of data (within each y) to copy values to matrix
    warn <- FALSE
    for (yy in y) for (n in 1:nrow(gdat)) {
      ## check if there was already a value
      if (yy == y[1] && !is.na(gmat[gdat[n, IDout], gdat[n, IDin], yy]))
        attr(gmat, "warn") <- attr(gmat, "warn") + 1L
      gmat[gdat[n, IDout], gdat[n, IDin], yy] <- gdat[n, yy]
    }
    gmat
  })
  #TODO: add diagnostics, or suggest using checkDupUni()

  ## check list for multiple values in a cell
  multiVals <- sapply(outList, attr, which = "warn")
  if (!all(multiVals == 0L)) {
    for (G in which(multiVals > 0L)) {
      attr(outList[[G]], "warn") <- paste('  Group', G, 'had', multiVals[G],
                                          'cell(s) with > 1 elligible value.')
    }
    warning('\n\n"data" contained multiple values for some cells in a matrix:\n',
            paste(sapply(outList[multiVals > 0L], attr, which = "warn"),
                  collapse = "\n"),
            '\nIn each instance, the last value was kept. Please check that ',
            'each subject only provides one datum about each other subject in ',
            'data[ , y].  Store repeated measures as additional variables.\n\n')
  }
  for (G in which(multiVals == 0L)) attr(outList[[G]], "warn") <- NULL

  ## convert arrays to 2-D matrices if univariate
  if (length(y) == 1L) for (G in seq_along(outList)) outList[[G]] <- outList[[G]][ , , 1]

  ## return result (remove list if 1 group)
  if (is.null(group) || !returnList) {
    return(outList[[1]])
  } else {
    names(outList) <- levels(gFac)
  }
  outList
}


## transform from long format (bivariate outcome) to square format,
## like an adjacency matrix / sociometric data
longBi2square <- function(Ys, data, ID1, ID2, IDdyad = NULL, group = NULL) {

}


## transform from square (adjacency matrix) to long (bivariate outcome) format
## - mat = a 2-D matrix or 3-D array (third dim = variable)
## - dropNA indicates BOTH dyadic observations must be NA to drop that row
##          MUST BE IGNORED when "mat" is a 3-D array (patterns can differ)
## - label = the round-robin variable name, used as a prefix for both columns
##           when mat is a 2-D matrix (otherwise, use names of third dimension)
## - suffix = labels for first and second cases in dyad, used for ID names
##            and to distinguish between columns of round-robin variable
## - group = label used for the grouping variable when is.list(mat)
## - groupMC = whether to group mean center each RR variable
##' @importFrom utils combn
square2longBi <- function(mat, dropNA = FALSE, label = "y", group = "group",
                          suffix = c("i","j"), groupMC = FALSE) {
  if (is.list(mat)) {
    isMat <- sapply(mat, function(m) is.array(m) && is.numeric(m))
    eqDim <- sapply(mat, function(m) diff(dim(m)[1:2]) == 0)
    eqLab <- sapply(mat, function(m) all.equal(dimnames(m)[[1]], dimnames(m)[[2]]))
    if (!all(isMat)) stop('Every element in the list "mat" must be a numeric matrix/array')
    if (!all(eqDim)) stop('Every element in the list "mat" must be a square matrix,',
                          ' or the first 2 dimensions of each 3-D array must match')
    if (!all(eqLab)) stop('row and column names must match in each square matrix')
    multigroup <- TRUE
  } else {
    if (!(is.array(mat) && is.numeric(mat)))
      stop('"mat" must be a numeric matrix/array (or list of them)')
    if (diff(dim(mat)[1:2]) != 0)
      stop('"mat" must be a square matrix (or list of them), or the first 2',
           ' dimensions of a 3-D array must match')
    if ( !all.equal(dimnames(mat)[[1]], dimnames(mat)[[2]]) )
      stop('rownames(mat) must == colnames(mat)')
    mat <- list(mat)
    multigroup <- FALSE
  }
  if (!is.character(suffix) || length(suffix) != 2L)
    stop('"suffix" must be a length-2 character vector')

  outList <- sapply(mat, function(M) {
    ## save each pair of IDs
    RN <- rownames(M)
    if (is.null(RN)) RN <- 1:nrow(M)
    IDs <- combn(RN, m = 2, simplify = FALSE,
                 FUN = function(x) c(ID_A = x[1], ID_B = x[2]))

    ## save each dyadic (i.e., bivariate) outcome for each pair of IDs
    if (length(dim(M)) == 2L) {
      nYs <- 1L # 1 round-robin variable
      label <- as.character(label[1]) # use "label=" argument for variable names
      yLabels <- paste0(label, "_", suffix[1:2], suffix[2:1])

      ## (group-)mean center?
      if (groupMC) {
        Mcopy <- M
        diag(Mcopy) <- NA # ignore any self-reports in diagonal
        M <- M - mean(Mcopy, na.rm = TRUE)
      }

      ## for each dyad (i != j)
      yList <- utils::combn(RN, m = 2, simplify = FALSE,
                            FUN = function(x) c(AB = M[ x[1] , x[2] ],
                                                BA = M[ x[2] , x[1] ]))
      ## check for missing?
      if (dropNA) {
        bothNA <- sapply(yList, function(x) all(is.na(x)))
        IDs <- IDs[which(!bothNA)]
        yList <- yList[which(!bothNA)]
      }
      Ys <- do.call(rbind, yList)


      ## if multivariate, loop over third dimension
    } else if (length(dim(M)) == 3L) {
      nYs <- dim(M)[[3]] # how many round-robin variables?
      if (nYs == length(label)) {
        ## use "label=" argument for variable names
        yLabels <- paste0(rep(as.character(label), each = 2),
                          "_", suffix[1:2], suffix[2:1])

      } else if (!is.null(dimnames(M)[[3]])) {
        ## use list names to make variable names
        yLabels <- paste0(rep(dimnames(M)[[3]], each = 2),
                          "_", suffix[1:2], suffix[2:1])

      } else {
        ## use the first "label=" argument and add numbers
        yLabels <- paste0(as.character(label[1]), rep(1:nYs, each = 2),
                          "_", suffix[1:2], suffix[2:1])

      }

      yList <- list()
      if (is.null(dimnames(M)[[3]])) dimnames(M)[[3]] <- 1:(dim(M)[3])
      for (y in dimnames(M)[[3]]) {
        M2 <- M[ , , y, drop = TRUE] # extract slice for this RR variable

        ## (group-)mean center?
        if (groupMC) {
          Mcopy <- M2
          diag(Mcopy) <- NA # ignore any self-reports in diagonal
          M2 <- M2 - mean(Mcopy, na.rm = TRUE)
        }

        ## for each dyad (i != j)
        tempY <- combn(RN, m = 2, simplify = FALSE,
                       FUN = function(x) c(AB = M2[ x[1] , x[2] ],
                                           BA = M2[ x[2] , x[1] ]))
        yList[[y]] <- do.call(rbind, tempY)
      }
      Ys <- do.call(cbind, yList)

    } else stop('mat= argument must be a (list of) 2- or 3-dimensional array(s).')

    ## combine each dyad's information into a data.frame
    out <- data.frame(do.call(rbind, IDs), Ys, stringsAsFactors = FALSE)
    ## assemble requested variable names
    colnames(out)[1:2] <- paste("ID", suffix, sep = "_")
    colnames(out)[1:(2*nYs) + 2L] <- yLabels
    ## return result
    out
  }, simplify = FALSE)

  if (multigroup) {
    ## add group labels
    gNames <- names(mat)
    if (is.null(gNames)) gNames <- seq_along(mat)
    for (g in gNames) outList[[g]][ , group] <- g
  }

  result <- do.call(rbind, outList)
  if (multigroup) result <- result[ , c(group, setdiff(colnames(result), group))]
  rownames(result) <- NULL
  result
}


## TODO: extract case-level data from diagonal of each RR-variable matrix
square2case <- function(mat) {

}

## transform from square (adjacency matrix) to long-long (univariate) format
square2longUni <- function(mat, dropNA = TRUE) {

}


## transform from bivariate (long) to univariate (long-long) format
longBi2longUni <- function(Ys, data, ID1, ID2, IDdyad = NULL, group = NULL) {

}


## transform from univariate (long-long) to bivariate (long) format in 2 steps:
## longUni2square() then square2longBi()
longUni2longBi <- function(y, data, IDout, IDin, group = NULL, groupMC = FALSE,
                           returnList = TRUE, dropNA = FALSE,
                           label = NULL, suffix = c("i","j")) {
  ## convert to squares
  SQ <- longUni2square(y = y, data = data, IDout = IDout, IDin = IDin,
                       group = group, returnList = returnList)
  ## group ID lost, so recover the name from original data
  if (is.numeric(group)) group <- colnames(data)[group]
  ## convert to long (bivariate) format
  square2longBi(mat = SQ, dropNA = dropNA, suffix = suffix,
                group = group, groupMC = groupMC,
                label = if (is.null(label)) y else label)
}


