### Terrence D. Jorgensen
### Last updated: 13 February 2023
### Class and Methods for mvSRM object


##' Class for a Multivariate SRM using Stan
##'
##' This class extends the \code{\linkS4class{stanfit}} class, created by
##' fitting a multivariate social relations model (mvSRM; Nestler,
##' [2018](https://doi.org/10.3102/1076998617741106)) using the Stan
##' (Carpenter et al., [2017](https://doi.org/10.18637/jss.v076.i01)) software.
##'
##' @name mvSRM-class
##' @importClassesFrom rstan stanfit
##' @aliases mvSRM-class show,mvSRM-method
#   bayes_R2,mvSRM-method
#   posterior_interval,mvSRM-method
##' @docType class
##'
##' @slot call A copy of the `call` to the [mvsrm()] function that generated
##'   the object.
##' @slot nobs A vector of group-level (when relevant), case-level, and
##'   dyad-level sample sizes.
#TODO: make case-level nested for blocks?  Can be random across RR groups...
##' @slot IDs `list` of `data.frame`s containing ID variables from data sets at
##'   each level of analysis.  Stored as integers with [attributes()] that are
##'   the original labels in the user's `data=` (and optionally `case_data=`
##'   and `group_data=`) passed to [mvsrm()].
##' @slot varNames `list` of variable names at each level of analysis.  This is
##'   used to ease the writing of methods to extract information by variable,
##'   since each variable has 2 person-level and 2 dyad-level components.
##' @slot parNames `list` of parameters, organized in sets:
##'   - `$mu` contains `"Mvec"` (the mean vector, if estimated)
##'   - `$sigma` contains `c("s_rr","S_p")` (also `"S_g"`, if estimated). These
##'     are the *SD*s of relationship-level and case-level components of
##'     round-robin variables (and *SD*s of covariates at those levels), as well
##'     as group-level components when relevant.
##'   - `$corr` contains `c("Rd2","Rp")` (also `"Rg"`, if estimated). These are
##'     correlation matrices among relationship-level and case-level components
##'     of round-robin variables (and correlations with covariates at those
##'     levels), as well as group-level components when relevant.
##'   - `$derived` contains `"Rsq"`: a matrix with 3-4 columns and 1 row per
##'     round-robin variable listed in `object@varNames$RR`. When group-level
##'     variances are estimated, the first of 4 columns contains the proportion
##'     of variance (in each round-robin variable) accounted for by group-mean
##'     differences.  The remaining columns are proportions of variance
##'     accounted for by outgoing (e.g., actor, perceiver, ego) effects, by
##'     incoming (e.g., partner, target, alter) effects, and by relationship
##'     effects (confounded with any dyad-level measurement error).
##'   - `$components` contains the relationship-level and case-level components
##'     themselves when `object@call$saveComp=TRUE`. This makes the object take
##'     much more computer memory, but it would enable the plausible-values
##'     approach described by Lüdtke et al. (2018).
##'
##' @param object,x An object of class `mvSRM`
##' @param ... Additional arguments passed to or from other methods
##' @param component `character` specifying which SRM component to report
##'   estimated parameters for. Can be `%in% c("case", "dyad")`, as well as
##'   `"group"` when group effects are modeled.  Ignored when `stat="mean"`
##'   because means are only a group-level statistic.
##'   Multiple options can be passed to `summary()` method, but `as.matrix()`
##'   only uses the first.
##' @param stat `character` specifying which summary statistic to report for
##'   `component`. Can be `%in% c("sd", "cor", "cov")`, as well as `"mean"` when
##'   `fixed.groups=FALSE`. The `summary()` method will always return the
##'   means as a group-level statistic, even for a single round-robin group.
##' @param point For `as.matrix()`, a length-1 `character` or `numeric`
##'   indicating which posterior summary statistic should be returned.  Can be
##'   `%in% c("mean", "median", "mode", "min", "max", "hdi")`.
##'
##' * `point="mode"` calls [modeest::mlv()], whose additional arguments can be
##'   passed via `...` (recommended to avoid warnings about default `method=`).
##' * `point="hdi"` calls [HDInterval::hdi()], whose additional arguments can
##'   be passed via `...` (e.g., to overwrite the default `credMass=.95`).
##'   Unlike any other `point=` options for `as.matrix()`, `"hdi"` returns two
##'   values per parameter: the lower and upper credible-interval limits.
##'   This affects how the output is formatted.
##' * `point=numeric(1)` between 0 and 1 (inclusive) is accepted by
##'   `as.matrix()`, passed to [stats::quantile()] as the `probs=` argument;
##'   other arguments can be passed via `...`.
##' * For the `summary()` method, multiple options can be passed, but must be
##'   `%in% c("mean", "median", "mode")` or set to `NULL` to avoid returning
##'   point estimates.
##'
##' @param interval `character` indicating the type of uncertainty/credible
##'   interval estimate.  Can be `%in% c("central", "hdi")` or set to
##'   `NULL` to avoid returning interval estimates.
##' @param credMass `numeric` passed to [HDInterval::hdi()], or used to request
##'   `probs=` from [stats::quantile()] when `interval="central"`.
##' @param meanstructure `logical` indicating whether the `vcov()` method
##'   includes posterior variability of the mean estimates (only when
##'   `component="group"`).
##'
##' @return
##' * `show()` simply prints a message stating the class and inheritance.
##' * `as.matrix()` returns a vector or matrix of the requested posterior
##'   summary for the specified SRM parameter and level of analysis.
##' * `summary()` loops over `as.matrix()` to return a `list()` of potentially
##'   multiple summaries of potentially multiple SRM parameters and components.
##' * `vcov()` provides the posterior covariance matrix of sampled (co)variance
##'   parameters at the requested level of analysis.  This is only useful for
##'   the 2-stage estimation provided using [lavaan::lavaan()].
##'
##' @section Objects from the Class: See the [mvsrm()] function for details.
##'
##' @author Terrence D. Jorgensen (University of Amsterdam;
##'   \email{TJorgensen314@@gmail.com})
##'
##' @references
##'  Carpenter, B., Gelman, A., Hoffman, M. D., Lee, D., Goodrich, B.,
##'  Betancourt, M., ... & Riddell, A. (2017).
##'  Stan: A probabilistic programming language.
##'  *Journal of Statistical Software, 76*(1).
##'  <https://doi.org/10.18637/jss.v076.i01>
##'
##'  Lüdtke, O., Robitzsch, A., & Trautwein, U. (2018). Integrating covariates
##'  into social relations models: A plausible values approach for handling
##'  measurement error in perceiver and target effects.
##'  *Multivariate Behavioral Research, 53*(1), 102--124.
##'  <https://doi.org/10.1080/00273171.2017.140679>
##'
##'  Nestler, S. (2018). Likelihood estimation of the multivariate social
##'  relations model. *Journal of Educational and Behavioral Statistics, 43*(4),
##'  387--406. <https://doi.org/10.3102/1076998617741106>
##'
##' @seealso [mvsrm()]
##'
## @examples
##'
##' @importFrom methods setClass
##' @export
setClass("mvSRM", contains = "stanfit",
         slots = c(call     = "call",    # mvsrm() call
                   N        = "integer", # level-specific sample sizes
                   IDs      = "list",    # ID variables in each level's data set
                   varNames = "list",    # variable names at each level
                   parNames = "list"))   # names of sets of estimated parameters



##' @name mvSRM-class
##' @aliases show,mvSRM-method
##' @importFrom methods show
##' @export
setMethod("show", "mvSRM", function(object) {
  cat('This object contains results of a multivariate SRM fit using (r)Stan.',
      'See class?stanfit help page for general methods for rstan output,',
      'and class?mvSRM for additional methods relevant for a mvSRM model.',
      sep = "\n")
  invisible(object)
})


## @importFrom?
summary.mvSRM <- function(object, component = c("case","dyad"),
                          stat = c("sd","cor"), point = "mean",
                          interval = "central", credMass = .95,
                          #TODO? to.data.frame = FALSE,
                          #TODO: rsquare = TRUE,
                          ...) {
  component <- intersect(tolower(component), c("group", "case", "dyad"))
  if (!length(component)) stop('No valid choice of component= was found')

  stat <- intersect(tolower(stat), c("mean", "sd", "cor", "cov"))
  if (!length(stat)) stop('No valid choice of stat= was found')

  ## can be NULL
  point <- intersect(tolower(point), c("mean", "median", "mode"))
  interval <- intersect(tolower(interval), c("central", "hdi"))

  output <- list()

  ## add mean structure to output?  Always in "group" component
  if ("mean" %in% stat) {
    if (length(point))  for (p in point) {
      output$group$mean[[p]] <- as.matrix.mvSRM(x = object, stat = "mean",
                                                point = p, ...)
    }
    if (length(interval)) {
      if ("central" %in% interval) {
        probs <- abs(0:1 - (1 - credMass)/2)
        lower <- as.matrix.mvSRM(object, stat = "mean", point = probs[1], ...)
        upper <- as.matrix.mvSRM(object, stat = "mean", point = probs[2], ...)
        output$group$mean$central <- rbind(lower = lower, upper = upper)
      }
      if ("hdi" %in% interval) {
        probs <- abs(0:1 - (1 - credMass)/2)
        output$group$mean$hdi <- as.matrix.mvSRM(object, stat = "mean",
                                                 point = "hdi",
                                                 credMass = credMass, ...)
      }
    }

    if (length(stat) == 1L) {
      ## nothing else to return
      return(output)

      ## otherwise, drop "mean" from stat=
    } else stat <- stat[-which(stat == "mean")]
  }


  ## loop over components for other SRM stats
  for (COMP in component) {
    if (COMP == "group" && !any(object@parNames$sigma == "S_g")) {
      message('group level not modeled (IDgroup=NULL or fixed.groups=TRUE)')
      next
    }
    for (STAT in stat) {
      ## put SDs first
      if ("sd" %in% stat) {

        if (length(point))  for (p in point) {
          output[[COMP]]$sd[[p]] <- as.matrix.mvSRM(x = object, stat = "sd",
                                                    component = COMP, point = p)
        }
        ## add intervals?
        if (length(interval)) {
          if ("central" %in% interval) {
            probs <- abs(0:1 - (1 - credMass)/2)
            lower <- as.matrix.mvSRM(object, stat = "sd", component = COMP,
                                     point = probs[1], ...)
            upper <- as.matrix.mvSRM(object, stat = "sd", component = COMP,
                                     point = probs[2], ...)
            output[[COMP]]$sd$central <- rbind(lower = lower, upper = upper)
          }
          if ("hdi" %in% interval) {
            probs <- abs(0:1 - (1 - credMass)/2)
            output$group$sd$hdi <- as.matrix.mvSRM(object, stat = "sd",
                                                   point = "hdi",
                                                   credMass = credMass, ...)
          }
        }

      }

      ## then correlations
      if ("cor" %in% stat) {

        if (length(point))  for (p in point) {
          output[[COMP]]$cor[[p]] <- as.matrix.mvSRM(x = object, stat = "cor",
                                                     component = COMP, point = p)
        }
        ## add intervals?
        if (length(interval)) {
          if ("central" %in% interval) {
            probs <- abs(0:1 - (1 - credMass)/2)
            output[[COMP]]$cor$central$lower <- as.matrix.mvSRM(x = object,
                                                                stat = "cor",
                                                                component = COMP,
                                                                point = probs[1],
                                                                ...)
            output[[COMP]]$cor$central$upper <- as.matrix.mvSRM(x = object,
                                                                stat = "cor",
                                                                component = COMP,
                                                                point = probs[2],
                                                                ...)
          }
          if ("hdi" %in% interval) {
            probs <- abs(0:1 - (1 - credMass)/2)
            output$group$cor$hdi <- as.matrix.mvSRM(object, stat = "cor",
                                                    point = "hdi",
                                                    credMass = credMass, ...)
          }
        }

      }

      #TODO: add dyadic correlations to stat = c("inter","intra")
      #      when both, put in lower/upper.tri()

      ## then covariances (dyadic available?)
      if ("cov" %in% stat) {

        if (length(point))  for (p in point) {
          output[[COMP]]$cov[[p]] <- as.matrix.mvSRM(x = object, stat = "cov",
                                                     component = COMP, point = p)
        }
        ## add intervals?
        if (length(interval)) {
          if ("central" %in% interval) {
            probs <- abs(0:1 - (1 - credMass)/2)
            output[[COMP]]$cov$central$lower <- as.matrix.mvSRM(x = object,
                                                                stat = "cov",
                                                                component = COMP,
                                                                point = probs[1],
                                                                ...)
            output[[COMP]]$cov$central$upper <- as.matrix.mvSRM(x = object,
                                                                stat = "cov",
                                                                component = COMP,
                                                                point = probs[2],
                                                                ...)
          }
          if ("hdi" %in% interval) {
            probs <- abs(0:1 - (1 - credMass)/2)
            output$group$cov$hdi <- as.matrix.mvSRM(object, stat = "cov",
                                                    point = "hdi",
                                                    credMass = credMass, ...)
          }
        }

      }

    } # end STAT

  } # end COMP


  #TODO: if (rsquare) {}
  #      get R-squared from bayes_R2() method

  output
}
##' @name mvSRM-class
##' @aliases summary,mvSRM-method
##' @export
setMethod("summary", "mvSRM", summary.mvSRM)


#TODO: add stat = "dyadic.recip" (vector), "intra.inter.cor" (2-tri matrix)
#TODO: optionally add.header=FALSE for summary() to set TRUE (its default)
##' @importFrom stats quantile
##' @importFrom utils combn
as.matrix.mvSRM <- function(x, stat, component, point = "mean", ...) {
  if (missing(stat))
    stop('stat= argument must be specified. See class?mvSRM ',
         'help page for descriptions.')
  stat      <- tolower(stat[1])
  stopifnot(stat %in% c("mean", "sd", "cor", "cov"))

  point <- point[1]
  if (is.numeric(point)) {
    stopifnot(point >= 0 && point <= 1)
  } else if (tolower(point) == "min") {
    point <- 0
  } else if (tolower(point) == "median") {
    point <- .5
  } else if (tolower(point) == "max") {
    point <- 1
  } else {
    point <- tolower(point)
    stopifnot(point %in% c("mean", "mode", "hdi"))

    if (point == "mode") {
      if (requireNamespace("modeest")) {
        if (!"package:modeest" %in% search()) attachNamespace("modeest")
      } else stop('Install the modeest package to obtain point="mode".')
      point <- "mlv"

    } else if (point == "hdi") {
      if (requireNamespace("HDInterval")) {
        if (!"package:HDInterval" %in% search()) attachNamespace("HDInterval")
      } else stop('Install the HDInterval package to obtain point="HDI".')
    }

  }

  if (stat == "mean") {
    if (!"mu" %in% names(x@parNames))
      stop('Means were not included in this model (i.e., fixed.groups = TRUE)')

    ## extract samples
    paramMat <- do.call(rbind, rstan::As.mcmc.list(x, pars = "Mvec"))
    ## obtain summary statistic
    if (is.numeric(point)) {
      out <- apply(paramMat, MARGIN = 2, stats::quantile,
                   probs = point, ...)
    } else out <- apply(paramMat, MARGIN = 2, point, ...)

    #TODO: c() additional names of case/group_data or dyad-constant variables
    if (point == "hdi") {
      colnames(out) <- names(x@varNames$RR)
    } else names(out) <- names(x@varNames$RR)

    ## nothing else to do for means
    return(out)


  } else if (missing(component)) {
    stop('Unless stat="mean", the component= argument must be specified. ',
         'See class?mvSRM help page for descriptions.')
  }
  component <- tolower(component[1])
  stopifnot(component %in% c("group", "case", "dyad"))


  if (component == "group") {
    if (!any(x@parNames$sigma == "S_g"))
      stop('group level not modeled (IDgroup=NULL or fixed.groups=TRUE)')
    PARS <- ifelse(stat == "sd", "S_g", "Rg")
    NAMES <- names(x@varNames$RR)

  } else if (component == "case") {
    PARS <- ifelse(stat == "sd", "S_p", "Rp")
    NAMES <- paste(rep(names(x@varNames$RR), each = 2),
                   c("out","in"), sep = "_")

  } else if (component == "dyad") {
    PARS <- ifelse(stat == "sd", "s_rr", "Rd2")
    if (stat == "sd") {
      NAMES <- names(x@varNames$RR)
    } else {
      ## extended names for "cor" or "cov"
      NAMES <- paste(rep(names(x@varNames$RR), each = 2),
                     c("ij","ji"), sep = "_")
    }

  }


  if (stat %in% c("sd", "cor")) {
    ## extract samples
    paramMat <- do.call(rbind, rstan::As.mcmc.list(x, pars = PARS))

    ## obtain summary statistic
    if (is.numeric(point)) {
      out <- apply(paramMat, MARGIN = 2, stats::quantile,
                   probs = point, ...)
    } else out <- apply(paramMat, MARGIN = 2, point, ...)

    if (stat == "sd") {
      ## attach variable names
      if (point == "hdi") {
        colnames(out)   <- NAMES #TODO: c() case/group_data names
      } else names(out) <- NAMES
      return(out) # nothing else to do for SDs

    } else if (point == "hdi") {
      lower <- out[1,]
      upper <- out[2,]
      ## store correlations in a list of 2 matrices ($lower and $upper)
      assign(PARS, matrix(0, nrow = length(NAMES), ncol = length(NAMES),
                          dimnames = list(NAMES, NAMES)))
      for (i in names(lower)) {
        eval(parse(text = paste(i, "<-", lower[i]) ))
      }
      LOWER <- eval(as.name(PARS))
      class(LOWER) <- c("lavaan.matrix.symmetric","matrix")
      for (i in names(upper)) {
        eval(parse(text = paste(i, "<-", upper[i]) ))
      }
      UPPER <- eval(as.name(PARS))
      class(UPPER) <- c("lavaan.matrix.symmetric","matrix")
      ## nothing else to do for correlation limits
      return(list(lower = LOWER, upper = UPPER))

    } else {
      ## store correlations in a matrix
      assign(PARS, matrix(0, nrow = length(NAMES), ncol = length(NAMES),
                          dimnames = list(NAMES, NAMES)))
      for (i in names(out)) {
        eval(parse(text = paste(i, "<-", out[i]) ))
      }
      out <- eval(as.name(PARS))
      class(out) <- c("lavaan.matrix.symmetric","matrix")

      return(out) # nothing else to do for correlations
    }

  }
  ## else stat == "cov"

  ## scale R with SDs per iteration, then assemble matrix

  ## extract samples of both SDs and correlations
  if (component == "group") {
    PARS <- c("S_g", "Rg")
  } else if (component == "case") {
    PARS <- c("S_p", "Rp")
  } else if (component == "dyad") {
    PARS <- c("s_rr", "Rd2")
  }
  sdMat  <- do.call(rbind, rstan::As.mcmc.list(x, pars = PARS[1]))
  corMat <- do.call(rbind, rstan::As.mcmc.list(x, pars = PARS[2]))


  ## store correlations in a matrix (per posterior sample)
  corList <- apply(corMat, MARGIN = 1, FUN = function(m) {
    assign(PARS[2], matrix(0, nrow = length(NAMES), ncol = length(NAMES),
                           dimnames = list(NAMES, NAMES)))
    for (i in names(m)) eval(parse(text = paste(i, "<-", m[i]) ))
    return( eval(as.name(PARS[2])) )
  }, simplify = FALSE)

  ## store SDs in a diagonal matrix (per posterior sample)
  sdList <- apply(sdMat, MARGIN = 1, FUN = function(m) {
    if (PARS[1] == "s_rr") {
           assign(PARS[1], numeric(length(NAMES)/2)) # equality constraints
    } else assign(PARS[1], numeric(length(NAMES)  ))

    for (i in names(m)) eval(parse(text = paste(i, "<-", m[i]) ))

    if (PARS[1] == "s_rr") {
      SD <- diag(rep(eval(as.name(PARS[1])), each = 2)) # repeat equal variances
    } else SD <- diag( eval(as.name(PARS[1])) )

    dimnames(SD) <- list(NAMES, NAMES)
    SD
  }, simplify = FALSE)

  ## scale correlations to covariance matrices
  SigmaList <- mapply(function(R, SD) SD %*% R %*% SD,
                       R = corList, S = sdList, SIMPLIFY = FALSE)

  Sigma <- Reduce("+", SigmaList) / length(SigmaList)
  class(Sigma) <- c("lavaan.matrix.symmetric","matrix")
  if (point == "mean") return(Sigma)


  ## otherwise, get point= for each cell in Sigma (copy to retain class)
  if (point == "hdi") {
    SigmaSummaryStat <- list(lower = Sigma, upper = Sigma)
    ## obtain diagonal entries first
    for (RR in NAMES) {
      sigmaHDI <- HDInterval::hdi(sapply(SigmaList, "[", i = RR, j = RR), ...)
      SigmaSummaryStat$lower[RR, RR] <- sigmaHDI[["lower"]]
      SigmaSummaryStat$upper[RR, RR] <- sigmaHDI[["upper"]]
    }
    ## loop over (pairs of) variables
    PAIRS <- combn(NAMES, m = 2)
    for (RR in PAIRS[1,]) {
      for (CC in PAIRS[2,]) {
        sigmaHDI <- HDInterval::hdi(sapply(SigmaList, "[", i = RR, j = CC), ...)

        SigmaSummaryStat$lower[RR, CC] <- sigmaHDI[["lower"]]
        SigmaSummaryStat$lower[CC, RR] <- sigmaHDI[["lower"]]

        SigmaSummaryStat$upper[RR, CC] <- sigmaHDI[["upper"]]
        SigmaSummaryStat$upper[CC, RR] <- sigmaHDI[["upper"]]
      }
    }

  } else {
    ## only scalars returned
    SigmaSummaryStat <- Sigma
    ## obtain diagonal entries first
    for (RR in NAMES) {
      sigmaVec <- sapply(SigmaList, "[", i = RR, j = RR)
      if (is.numeric(point)) {
        SigmaSummaryStat[RR, RR] <- quantile(sigmaVec, probs = point, ...)
      } else SigmaSummaryStat[RR, RR] <- eval(as.name(point))(sigmaVec, ...)
    }
    ## loop over (pairs of) variables
    PAIRS <- combn(NAMES, m = 2)
    for (RR in PAIRS[1,]) {
      for (CC in PAIRS[2,]) {
        if (is.numeric(point)) {
          SigmaSummaryStat[RR, CC] <- quantile(sigmaVec, probs = point, ...)
          SigmaSummaryStat[CC, RR] <- SigmaSummaryStat[RR, CC]
        } else {
          SigmaSummaryStat[RR, CC] <- eval(as.name(point))(sigmaVec, ...)
          SigmaSummaryStat[CC, RR] <- SigmaSummaryStat[RR, CC]
        }
      }
    }

  }

  SigmaSummaryStat
}
##' @name mvSRM-class
##' @aliases as.matrix,mvSRM-method
##' @export
setMethod("as.matrix", "mvSRM", as.matrix.mvSRM)


#TODO: vcov method to assemble uncertainty?
#      This might not be useful for users, just to make lavMoments.
#      But could be used for Wald tests based on MCMC-estimated ACOV.
##' @importFrom stats cov
vcov.mvSRM <- function(object, component, meanstructure = FALSE, ...) {
  categorical <- FALSE #TODO: add threshold models to stan scripts
  #TODO: robust options?  (e.g., Spearman rank cor, scaled by median abs dev)
  component <- tolower(component[1])
  stopifnot(component %in% c("group", "case", "dyad"))

  if (meanstructure) {
    if (isTRUE(object@call$fixed.groups)) stop('means not modeled when fixed.groups=TRUE')
    if (component != "group") stop('Means can only be modeled at the group level.')
  }

  if (component == "group") {
    if (!any(object@parNames$sigma == "S_g"))
      stop('group level not modeled (IDgroup=NULL or fixed.groups=TRUE)')

    NAMES <- names(object@varNames$RR)
    PARS <- c("S_g", "Rg")

  } else if (component == "case") {

    PARS <- c("S_p", "Rp")
    NAMES <- paste(rep(names(object@varNames$RR), each = 2),
                   c("out","in"), sep = "_")

  } else if (component == "dyad") {

    PARS <- c("s_rr", "Rd2")
    NAMES <- paste(rep(names(object@varNames$RR), each = 2),
                   c("ij","ji"), sep = "_")
  }


  ## store correlations in a matrix (per posterior sample)
  corMat <- do.call(rbind, rstan::As.mcmc.list(object, pars = PARS[2]))
  corList <- apply(corMat, MARGIN = 1, FUN = function(m) {
    assign(PARS[2], matrix(0, nrow = length(NAMES), ncol = length(NAMES),
                           dimnames = list(NAMES, NAMES)))
    for (i in names(m)) eval(parse(text = paste(i, "<-", m[i]) ))
    return( eval(as.name(PARS[2])) )
  }, simplify = FALSE)

  ## store SDs in a diagonal matrix (per posterior sample)
  sdMat  <- do.call(rbind, rstan::As.mcmc.list(object, pars = PARS[1]))
  sdList <- apply(sdMat, MARGIN = 1, FUN = function(m) {
    if (PARS[1] == "s_rr") {
      assign(PARS[1], numeric(length(NAMES)/2)) # equality constraints
    } else assign(PARS[1], numeric(length(NAMES)  ))

    for (i in names(m)) eval(parse(text = paste(i, "<-", m[i]) ))

    if (PARS[1] == "s_rr") {
      SD <- diag(rep(eval(as.name(PARS[1])), each = 2)) # repeat equal variances
    } else SD <- diag( eval(as.name(PARS[1])) )

    dimnames(SD) <- list(NAMES, NAMES)
    SD
  }, simplify = FALSE)

  ## scale correlations to covariance matrices
  SigmaList <- mapply(function(R, SD) SD %*% R %*% SD,
                      R = corList, S = sdList, SIMPLIFY = FALSE)

  ## posterior variability
  if (meanstructure) {
    #TODO

  }
  postVar <- lapply(SigmaList, function(x) {
    if (categorical) {
      wls.obs <- c(NULL, # 1. interleaved thresholds + (negative) means, if any
                   # 2. slopes (if any)
                   # 3. variances (of continuous variables, if any)
                   as.numeric(diag(x)),
                   # 4. lower.tri covariances/correlations
                   lavaan::lav_matrix_vech(x, diagonal = FALSE))
    } else {
       wls.obs <- c(NULL, # 1. means (if any)
                    # 2. slopes (if any)
                    # 3. lower.tri (co)variance matrix
                    lavaan::lav_matrix_vech(x, diagonal = TRUE))
    }
    wls.obs
  })
  NACOV <- (length(SigmaList) - 1L) * cov(do.call(rbind, postVar))
  class(NACOV) <- c("lavaan.matrix.symmetric", "matrix", "array")
  NACOV
}
##' @name mvSRM-class
##' @aliases vcov,mvSRM-method
##' @importFrom stats vcov
##' @export
setMethod("vcov", "mvSRM", vcov.mvSRM)



#TODO: @importFrom stats
#      nobs = vector of level-specific Ns
#      resid(uals) = random effects (or direct to semTools::plausibleValues)
#      fitted fitted.values
#      coef = call as.matrix(..., point = c("mean","median","mode"))
#      confint = central or HDI
#      vcov = (N)ACOV


#TODO: as.lavMoments, or define inherting lavSRMoments class?
#      or just stan2lavData(), no real benefit of publicizing class/methods

#TODO: as.stanfit() to enable using stanfit-class methods


