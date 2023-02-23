### Terrence D. Jorgensen
### Last updated: 23 February 2023
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
##'   `"group"` when group effects are modeled.  Ignored when `srm.param="mean"`
##'   because means are only a group-level statistic.
##'   Multiple options can be passed to `summary()` method, but `as.matrix()`
##'   only uses the first.
##' @param srm.param `character` specifying the type of SRM parameter to provide
##'   for `component`. Can be `%in% c("sd", "cor", "cov")`, as well as `"mean"`
##'   when `fixed.groups=FALSE`. The `summary()` method will always return the
##'   means as a group-level statistic, even for a single round-robin group.
##' @param posterior.est For `as.matrix()`, a length-1 `character` or `numeric`
##'   indicating which posterior summary statistic should be used as the point
##'   estimate.  Can be `%in% c("mean", "median", "mode", "min", "max", "hdi")`,
##'   or `c("EAP","MAP")` as alias for `c("mean","mode")`, respectively).
##'
##' * `posterior.est="mode"` calls [modeest::mlv()]. Pass additional arguments
##'   via `...` (e.g., to avoid warnings about default `method=`).
##' * `posterior.est="hdi"` calls [HDInterval::hdi()]. Pass additional arguments
##'   via `...` (e.g., to overwrite the default `credMass=.95`).
##'   Unlike any other `posterior.est=` options for `as.matrix()`, `"hdi"`
##'   returns two values per parameter: the lower and upper credible-interval
##'   limits.  This affects how the output is formatted.
##' * `posterior.est=numeric(1)` between 0 and 1 (inclusive) is accepted by
##'   `as.matrix()`, passed to [stats::quantile()] as the `probs=` argument;
##'   other arguments can be passed via `...` (e.g., `type=`).
##' * For the `summary()` method, `posterior.est` can be a vector with multiple
##'   options `%in% c("mean", "EAP", "median", "mode", "MAP")`, or can be set to
##'   `NULL` to return no point estimates (e.g., only `interval=` estimates).
##'
##' @param interval `character` indicating the type of uncertainty/credible
##'   interval estimate.  Can be `%in% c("central", "hdi")` or set to
##'   `NULL` to avoid returning interval estimates.
##' @param credMass `numeric` passed to [HDInterval::hdi()], or used to request
##'   `probs=` from [stats::quantile()] when `interval="central"`.
##' @param as.stanfit `logical` indicating whether to use a method (e.g.,
##'   `summary`) defined for a \code{\linkS4class{stanfit}} object.  Useful for
##'   obtaining diagnostics (e.g., R-hat, effective sample size) about MCMC
##'   results.  If `TRUE`, further arguments to
##'   \link[rstan]{summary,stanfit-method} can be passed via `...`
#TODO: @param from.stanfit (pars= via ...) to assign c("n_eff","Rhat") to rownames($summary)
##' @param meanstructure `logical` indicating whether the `vcov()` method
##'   includes posterior variability of the mean estimates (only when
##'   `component="group"`).
##' @param keep,drop Variables to include (`keep=`) or exclude (`drop=`) from
##'   the result.  `drop=` is ignored when `keep=` is specified.  For case- or
##'   dyad-level variables, `keep`/`drop` can be any combination of the original
##'   names of round-robin variables passed to [mvsrm()] or the names of
##'   specific components, which include `c("out","in")` suffixes for case-level
##'   components or `c("ij","ji")` suffixes for dyad/relationship-level
##'   components.
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
                   nobs     = "integer", # level-specific sample sizes
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


##' @importFrom methods getMethod
#TODO: add EAP and MAP as options
summary.mvSRM <- function(object, component = c("case","dyad"),
                          srm.param = c("sd","cor"),
                          posterior.est = "mean",
                          interval = "central", credMass = .95,
                          #TODO? to.data.frame = FALSE,
                          #TODO: rsquare = TRUE,
                          #TODO: from.stanfit = c("n_eff","Rhat"),
                          as.stanfit = FALSE, ...) {
  if (as.stanfit) {
    return(getMethod("summary", "stanfit")(object, ...))
  }

  component <- intersect(tolower(component), c("group", "case", "dyad"))
  if (!length(component)) stop('No valid choice of component= was found')

  srm.param <- intersect(tolower(srm.param), c("mean", "sd", "cor", "cov"))
  if (!length(srm.param)) stop('No valid choice of srm.param= was found')

  ## can be NULL
  posterior.est <- intersect(posterior.est,
                             c("mean", "median", "mode", "EAP", "MAP"))
  interval <- intersect(tolower(interval), c("central", "hdi"))

  output <- list()

  ## add mean structure to output?  Always in "group" component
  if ("mean" %in% srm.param) {
    if (length(posterior.est))  for (p in posterior.est) {
      output$group$mean[[p]] <- as.matrix.mvSRM(x = object, srm.param = "mean",
                                                posterior.est = p, ...)
    }
    if (length(interval)) {
      if ("central" %in% interval) {
        probs <- abs(0:1 - (1 - credMass)/2)
        lower <- as.matrix.mvSRM(object, srm.param = "mean",
                                 posterior.est = probs[1], ...)
        upper <- as.matrix.mvSRM(object, srm.param = "mean",
                                 posterior.est = probs[2], ...)
        output$group$mean$central <- rbind(lower = lower, upper = upper)
      }
      if ("hdi" %in% interval) {
        probs <- abs(0:1 - (1 - credMass)/2)
        output$group$mean$hdi <- as.matrix.mvSRM(object, srm.param = "mean",
                                                 posterior.est = "hdi",
                                                 credMass = credMass, ...)
      }
    }

    if (length(srm.param) == 1L) {
      ## nothing else to return
      return(output)

      ## otherwise, drop "mean" from srm.param=
    } else srm.param <- srm.param[-which(srm.param == "mean")]
  }


  ## loop over components to obtain their other SRM parameters
  for (COMP in component) {
    if (COMP == "group" && !any(object@parNames$sigma == "S_g")) {
      message('group level not modeled (IDgroup=NULL or fixed.groups=TRUE)')
      next
    }
    for (STAT in srm.param) {
      ## put SDs first
      if ("sd" %in% srm.param) {

        if (length(posterior.est))  for (p in posterior.est) {
          output[[COMP]]$sd[[p]] <- as.matrix.mvSRM(x = object,
                                                    srm.param = "sd",
                                                    component = COMP,
                                                    posterior.est = p)
        }
        ## add intervals?
        if (length(interval)) {
          if ("central" %in% interval) {
            probs <- abs(0:1 - (1 - credMass)/2)
            lower <- as.matrix.mvSRM(object, srm.param = "sd", component = COMP,
                                     posterior.est = probs[1], ...)
            upper <- as.matrix.mvSRM(object, srm.param = "sd", component = COMP,
                                     posterior.est = probs[2], ...)
            output[[COMP]]$sd$central <- rbind(lower = lower, upper = upper)
          }
          if ("hdi" %in% interval) {
            probs <- abs(0:1 - (1 - credMass)/2)
            output$group$sd$hdi <- as.matrix.mvSRM(object, srm.param = "sd",
                                                   posterior.est = "hdi",
                                                   credMass = credMass, ...)
          }
        }

      }

      ## then correlations
      if ("cor" %in% srm.param) {

        if (length(posterior.est))  for (p in posterior.est) {
          output[[COMP]]$cor[[p]] <- as.matrix.mvSRM(x = object,
                                                     srm.param = "cor",
                                                     component = COMP,
                                                     posterior.est = p)
        }
        ## add intervals?
        if (length(interval)) {
          if ("central" %in% interval) {
            probs <- abs(0:1 - (1 - credMass)/2)
            output[[COMP]]$cor$central$lower <- as.matrix.mvSRM(x = object,
                                                                srm.param = "cor",
                                                                component = COMP,
                                                                posterior.est = probs[1],
                                                                ...)
            output[[COMP]]$cor$central$upper <- as.matrix.mvSRM(x = object,
                                                                srm.param = "cor",
                                                                component = COMP,
                                                                posterior.est = probs[2],
                                                                ...)
          }
          if ("hdi" %in% interval) {
            probs <- abs(0:1 - (1 - credMass)/2)
            output$group$cor$hdi <- as.matrix.mvSRM(object, srm.param = "cor",
                                                    posterior.est = "hdi",
                                                    credMass = credMass, ...)
          }
        }

      }

      #TODO: add dyadic correlations to srm.param = c("inter","intra")
      #      when both, put in lower/upper.tri()

      ## then covariances (dyadic available?)
      if ("cov" %in% srm.param) {

        if (length(posterior.est))  for (p in posterior.est) {
          output[[COMP]]$cov[[p]] <- as.matrix.mvSRM(x = object,
                                                     srm.param = "cov",
                                                     component = COMP,
                                                     posterior.est = p)
        }
        ## add intervals?
        if (length(interval)) {
          if ("central" %in% interval) {
            probs <- abs(0:1 - (1 - credMass)/2)
            output[[COMP]]$cov$central$lower <- as.matrix.mvSRM(x = object,
                                                                srm.param = "cov",
                                                                component = COMP,
                                                                posterior.est = probs[1],
                                                                ...)
            output[[COMP]]$cov$central$upper <- as.matrix.mvSRM(x = object,
                                                                srm.param = "cov",
                                                                component = COMP,
                                                                posterior.est = probs[2],
                                                                ...)
          }
          if ("hdi" %in% interval) {
            probs <- abs(0:1 - (1 - credMass)/2)
            output$group$cov$hdi <- as.matrix.mvSRM(object, srm.param = "cov",
                                                    posterior.est = "hdi",
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


#TODO: add srm.param = "dyadic.recip" (vector), "intra.inter.cor" (2-tri matrix)
#TODO: optionally add.header=FALSE for summary() to set TRUE (its default)
##' @importFrom methods as
##' @importFrom stats quantile
##' @importFrom utils combn
as.matrix.mvSRM <- function(x, component, srm.param, posterior.est = "mean",
                            as.stanfit = FALSE, ...) {
  if (as.stanfit) {
    return(as.matrix(as(object = x, Class = "stanfit", strict = TRUE), ...))
  }

  if (missing(srm.param))
    stop('srm.param= argument must be specified. See class?mvSRM ',
         'help page for descriptions.')
  srm.param <- tolower(srm.param[1])
  stopifnot(srm.param %in% c("mean", "sd", "cor", "cov"))

  posterior.est <- posterior.est[1]
  if (is.numeric(posterior.est)) {
    stopifnot(posterior.est >= 0 && posterior.est <= 1)
  } else if (tolower(posterior.est) == "min") {
    posterior.est <- 0
  } else if (tolower(posterior.est) == "median") {
    posterior.est <- .5
  } else if (tolower(posterior.est) == "max") {
    posterior.est <- 1
  } else {
    posterior.est <- tolower(posterior.est)
    stopifnot(posterior.est %in% c("mean", "eap", "mode", "map", "hdi"))

    if (posterior.est %in% c("mode","map")) {
      if (requireNamespace("modeest")) {
        if (!"package:modeest" %in% search()) attachNamespace("modeest")
      } else stop('Install package `modeest` to estimate the posterior mode.')
      posterior.est <- "mlv"

    } else if (posterior.est == "hdi") {
      if (requireNamespace("HDInterval")) {
        if (!"package:HDInterval" %in% search()) attachNamespace("HDInterval")
      } else stop('Install package `HDInterval` to obtain posterior.est="HDI".')

    } else if (posterior.est %in% c("mean","eap")) posterior.est <- "mean"

  }

  if (srm.param == "mean") {
    if (!"mu" %in% names(x@parNames))
      stop('Means were not included in this model (i.e., fixed.groups = TRUE)')

    ## extract samples
    paramMat <- do.call(rbind, rstan::As.mcmc.list(x, pars = "Mvec"))
    ## obtain summary statistic
    if (is.numeric(posterior.est)) {
      out <- apply(paramMat, MARGIN = 2, stats::quantile,
                   probs = posterior.est, ...)
    } else out <- apply(paramMat, MARGIN = 2, posterior.est, ...)

    #TODO: c() additional names of case/group_data or dyad-constant variables
    if (posterior.est == "hdi") {
      colnames(out) <- names(x@varNames$RR)
    } else names(out) <- names(x@varNames$RR)

    ## nothing else to do for means
    return(out)


  } else if (missing(component)) {
    stop('Unless srm.param="mean", the component= argument must be specified. ',
         'See class?mvSRM help page for descriptions.')
  }
  component <- tolower(component[1])
  stopifnot(component %in% c("group", "case", "dyad"))


  if (component == "group") {
    if (!any(x@parNames$sigma == "S_g"))
      stop('group level not modeled (IDgroup=NULL or fixed.groups=TRUE)')
    PARS <- ifelse(srm.param == "sd", "S_g", "Rg")
    NAMES <- names(x@varNames$RR)

  } else if (component == "case") {
    PARS <- ifelse(srm.param == "sd", "S_p", "Rp")
    NAMES <- paste(rep(names(x@varNames$RR), each = 2),
                   c("out","in"), sep = "_")

  } else if (component == "dyad") {
    PARS <- ifelse(srm.param == "sd", "s_rr", "Rd2")
    if (srm.param == "sd") {
      NAMES <- names(x@varNames$RR)
    } else {
      ## extended names for "cor" or "cov"
      NAMES <- paste(rep(names(x@varNames$RR), each = 2),
                     c("ij","ji"), sep = "_")
    }

  }


  if (srm.param %in% c("sd", "cor")) {
    ## extract samples
    paramMat <- do.call(rbind, rstan::As.mcmc.list(x, pars = PARS))

    ## obtain summary statistic
    if (is.numeric(posterior.est)) {
      out <- apply(paramMat, MARGIN = 2, stats::quantile,
                   probs = posterior.est, ...)
    } else out <- apply(paramMat, MARGIN = 2, posterior.est, ...)

    if (srm.param == "sd") {
      ## attach variable names
      if (posterior.est == "hdi") {
        colnames(out)   <- NAMES #TODO: c() case/group_data names
      } else names(out) <- NAMES
      return(out) # nothing else to do for SDs

    } else if (posterior.est == "hdi") {
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
  ## else srm.param == "cov"

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
  if (posterior.est == "mean") return(Sigma)


  ## otherwise, get posterior.est= for each cell in Sigma (copy to retain class)
  if (posterior.est == "hdi") {
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
      if (is.numeric(posterior.est)) {
        SigmaSummaryStat[RR, RR] <- quantile(sigmaVec, probs = posterior.est, ...)
      } else SigmaSummaryStat[RR, RR] <- eval(as.name(posterior.est))(sigmaVec, ...)
    }
    ## loop over (pairs of) variables
    PAIRS <- combn(NAMES, m = 2)
    for (RR in PAIRS[1,]) {
      for (CC in PAIRS[2,]) {
        if (is.numeric(posterior.est)) {
          SigmaSummaryStat[RR, CC] <- quantile(sigmaVec, probs = posterior.est, ...)
          SigmaSummaryStat[CC, RR] <- SigmaSummaryStat[RR, CC]
        } else {
          SigmaSummaryStat[RR, CC] <- eval(as.name(posterior.est))(sigmaVec, ...)
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



##' @importFrom stats cov
vcov.mvSRM <- function(object, component, meanstructure = FALSE,
                       keep, drop, add.names.attr = FALSE, ...) {
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

    PARS <- c("S_g", "Rg")
    NAMES <- names(object@varNames$RR) #TODO: add covariates

    if (!missing(keep)) {
      SUBSET <- intersect(keep, NAMES)
      if (!length(SUBSET)) stop('keep=', keep, ' leaves no variables to return')
    } else if (!missing(drop)) {
      SUBSET <- setdiff(NAMES, drop)
      if (!length(SUBSET)) stop('drop=', drop, ' leaves no variables to return')
    } else SUBSET <- NAMES

    if (meanstructure) {
      MuSamples <- do.call(rbind, rstan::As.mcmc.list(object, pars = "Mvec"))
      MuList <- apply(MuSamples, MARGIN = 1, FUN = function(m) {
        Mvec <- setNames(numeric(length(NAMES)), nm = NAMES)
        for (i in names(m)) eval(parse(text = paste(i, "<-", m[i]) ))
        Mvec[SUBSET]
      }, simplify = FALSE)

    }


  } else if (component == "case") {

    PARS <- c("S_p", "Rp")
    NAMES <- paste(rep(names(object@varNames$RR), each = 2),
                   c("out","in"), sep = "_")

    if (!missing(keep)) {
      KEEP <- do.call(c, sapply(keep, function(k) {
        if (k %in% names(object@varNames$RR)) {
          return(paste0(k, c("_out", "_in")))
        } else if (k %in% object@varNames$case) {
          return(k)
        }
        NULL
      }, simplify = FALSE))
      SUBSET <- intersect(KEEP, NAMES) #FIXME: necessary? KEEP should be SUBSET
      if (!length(SUBSET)) stop('keep=', keep, ' leaves no variables to return')

    } else if (!missing(drop)) {
      DROP <- do.call(c, sapply(drop, function(d) {
        if (d %in% names(object@varNames$RR)) {
          return(paste0(d, c("_out", "_in")))
        } else if (d %in% object@varNames$case) {
          return(d)
        }
        NULL
      }, simplify = FALSE))
      SUBSET <- setdiff(NAMES, DROP)
      if (!length(SUBSET)) stop('drop=', drop, ' leaves no variables to return')

    } else SUBSET <- NAMES


  } else if (component == "dyad") {

    PARS <- c("s_rr", "Rd2")
    NAMES <- paste(rep(names(object@varNames$RR), each = 2),
                   c("ij","ji"), sep = "_")

    if (!missing(keep)) {
      KEEP <- do.call(c, sapply(keep, function(k) {
        if (k %in% names(object@varNames$RR)) {
          return(paste0(k, c("_ij", "_ji")))
        } else if (k %in% do.call(c, object@varNames$RR)) {
          return(k)
        }
        NULL
      }, simplify = FALSE))
      SUBSET <- intersect(KEEP, NAMES) #FIXME: necessary? KEEP should be SUBSET
      if (!length(SUBSET)) stop('keep=', keep, ' leaves no variables to return')

    } else if (!missing(drop)) {
      DROP <- do.call(c, sapply(drop, function(d) {
        if (d %in% names(object@varNames$RR)) {
          return(paste0(d, c("_ij", "_ji")))
        } else if (d %in% do.call(c, object@varNames$RR)) {
          return(d)
        }
        NULL
      }, simplify = FALSE))
      SUBSET <- setdiff(NAMES, DROP)
      if (!length(SUBSET)) stop('drop=', drop, ' leaves no variables to return')

    } else SUBSET <- NAMES
  }


  ## store correlations in a matrix (per posterior sample)
  corMat <- do.call(rbind, rstan::As.mcmc.list(object, pars = PARS[2]))
  corList <- apply(corMat, MARGIN = 1, FUN = function(m) {
    assign(PARS[2], matrix(0, nrow = length(NAMES), ncol = length(NAMES),
                           dimnames = list(NAMES, NAMES)))
    for (i in names(m)) eval(parse(text = paste(i, "<-", m[i]) ))
    return( eval(as.name(PARS[2]))[SUBSET, SUBSET] )
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
    SD[SUBSET, SUBSET]
  }, simplify = FALSE)

  ## scale correlations to covariance matrices
  SigmaList <- mapply(function(R, SD) SD %*% R %*% SD,
                      R = corList, S = sdList, SIMPLIFY = FALSE)

  if (!exists("MuList")) MuList <- list(NULL)
  ## posterior variability
  postVar <- mapply(function(S, M = NULL) {
    if (categorical) {
      wls.obs <- c(M, # 1. interleaved thresholds + (negative) means, if any
                   # 2. slopes (if any)
                   # 3. variances (of continuous variables, if any)
                   as.numeric(diag(S)),
                   # 4. lower.tri covariances/correlations
                   lavaan::lav_matrix_vech(S, diagonal = FALSE))
    } else {
       wls.obs <- c(M, # 1. means (if any)
                    # 2. slopes (if any)
                    # 3. lower.tri (co)variance matrix
                    lavaan::lav_matrix_vech(S, diagonal = TRUE))
    }
    wls.obs
  }, S = SigmaList, M = MuList, SIMPLIFY = FALSE)
  #FIXME: multiply ACOV by N instead of N-1 ?
  NACOV <- (object@nobs[[component]] - 1L) * cov(do.call(rbind, postVar))
  class(NACOV) <- c("lavaan.matrix.symmetric", "matrix", "array")
  if (add.names.attr) attr(NACOV, "subset") <- SUBSET # used in srm2lavData()
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
#      coef = call as.matrix(..., posterior.est = c("mean","median","mode"))
#      confint = central or HDI


#TODO: as.lavMoments, or define inherting lavSRMoments class?
#      or just srm2lavData(), no real benefit of publicizing class/methods

#TODO? as.stanfit() to enable using stanfit-class methods
#      Methods already work, unless I defined new ones (e.g., summary).
#      Better to keep including as.stanfit= argument to internally getMethod()


