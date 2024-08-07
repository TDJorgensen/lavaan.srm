### Terrence D. Jorgensen
### Last updated: 11 July 2024
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
##' @aliases mvSRM-class
##'          show,mvSRM-method
##'          summary,mvSRM-method
##'          as.matrix,mvSRM-method
##'          vcov,mvSRM-method
#   bayes_R2,mvSRM-method
#   posterior_interval,mvSRM-method
##' @docType class
##'
##' @slot version Named `character` vector indicating the `stan`, `rstan`, and
##'   `lavaan.srm` version numbers.
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
##'   since each round-robin variable has 2 person-level and 2 dyad-level components.
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
##'     The dyad-level correlations are stored in a compact way to assist
##'     evaluating nonredundant elements: each variable's dyadic reciprocity is
##'     stored, and for each pair of variables, the inter/intra-personal
##'     correlations are above/below the diagonal, respectively.
##'   - `$derived` contains `"Rsq"`: a matrix with 3-4 columns and 1 row per
##'     round-robin variable listed in `object@varNames$RR`. When group-level
##'     variances are estimated, the first of 4 columns contains the proportion
##'     of variance (in each round-robin variable) accounted for by group-mean
##'     differences.  The remaining columns are proportions of variance
##'     accounted for by outgoing (e.g., actor, perceiver, ego) effects, by
##'     incoming (e.g., partner, target, alter) effects, and by relationship
##'     effects (confounded with any dyad-level measurement error).
##'     `$derived` also contains covariance matrices, which are simply the
##'     correlation matrices in `$corr` scaled by the standard deviations
##'     in `$sigma`. Relationship-level (`"dSigma"`) and case-level (`"pSigma"`)
##'     covariance matrices are always included, as is the group-level matrix
##'     (`"gSigma"`) when estimated.
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
##' @param srm.param Optional `character` vector (or `NULL`) specifying the type
##'   of SRM parameter to provide for `component`.
##'   Can be `%in% c("sd", "cor", "cov")`, as well as `"mean"`
##'   when `fixed.groups=FALSE`. The `summary()` method will always return the
##'   means as a group-level statistic, even for a single round-robin group.
##' @param cor.full `logical` indicating whether the full dyad-level correlation
##'   matrix is returned.  Ignored when `srm.param` does not include `"cor"`.
##'   The default (FALSE) returns only nonredundant dyad-level correlations, in
##'   a matrix with 1 row/column per round-robin variables:
##'
##'   - dyadic reciprocity on the diagonal
##'   - intra-case (e.g., within-person) correlations below the diagonal
##'   - inter-case (e.g., between-person) correlations above the diagonal
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
##'   `NULL` (the default) to avoid returning interval estimates.
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
##'   parameters at the requested level of analysis, multiplied by that level's
##'   sample size.  This is only useful for second-stage SEM estimation, which
##'   requires passing this to the `NACOV=` argument in [lavaan::lavaan()].
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
         slots = c(version  = "vector",  # stan, rstan, and lavaan.srm versions
                   call     = "call",    # mvsrm() call
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
#TODO: for cov/cor matrices, return single matrix for interval estimates
#      (lower / upper limits can go below / above the diagonal)
summary_mvSRM <- function(object, component = c("case","dyad"),
                          srm.param = c("sd","cor"), cor.full = FALSE,
                          posterior.est = "mean",
                          interval = NULL, credMass = .95,
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
  if ("cor" %in% srm.param) {
    stopifnot(is.logical(cor.full))
    cor.full <- isTRUE(cor.full[1])
  }

  ## can be NULL
  posterior.est <- intersect(posterior.est,
                             c("mean", "median", "mode", "EAP", "MAP"))
  interval <- intersect(tolower(interval), c("central", "hdi"))

  output <- list()

  ## add mean structure to output?  Always in "group" component
  if ("mean" %in% srm.param) {
    if (length(posterior.est))  for (p in posterior.est) {
      output$group$mean[[p]] <- as_matrix_mvSRM(x = object, srm.param = "mean",
                                                posterior.est = p, ...)
    }
    if (length(interval)) {
      if ("central" %in% interval) {
        probs <- abs(0:1 - (1 - credMass)/2)
        lower <- as_matrix_mvSRM(object, srm.param = "mean",
                                 posterior.est = probs[1], ...)
        upper <- as_matrix_mvSRM(object, srm.param = "mean",
                                 posterior.est = probs[2], ...)
        output$group$mean$central <- rbind(lower = lower, upper = upper)
      }
      if ("hdi" %in% interval) {
        probs <- abs(0:1 - (1 - credMass)/2)
        output$group$mean$hdi <- as_matrix_mvSRM(object, srm.param = "mean",
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
          output[[COMP]]$sd[[p]] <- as_matrix_mvSRM(x = object,
                                                    srm.param = "sd",
                                                    component = COMP,
                                                    posterior.est = p, ...)
        }
        ## add intervals?
        if (length(interval)) {
          if ("central" %in% interval) {
            probs <- abs(0:1 - (1 - credMass)/2)
            lower <- as_matrix_mvSRM(object, srm.param = "sd", component = COMP,
                                     posterior.est = probs[1], ...)
            upper <- as_matrix_mvSRM(object, srm.param = "sd", component = COMP,
                                     posterior.est = probs[2], ...)
            output[[COMP]]$sd$central <- rbind(lower = lower, upper = upper)
          }
          if ("hdi" %in% interval) {
            output$group$sd$hdi <- as_matrix_mvSRM(object, srm.param = "sd",
                                                   component = COMP,
                                                   posterior.est = "hdi",
                                                   credMass = credMass, ...)
          }
        }

      }

      ## then correlations
      if ("cor" %in% srm.param) {

        if (length(posterior.est))  for (p in posterior.est) {
          output[[COMP]]$cor[[p]] <- as_matrix_mvSRM(x = object,
                                                     srm.param = "cor",
                                                     cor.full = cor.full,
                                                     component = COMP,
                                                     posterior.est = p, ...)
        }
        ## add intervals?
        if (length(interval)) {
          if ("central" %in% interval) {
            probs <- abs(0:1 - (1 - credMass)/2)
            output[[COMP]]$cor$central$lower <- as_matrix_mvSRM(x = object,
                                                                srm.param = "cor",
                                                                cor.full = cor.full,
                                                                component = COMP,
                                                                posterior.est = probs[1],
                                                                ...)
            output[[COMP]]$cor$central$upper <- as_matrix_mvSRM(x = object,
                                                                srm.param = "cor",
                                                                cor.full = cor.full,
                                                                component = COMP,
                                                                posterior.est = probs[2],
                                                                ...)
          }
          if ("hdi" %in% interval) {
            output$group$cor$hdi <- as_matrix_mvSRM(object, srm.param = "cor",
                                                    cor.full = cor.full,
                                                    component = COMP,
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
          output[[COMP]]$cov[[p]] <- as_matrix_mvSRM(x = object,
                                                     srm.param = "cov",
                                                     component = COMP,
                                                     posterior.est = p, ...)
        }
        ## add intervals?
        if (length(interval)) {
          if ("central" %in% interval) {
            probs <- abs(0:1 - (1 - credMass)/2)
            output[[COMP]]$cov$central$lower <- as_matrix_mvSRM(x = object,
                                                                srm.param = "cov",
                                                                component = COMP,
                                                                posterior.est = probs[1],
                                                                ...)
            output[[COMP]]$cov$central$upper <- as_matrix_mvSRM(x = object,
                                                                srm.param = "cov",
                                                                component = COMP,
                                                                posterior.est = probs[2],
                                                                ...)
          }
          if ("hdi" %in% interval) {
            output$group$cov$hdi <- as_matrix_mvSRM(object, srm.param = "cov",
                                                    component = COMP,
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
setMethod("summary", "mvSRM", summary_mvSRM)


#TODO: add srm.param = "dyadic.recip" (vector), "intra.inter.cor" (2-tri matrix)
#TODO: optionally add.header=FALSE for summary() to set TRUE (its default)
##' @importFrom methods as
##' @importFrom stats quantile
##' @importFrom utils combn
as_matrix_mvSRM <- function(x, component, srm.param, posterior.est = "mean",
                            as.stanfit = FALSE, cor.full = FALSE, ...) {
  #TODO: Test this EXPERIMENTAL idea, implement with public argument?
  if (is.list(x)) {
    ## Assume all are valid object= arguments with the same structure.
    stopifnot(all(sapply(x, inherits, what = "mvSRM")))
    ## Assign the first element to object
    objectList <- x
    x <- objectList[[1]]
  } else objectList <- NULL


  if (as.stanfit) {
    if (!is.null(objectList)) {
      message("Argument x= not provided. Running as.matrix() method for the ",
              "first stanfit-class element in objectList=")
    }
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
    if (is.null(objectList)) {
      paramMat <- do.call(rbind, rstan::As.mcmc.list(x, pars = "Mvec"))
    } else {
      paramMat <- do.call(rbind, sapply(objectList, rstan::As.mcmc.list,
                                         pars = "Mvec"))
    }
    ## obtain summary statistic
    if (is.numeric(posterior.est)) {
      out <- apply(paramMat, MARGIN = 2, stats::quantile,
                   probs = posterior.est, ...)
    } else out <- apply(paramMat, MARGIN = 2, posterior.est, ...)

    #TODO: c() additional names of case/group_data or dyad-constant variables
    if (posterior.est == "hdi") {
      colnames(out) <- names(x@varNames$RR)
      class(out) <- c("lavaan.matrix", "matrix")
    } else {
      names(out) <- names(x@varNames$RR)
      class(out) <- c("lavaan.vector", "numeric")
    }
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
    if (srm.param == "sd") {
      PARS <- "S_g"
    } else if (srm.param == "cor") {
      PARS <- "Rg"
    } else if (srm.param == "cov") {
      PARS <- "gSigma"
    }
    NAMES <- x@varNames$group

  } else if (component == "case") {
    if (srm.param == "sd") {
      PARS <- "S_p"
    } else if (srm.param == "cor") {
      PARS <- "Rp"
    } else if (srm.param == "cov") {
      PARS <- "pSigma"
    }
    NAMES <- x@varNames$case

  } else if (component == "dyad") {
    if (srm.param == "sd") {
      PARS <- "s_rr"
      NAMES <- c(names(x@varNames$RR), x@varNames$dyad)
    } else if (srm.param == "cor") {
      PARS <- ifelse(cor.full, "Rd2", "r_d2")
      NAMES <- c(names(x@varNames$RR), x@varNames$dyad)
    } else if (srm.param == "cov") {
      PARS <- "dSigma"
      ## extended names for covariance matrix
      NAMES <- c(c(x@varNames$RR, recursive = TRUE, use.names = FALSE),
                 x@varNames$dyad)
    }

  }


  ## extract samples
  if (is.null(objectList)) {
    paramMat <- do.call(rbind, rstan::As.mcmc.list(x, pars = PARS))
  } else {
    paramMat <- do.call(rbind, sapply(objectList, rstan::As.mcmc.list,
                                      pars = PARS))
  }


  ## obtain summary statistic
  if (is.numeric(posterior.est)) {
    out <- apply(paramMat, MARGIN = 2, stats::quantile,
                 probs = posterior.est, ...)
  } else out <- apply(paramMat, MARGIN = 2, posterior.est, ...)

  if (srm.param == "sd") {
    ## attach variable names
    if (posterior.est == "hdi") {
      colnames(out) <- NAMES #TODO: c() case/group_data names
      class(out) <- c("lavaan.matrix", "matrix")
    } else {
      names(out) <- NAMES
      class(out) <- c("lavaan.vector", "numeric")
    }

  } else if (posterior.est == "hdi") {
    lower <- out[1,]
    upper <- out[2,]
    ## store correlations/covariances in a list of 2 matrices ($lower & $upper)
    assign(PARS, matrix(0, nrow = length(NAMES), ncol = length(NAMES),
                        dimnames = list(NAMES, NAMES)))
    for (i in names(lower)) {
      eval(parse(text = paste(i, "<-", lower[i]) ))
    }
    LOWER <- eval(as.name(PARS))
    class(LOWER) <- c(ifelse(PARS == "r_d2", "lavaan.matrix",
                             "lavaan.matrix.symmetric"), "matrix")
    for (i in names(upper)) {
      eval(parse(text = paste(i, "<-", upper[i]) ))
    }
    UPPER <- eval(as.name(PARS))
    class(UPPER) <- c(ifelse(PARS == "r_d2", "lavaan.matrix",
                             "lavaan.matrix.symmetric"), "matrix")
    ## transform raw beta parameters to correlation metric?
    if (PARS == "r_d2") {
      LOWER <- LOWER*2 - 1
      UPPER <- UPPER*2 - 1
    }
    ## nothing else to do for correlation/covariance limits
    out <- list(lower = LOWER, upper = UPPER)

  } else {
    ## store correlations/covariances in a matrix
    assign(PARS, matrix(0, nrow = length(NAMES), ncol = length(NAMES),
                        dimnames = list(NAMES, NAMES)))
    for (i in names(out)) {
      eval(parse(text = paste(i, "<-", out[i]) ))
    }
    out <- eval(as.name(PARS))
    class(out) <- c(ifelse(PARS == "r_d2", "lavaan.matrix",
                           "lavaan.matrix.symmetric"), "matrix")
    ## transform raw beta parameters to correlation metric?
    if (PARS == "r_d2") out <- out*2 - 1
  }

  out
}
##' @name mvSRM-class
##' @aliases as.matrix,mvSRM-method
##' @export
setMethod("as.matrix", "mvSRM", as_matrix_mvSRM)



##' @importFrom stats cov
##' @importFrom stats setNames
vcov_mvSRM <- function(object, component, meanstructure = FALSE,
                       keep, drop, ...) {
  #TODO: Test this EXPERIMENTAL idea, implement with public argument?
  if (is.list(object)) {
    ## Assume all are valid object= arguments with the same structure.
    stopifnot(all(sapply(object, inherits, what = "mvSRM")))
    ## Assign the first element to object
    objectList <- object
    object <- objectList[[1]]
  } else objectList <- NULL

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

    PARS <- "gSigma"
    NAMES <- object@varNames$group

    if (!missing(keep)) {
      SUBSET <- intersect(keep, NAMES)
      if (!length(SUBSET)) stop('keep=', keep, ' leaves no variables to return')
    } else if (!missing(drop)) {
      SUBSET <- setdiff(NAMES, drop)
      if (!length(SUBSET)) stop('drop=', drop, ' leaves no variables to return')
    } else SUBSET <- NAMES

    if (meanstructure) {
      if (is.null(objectList)) {
        MuSamples <- do.call(rbind, rstan::As.mcmc.list(object, pars = "Mvec"))
      } else {
        MuSamples <- do.call(rbind, sapply(objectList, rstan::As.mcmc.list,
                                           pars = "Mvec"))
      }

      MuList <- apply(MuSamples, MARGIN = 1, FUN = function(m) {
        Mvec <- setNames(numeric(length(NAMES)), nm = NAMES)
        for (i in names(m)) eval(parse(text = paste(i, "<-", m[i]) ))
        Mvec[SUBSET]
      }, simplify = FALSE)

    }


  } else if (component == "case") {

    PARS <- "pSigma"
    NAMES <- object@varNames$case

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

    PARS <- "dSigma"
    NAMES <- c(c(object@varNames$RR, recursive = TRUE, use.names = FALSE),
               object@varNames$dyad)

    if (!missing(keep)) {
      KEEP <- do.call(c, sapply(keep, function(k) {
        if (k %in% names(object@varNames$RR)) {
          return(paste0(k, c("_ij", "_ji")))
        } else if (k %in% NAMES) {
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
        } else if (d %in% NAMES) {
          return(d)
        }
        NULL
      }, simplify = FALSE))
      SUBSET <- setdiff(NAMES, DROP)
      if (!length(SUBSET)) stop('drop=', drop, ' leaves no variables to return')

    } else SUBSET <- NAMES
  }


  ## store covariances in a matrix (per posterior sample)
  if (is.null(objectList)) {
    covMats <- do.call(rbind, rstan::As.mcmc.list(object, pars = PARS))
  } else {
    covMats <- do.call(rbind, sapply(objectList, rstan::As.mcmc.list,
                                     pars = PARS))
  }
  SigmaList <- apply(covMats, MARGIN = 1, FUN = function(m) {
    assign(PARS, matrix(0, nrow = length(NAMES), ncol = length(NAMES),
                        dimnames = list(NAMES, NAMES)))
    for (i in names(m)) eval(parse(text = paste(i, "<-", m[i]) ))
    return( eval(as.name(PARS))[SUBSET, SUBSET] )
  }, simplify = FALSE)

  if (!exists("MuList")) MuList <- list(NULL)
  ## posterior variability
  postList <- mapply(function(S, M = NULL) {
    namesMatrix <- outer(X = rownames(S), Y = colnames(S), paste, sep = "~~")
    if (categorical) {
      wls.obs <- c(M, # 1. interleaved thresholds + (negative) means, if any
                   # 2. slopes (if any)
                   # 3. variances (of continuous variables, if any)
                   as.numeric(diag(S)),
                   # 4. lower.tri covariances/correlations
                   lavaan::lav_matrix_vech(S, diagonal = FALSE))
      DIMNAMES <- c(diag(namesMatrix),
                    lavaan::lav_matrix_vech(namesMatrix, diagonal = FALSE))
    } else {
       wls.obs <- c(M, # 1. means (if any)
                    # 2. slopes (if any)
                    # 3. lower.tri (co)variance matrix
                    lavaan::lav_matrix_vech(S, diagonal = TRUE))
       DIMNAMES <- lavaan::lav_matrix_vech(namesMatrix, diagonal = TRUE)
    }
    if (!is.null(M)) DIMNAMES <- c(paste0(names(M), "~1"), DIMNAMES)
    names(wls.obs) <- DIMNAMES
    wls.obs
  }, S = SigmaList, M = MuList, SIMPLIFY = FALSE)
  postSamp <- do.call(rbind, postList)
  NACOV <- object@nobs[[component]] * cov(postSamp)
  class(NACOV) <- c("lavaan.matrix.symmetric", "matrix", "array")

  ## add an attribute used by srm2lavData()?
  if (isTRUE(list(...)$add.names.attr)) attr(NACOV, "subset") <- SUBSET
  NACOV
}
##' @name mvSRM-class
##' @aliases vcov,mvSRM-method
##' @importFrom stats vcov
##' @export
setMethod("vcov", "mvSRM", vcov_mvSRM)



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


