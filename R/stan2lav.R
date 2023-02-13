### Terrence D. Jorgensen
### Last updated: 13 February 2023
### (currently hidden) function to create a lavMoments-class object
### from a mvSRM-class object (inherits from stanfit-class)

## @param object An object of class `mvSRM`
## @param component `character` specifying which SRM component to provide
##   estimated parameters for. Can be `%in% c("case", "dyad")`, as well as
##   `"group"` when group effects are modeled.
## @param point `character` indicating choice of posterior point estimate.
##   Must be `%in% c("mean", "median", "mode")`.
## @param meanstructure `logical` indicating whether to include estimated means
##    (when available), only for `component="group"`.
## @param lavData An object of class `lavMoments`. When `NULL`, a new instance
##   of `lavMoments` will be created for the specified `component=`. If
##   provided, the new results (for currently specified `component=`) will be
##   concatenated with the existing results in `lavData=`.  No checks are
##   performed, so take care that it makes sense.  This can enable fitting
##   simultaneous models to different (types of) round-robin groups or to
##   different levels (e.g., concatenating case- and dyad-level results). It is
##   up to the user to create sensible [lavaan::model.syntax()] for models
##   that include multiple components; the variable names will differ between
##   case- and dyad-level components, so syntax blocks will be required, rather
##   than standard multigroup syntax (see [lavaan::Demo.twolevel()] for an
##   example of block syntax, but replace `level:` with `group:`).
## @param keep,drop Variables to include (`keep=`) or exclude (`drop=`) from
##   the result.  `drop=` is ignored when `keep=` is specified.  For case- or
##   dyad-level variables, `keep`/`drop` can be any combination of the original
##   names of round-robin variables passed to [mvsrm()] or the names of
##   specific components, which include `c("out","in")` suffixes for case-level
##   components or `c("ij","ji")` suffixes for dyad/relationship-level
##   components.
#TODO (if this becomes public): create a syntax example, verify blocks work
## @param ... Any [lavaan::lavOptions()] to be passed to [lavaan::lavaan()]
srm2lavData <- function(object, component, point = "mean", keep, drop,
                        meanstructure = FALSE, lavData = NULL, ...) {
  stopifnot(inherits(object, "mvSRM"))
  categorical <- FALSE #TODO: specify threshold model in Stan
  component <- tolower(component[1]) # one at a time



  COV <- summary.mvSRM(object, component = component, stat = "cov",
                       point = point, interval = NULL)[[component]]$cov[[point]]
  COV <- setNames(list(COV), nm = component)

  if (meanstructure && component == "group") {
    M <- summary.mvSRM(object, stat = "mean", point = point,
                       interval = NULL)$group$mean[[point]]
    M <- setNames(list(M), nm = component)
  } else M <- NULL

  NACOV <- vcov.mvSRM(object, component = component, keep = keep, drop = drop,
                      meanstructure = meanstructure)
  if (categorical) {
    #TODO: thresholds

    #FIXME?  Always needed, even for MLE, continuous?  ask Yves
    WLS.V <- diag(diag(solve(NACOV)))
    attributes(WLS.V) <- attributes(NACOV)
    WLS.V <- setNames(list(WLS.V), nm = component)
  }
  NACOV <- setNames(list(NACOV), nm = component)

  ## apply keep/drop to sample.stats
  if (!is.null(attr(NACOV, "subset"))) {
    SUBSET <- attr(NACOV, "subset")
    attr(NACOV, "subset") <- NULL

    COV <- COV[SUBSET, SUBSET]
    if (!is.null(M)) M <- M[SUBSET]
    #TODO: thresholds
  }

  if (is.null(lavData)) {
    ## make a new lavMoments object
    out <- list(sample.cov  = COV,
                sample.mean = M,
                sample.nobs = object@N[component],
                NACOV       = NACOV)
    if (categorical) {
      #TODO: thresholds
      out$WLS.V <- WLS.V
    }


  } else {
    ## concatenate after existing lavData
    out <- lavData

    out$sample.cov  <- c(out$sample.cov ,  COV )
    out$sample.mean <- c(out$sample.mean,   M  )
    out$sample.nobs <- c(out$sample.nobs, object@N[component])
    out$NACOV       <- c(out$NACOV      , NACOV)

    if (categorical) {
      #TODO: thresholds
      out$WLS.V <- c(out$WLS.V, WLS.V)
    }

  }


  ## set recommended arguments (estimator, se, test)
  out$lavOptions <- list(sample.cov.rescale = FALSE,
                         fixed.x            = FALSE, #TODO: model-specific NACOV based on ov.x
                         conditional.x      = FALSE,
                         estimator          = ifelse(categorical, "DWLS","ML"),
                         se                 = "robust.sem",
                         test               = "Browne.residual.adf")
  ## set class and return
  class(out) <- c("lavMoments", "list")
  out
}



## copied from semTools
#FIXME: add semTools dependency?  Already in Suggests for plausibleValues()?
## @exportS3Method print lavMoments
# print.lavMoments <- function(x, ...) {
#   nameX <- substitute(x)
#
#   cat('This lavMoments-class object contains summary statistics in a list ',
#       'consisting of the following elements:\n  ',
#       sep = '')
#   cat(names(x), sep = ', ')
#   cat('\nYou can view list elements with str(', nameX, ') or standard ',
#       'extraction methods (e.g., ', nameX, '$sample.cov).\n\n', sep = '')
#
#   cat('This object can be passed to lavaan() as the data= argument, in which ',
#       'case each element in this object (e.g., $sample.cov) will be internally ',
#       'passed to the corresponding lavaan() argument.\n', sep = '')
#
#   if ("lavOptions" %in% names(x)) {
#     cat('The following recommended lavOptions() will also be passed to lavaan():\n  ',
#         sep = '')
#     cat(paste0(names(x$lavOptions),
#                ifelse(sapply(x$lavOptions, is.character), ' = "', ' = '),
#                x$lavOptions,
#                ifelse(sapply(x$lavOptions, is.character), '"', '')),
#         sep = ',\n  ')
#     cat('\nYou can override these settings in your lavaan() call, but you ',
#         'should NOT set fixed.x=TRUE, which would yield incorrect standard ',
#         'errors and test statistics.\n\n',
#         sep = '')
#   }
#
#   return(invisible(x))
# }



