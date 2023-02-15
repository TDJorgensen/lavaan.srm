### Terrence D. Jorgensen
### Last updated: 15 February 2023
### pass model, lavMoments, and other arguments to lavaan()


##' Fit a `lavaan` Model to Multivariate SRM Results
##'
##' Estimate a structural equation model (SEM) parameters, using MCMC-estimated
##' summary statistics from a multivariate social relations model (mvSRM;
##' Nestler, [2018](https://doi.org/10.3102/1076998617741106)) as data.
##'
##'
##' @param model `character` specifying an SEM.  See [lavaan::model.syntax()].
##' @param data An object of class \code{\linkS4class{mvSRM}}
##              (or lavMoments, not documented)
## @param component `character` specifying which SRM component(s) the `model=`
##   syntax is specified for.  Can be `%in% c("case", "dyad")`, as well as
##   `"group"` when group effects are modeled.  Ignored when `stat="mean"`
##   because means are only a group-level statistic.
##' @param point `character` indicating choice of posterior point estimate.
##'   Must be `%in% c("mean", "median", "mode")`.
##' @param ... Arguments and [lavaan::lavOptions()] passed to [lavaan::lavaan()]
##'
##' @return An object of `class?lavaan`.
##'
##'
##' @seealso \code{\linkS4class{mvSRM}}
##' @examples
##'
##'
##' \dontrun{
##'
##' ## STAGE 1: Estimate the multivariate SRM
##'
##'
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
##' ## STAGE 2: Specify and fit an SEM
##'
##'
##' ## specify PERSON-level model
##' modP <- ' ## factor loadings
##'   f1_out =~ Wert1_out + Wert2_out + Wert3_out
##'   f1_in  =~ Wert1_in  + Wert2_in  + Wert3_in
##'
##' ## correlated residuals
##'   Wert1_out ~~ Wert1_in
##'   Wert2_out ~~ Wert2_in
##'   Wert3_out ~~ Wert3_in
##' '
##'
##' ## fit CFA to PERSON-level components,
##' ## using mvSRM-class object (Stage 1) as input data
##' fitP <- lavaan.srm(modP, data = srmOut, # more arguments passed to lavaan()
##'                    std.lv = TRUE, auto.var = TRUE, auto.cov.lv.x = TRUE)
##' ## fitP is a lavaan-class object, so any lavaan function is available
##' summary(fitP, std = TRUE, rsq = TRUE) # only use Browne's residual-based test
##'
##'
##'
##' ## specify DYAD-level model
##' modD <- ' # equal loadings
##'   f1_ij =~ L1*Wert1_ij + L2*Wert2_ij + L3*Wert3_ij
##'   f1_ji =~ L1*Wert1_ji + L2*Wert2_ji + L3*Wert3_ji
##'
##' ## equal variances
##'   f1_ij ~~ phi*f1_ij
##'   f1_ji ~~ phi*f1_ji
##'
##'   Wert1_ij ~~ v1*Wert1_ij
##'   Wert2_ij ~~ v2*Wert2_ij
##'   Wert3_ij ~~ v3*Wert3_ij
##'
##'   Wert1_ji ~~ v1*Wert1_ji
##'   Wert2_ji ~~ v2*Wert2_ji
##'   Wert3_ji ~~ v3*Wert3_ji
##'
##' ## correlated residuals
##'   Wert1_ij ~~ Wert1_ji
##'   Wert2_ij ~~ Wert2_ji
##'   Wert3_ij ~~ Wert3_ji
##' '
##'
##' ## fit CFA to DYAD-level components,
##' ## using mvSRM-class object (Stage 1) as input data
##' fitD <- lavaan.srm(modD, data = srmOut, # more arguments passed to lavaan()
##'                    std.lv = TRUE, auto.var = TRUE, auto.cov.lv.x = TRUE)
##' summary(fitD, std = TRUE, rsq = TRUE)
##'
##'
#TODO: model both at once
##' }
##'
##'
##'
##' @importFrom lavaan lavaan lavParseModelString lavInspect
##' @export
lavaan.srm <- function(model, data, point = "mean", ...) {
  MC <- match.call(expand.dots = TRUE) # to overwrite lavaan's default @call
  dots <- list(...)

  ## lavOptions that can be passed to lavParseModelString()
  WARN  <- dots$WARN
  DEBUG <- dots$debug
  ## if NULL, restore defaults
  if (is.null(WARN))  WARN  <- formals(lavParseModelString)$warn
  if (is.null(DEBUG)) DEBUG <- formals(lavParseModelString)$debug

  ## make template parTable to extract variable names
  if (is.character(model)) {
    PT <- lavParseModelString(model, as.data.frame. = TRUE,
                              warn = WARN, debug = DEBUG)
  } else if (is.list(model)) {
    PT <- as.data.frame(model) # no further checks, it will break if wrong
  } else {
    stop('model= must be a character string or paramter table (data.frame or ',
         'list). See ?model.syntax and ?parTable for details.')
  }

  lavars <- c(PT$lhs, PT$rhs)


  if (inherits(data, "lavMoments")) {
    srmMoments <- data

    ## check for group-level variables
    ov.group <- intersect(lavars, rownames(data$sample.cov$group))
    extra.group <- setdiff(rownames(data$sample.cov$group), ov.group)
    if (length(extra.group)) stop('There are some observed variables in the ',
                                  'data@sample.cov matrix that do not appear ',
                                  'in the model= syntax, so it will not be ',
                                  'consistent with data@NACOV:\n\t',
                                  paste(extra.group, collapse = ", "))

    ## check for case-level variables
    caseVars <- c(grep(pattern = "_out", x = lavars, fixed = TRUE, value = TRUE),
                  grep(pattern = "_in" , x = lavars, fixed = TRUE, value = TRUE))
    ov.case <- intersect(caseVars, rownames(data$sample.cov$case))
    extra.case <- setdiff(rownames(data$sample.cov$case), ov.case)
    if (length(extra.case)) stop('There are some observed variables in the ',
                                 'data@sample.cov matrix that do not appear ',
                                 'in the model= syntax, so it will not be ',
                                 'consistent with data@NACOV:\n\t',
                                 paste(extra.case, collapse = ", "))

    ## check for dyad-level variables
    dyadVars <- c(grep(pattern = "_ij", x = lavars, fixed = TRUE, value = TRUE),
                  grep(pattern = "_ji", x = lavars, fixed = TRUE, value = TRUE))
    ov.dyad <- intersect(dyadVars, rownames(data$sample.cov$dyad))
    extra.dyad <- setdiff(rownames(data$sample.cov$dyad), ov.dyad)
    if (length(extra.dyad)) stop('There are some observed variables in the ',
                                 'data@sample.cov matrix that do not appear ',
                                 'in the model= syntax, so it will not be ',
                                 'consistent with data@NACOV:\n\t',
                                 paste(extra.dyad, collapse = ", "))


  } else if (inherits(data, "mvSRM")) {
    ## must create a lavMoments object

    ## check for group-level variables
    ov.group <- intersect(lavars, data@varNames$group)

    ## check for case-level variables
    caseVars <- c(grep(pattern = "_out", x = lavars, fixed = TRUE, value = TRUE),
                  grep(pattern = "_in" , x = lavars, fixed = TRUE, value = TRUE))
    ov.case <- intersect(caseVars, data@varNames$case)

    ## check for dyad-level variables
    dyadVars <- c(grep(pattern = "_ij", x = lavars, fixed = TRUE, value = TRUE),
                  grep(pattern = "_ji", x = lavars, fixed = TRUE, value = TRUE))
    rrVars <- do.call(c, c(data@varNames$RR,  list(use.names = FALSE)))
    ov.dyad <- intersect(dyadVars, rrVars)


    #TODO: call summary() and vcov()
    srmMoments <- NULL # initialize
    srm2lavArgs <- list(srm2lavData, object = data, point = point)

    if (length(ov.group)) {
      srm2lavArgs$component <- "group"
      srm2lavArgs$keep <- ov.group
      if (any(PT$op == "~1")) srm2lavArgs$meanstructure <- TRUE
      srmMoments <- eval(as.call(srm2lavArgs))
    }
    if (length(ov.case)) {
      srm2lavArgs$component     <- "case"
      srm2lavArgs$lavData       <- srmMoments
      srm2lavArgs$keep          <- ov.case
      srm2lavArgs$meanstructure <- FALSE
      srmMoments <- eval(as.call(srm2lavArgs))
      if (any(PT$op == "~1")) {
        ## sample.mean must be a vector of zeros
        srmMoments$sample.mean <- c(srmMoments$sample.mean,
                                    case = setNames(numeric(length(ov.case)),
                                             nm = ov.case))
      }
    }
    if (length(ov.dyad)) {
      srm2lavArgs$component <- "dyad"
      srm2lavArgs$lavData <- srmMoments
      srm2lavArgs$keep          <- ov.dyad
      srm2lavArgs$meanstructure <- FALSE
      srmMoments <- eval(as.call(srm2lavArgs))
      if (any(PT$op == "~1")) {
        ## sample.mean must be a vector of zeros
        srmMoments$sample.mean <- c(srmMoments$sample.mean,
                                    case = setNames(numeric(length(ov.dyad)),
                                                    nm = ov.dyad))
      }
    }
    if (!length(c(ov.group, ov.case, ov.dyad))) {
      stop('No variables in syntax match variable names in data= argument.')
    }

  } else stop('data= must be an object of class?mvSRM (or a "lavMoments" list)')


  ## fit target model
  fit <- lavaan(model, data = srmMoments, ...)

  ## fit additional models?
  if (inherits(data, "lavMoments")) {
    warning('Cannot automatically specify h1 or baseline model unless ',
            'data= is a class?mvSRM object.')
    #FIXME: fail-safe way to check for pairs of RR names?
    #       Have to assume "_ij" only appears at the end.

  } else if (lavInspect(fit, "converged")) {
    ## check whether defaults are overwritten
    if (is.null(dots$h1)) dots$h1 <- lavaan::lavOptions("h1")$h1
    if (is.null(dots$baseline))
      dots$baseline <- lavaan::lavOptions("baseline")$baseline

    ## need custom dyad-level h1 or baseline model?
    RRnm <- names(data@varNames$RR)
    if ((dots$h1 || dots$baseline) && length(ov.dyad)) {

      dyad.baseline <- character(0) # initialize
      dyad.h1       <- character(0)

      ## start with all RR variables
      for (v1 in seq_along(RRnm)) {
        if (v1 %in% ov.dyad) next

        nm1_ij <- paste(RRnm[v1], "ij", sep = "_")
        nm1_ji <- paste(RRnm[v1], "ji", sep = "_")

        for (v2 in v1:length(RRnm)) {
          if (!v2 %in% ov.dyad) next

          nm2_ij <- paste(RRnm[v2], "ij", sep = "_")
          nm2_ji <- paste(RRnm[v2], "ji", sep = "_")

          if (v1 == v2) {
            ## equate variances
            if (nm1_ij %in% ov.dyad) {
              tmp <- paste0(nm1_ij, " ~~ var", v1, "*", nm1_ij)
              dyad.baseline <- c(dyad.baseline, tmp)
            }
            if (nm1_ji %in% ov.dyad) {
              tmp <- paste0(nm1_ji, " ~~ var", v1, "*", nm1_ji)
              dyad.baseline <- c(dyad.baseline, tmp)
            }

            ## dyadic reciprocity
            if (nm1_ij %in% ov.dyad && nm1_ji %in% ov.dyad) {
              tmp <- paste0(nm1_ij, " ~~ ", nm1_ji)
              dyad.baseline <- c(dyad.baseline, tmp)
            }

          } else {
            if (!dots$h1) next

            ## inter-cor
            if (nm1_ij %in% ov.dyad && nm2_ji %in% ov.dyad) {
              tmp <- paste0(nm1_ij, " ~~ inter", v1, v2, "*", nm2_ji)
              dyad.h1 <- c(dyad.h1, tmp)
            }
            if (nm2_ij %in% ov.dyad && nm1_ji %in% ov.dyad) {
              tmp <- paste0(nm2_ij, " ~~ inter", v1, v2, "*", nm1_ji)
              dyad.h1 <- c(dyad.h1, tmp)
            }

            ## intra-cor
            if (nm1_ij %in% ov.dyad && nm2_ij %in% ov.dyad) {
              tmp <- paste0(nm1_ij, " ~~ intra", v1, v2, "*", nm2_ij)
              dyad.h1 <- c(dyad.h1, tmp)
            }
            if (nm1_ji %in% ov.dyad && nm2_ji %in% ov.dyad) {
              tmp <- paste0(nm1_ji, " ~~ intra", v1, v2, "*", nm2_ji)
              dyad.h1 <- c(dyad.h1, tmp)
            }

          }

        } # end v2
      }   # end v1

      ## concatenate other levels?
      #TODO: add "group: 1" to top of dyad syntax before concatenating
      if (length(ov.case)) {
        sat.case <- outer(ov.case, ov.case, paste, sep = "~~")
        sat.case[lower.tri(sat.case, diag = TRUE)]
      }
      if (length(ov.group)) {
        sat.group <- outer(ov.group, ov.group, paste, sep = "~~")
        sat.group[lower.tri(sat.group, diag = TRUE)]
      }

    }

    ## need custom dyad-level baseline model?
    if (dots$baseline & length(ov.dyad)) {

    }

  }


  fit
}
