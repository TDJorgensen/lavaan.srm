### Terrence D. Jorgensen
### Last updated: 19 February 2023
### pass model, lavMoments, and other arguments to lavaan()


##' Fit a `lavaan` Model to Multivariate SRM Results
##'
##' Estimate a structural equation model (SEM) parameters, using MCMC-estimated
##' summary statistics from a multivariate social relations model (mvSRM;
##' Nestler, [2018](https://doi.org/10.3102/1076998617741106)) as data.
##'
##'
#TODO: DETAILS above how syntax works for multiple groups/components
##'
##'
##' @param model `character` specifying an SEM.  See [lavaan::model.syntax()].
##' @param data An object of class \code{\linkS4class{mvSRM}}
##              (or lavMoments, not documented)
##' @param component `character` specifying which SRM component(s) the `model=`
##'   syntax is specified for.  Can be `%in% c("case", "dyad")`, as well as
##'   `"group"` when group effects are modeled.  The order of multiple blocks
##'   in the `model=` syntax must match the order specified in `component=`.
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
##' fitP <- lavaan.srm(modP, data = srmOut, component = "case",
##'                    ## more arguments passed to lavaan()
##'                    std.lv = TRUE, auto.var = TRUE, auto.cov.lv.x = TRUE)
##' ## fitP is a lavaan-class object, so any lavaan function is available
##' summary(fitP, std = TRUE, rsq = TRUE, fit.measures = TRUE,
##'         ## Ignore standard test, only use Browne's residual-based test.
##'         ## Also use Browne's residual-based test to calculate fit indices:
##'         fm.args = list(standard.test = "browne.residual.adf"))
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
##' ## using wrapper analogous to lavaan::cfa()
##' fitD <- cfa.srm(modD, data = srmOut, component = "dyad", std.lv = TRUE)
##' summary(fitD, std = TRUE, rsq = TRUE, fit.measures = TRUE,
##'         ## Ignore standard test, only use Browne's residual-based test.
##'         ## Also use Browne's residual-based test to calculate fit indices:
##'         fm.args = list(standard.test = "browne.residual.adf"))
##'
##'
##'
##' ## Simultaneously specifying BOTH person & dyad-level models
##' ## requires block-structured syntax, analogous to multilevel SEM.
##'
##' ## specify a model with cross-level measurement equivalence
##' modPD <- ' group: 1   # why does "group: case" not work?
##'
##'   f1_out =~ L1*Wert1_out + L2*Wert2_out + L3*Wert3_out
##'   f1_in  =~ L1*Wert1_in  + L2*Wert2_in  + L3*Wert3_in
##'
##' ## correlated residuals
##'   Wert1_out ~~ Wert1_in
##'   Wert2_out ~~ Wert2_in
##'   Wert3_out ~~ Wert3_in
##'
##' ## free factor variances (use dyad-level = 1 as reference "group")
##'   f1_out ~~ var_perceiver*f1_out
##'   f1_in  ~~ var_target*f1_in
##'
##'
##' group: 2   # why does "group: dyad" not work?
##'
##'   f1_ij =~ L1*Wert1_ij + L2*Wert2_ij + L3*Wert3_ij
##'   f1_ji =~ L1*Wert1_ji + L2*Wert2_ji + L3*Wert3_ji
##'
##' ## equal variances
##'   f1_ij ~~ var_relationship*f1_ij
##'   f1_ji ~~ var_relationship*f1_ji
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
##'
##'
##' ## Identification constraint: Sum of factor variances == 1
##'   var_relationship == 1 - (var_perceiver + var_target)
##' ## Thus, each estimated variance == proportion of total (ICC)
##' '
##' ## order in componenent= must match order in model syntax above
##' fitPD <- lavaan.srm(modPD, data = srmOut, component = c("case","dyad"),
##'                     auto.var = TRUE, auto.cov.lv.x = TRUE)
##'
##' summary(fitPD, std = TRUE, rsq = TRUEfit.measures = TRUE,
##'         ## Ignore standard test, only use Browne's residual-based test.
##'         ## Also use Browne's residual-based test to calculate fit indices:
##'         fm.args = list(standard.test = "browne.residual.adf"))
##'
##'
##' }
##'
##'
##'
##' @importFrom lavaan lavaan lavParseModelString lavInspect lavOptions
##' @export
lavaan.srm <- function(model, data, component, point = "mean", ...) {
  stopifnot(all(component %in% c("group","case","dyad")),
            point %in% c("mean","median","mode"))

  MC <- match.call(expand.dots = TRUE) # to overwrite lavaan's default @call
  dots <- list(...)

  ## lavOptions that can be passed to lavParseModelString()
  WARN  <- dots$WARN
  DEBUG <- dots$debug
  ## if NULL, restore defaults
  if (is.null(WARN))  WARN  <- formals(lavParseModelString)$warn
  if (is.null(DEBUG)) DEBUG <- formals(lavParseModelString)$debug
  ## don't overwrite user's setting, but ...
  if (is.null(dots$h1))       dots$h1       <- lavOptions("h1")$h1
  if (is.null(dots$baseline)) dots$baseline <- lavOptions("baseline")$baseline
  fitMore <- any(c("case","dyad") %in% component) && (dots$h1 || dots$baseline)

  ## make template parTable to extract variable names
  if (is.character(model)) {
    PT <- lavParseModelString(model, as.data.frame. = TRUE,
                              warn = WARN, debug = DEBUG)
  } else if (is.list(model)) {
    PT <- as.data.frame(model)
    if (is.null(PT$block)) {
      if (is.null(PT$group)) {
        PT$block <- 1L
        PT$block[PT$op %in% c(":=","==",">","<")] <- 0L
      } else PT$block <- PT$group
    }
    ## no further checks on data.frame, it will break if wrong

  } else {
    stop('model= must be a character string or paramter table (data.frame or ',
         'list). See ?model.syntax and ?parTable for details.')
  }

  blocks <- unique(PT$block)
  ## if multigroup, require same components analyzed in each group
  stopifnot(length(blocks) %% length(component) == 0)


  if (inherits(data, "lavMoments")) {
    srmMoments <- data

    #FIXME: assumes no multi(RR)group models fit to same component
    lavars <- c(PT$lhs, PT$rhs) # no need to distinguish blocks?

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

    if (length(ov.case) || length(ov.dyad))
      warning('Cannot automatically specify h1 or baseline case- or dyad-level ',
              'model unless data= is a class?mvSRM object.')
    #FIXME: fail-safe way to check for pairs of RR names?
    #       Have to assume "_ij" only appears at the end.

  } else if (inherits(data, "list")) {
    stopifnot(all(sapply(data, inherits, what = "mvSRM")))

    ## if multiple components, require same components analyzed in each group
    stopifnot(length(blocks) %% length(data) == 0)

    #TODO: Do same shit below per RR-group type
    #      Require components nested in groups?
    # rep(component, times = length(data))


  } else if (inherits(data, "mvSRM")) {
    ## to define for model and reuse for h1 / baseline
    ov.names <- list()
    ## initialize h1 & baseline syntax
    sat.mod <- character(0)
    bas.mod <- character(0)

    ## loop over blocks to synchronize order of components with syntax
    for (b in seq_along(blocks)) {
      ## extract block-b rows of parTable
      PTb <- PT[PT$block == blocks[b], ]
      ## save all variable names ("manifest" and latent)
      stxVars <- c(PTb$lhs, PTb$rhs)
      ## determine which are "manifest" for the specified component
      if (component[b] == "dyad") {
        srmVars <- c(do.call(c, c(data@varNames$RR,  list(use.names = FALSE))),
                     ## include symmetric variables (constant within dyad)
                     data@varNames$dyad)
        ov.names[[b]] <- intersect(stxVars, srmVars)
      } else {
        ov.names[[b]] <- intersect(stxVars, data@varNames[[ component[b] ]])
      }


      ## prepare arguments to create a lavMoments object
      srm2lavArgs <- list(srm2lavData, object = data, component = component[b],
                          point = point, keep = ov.names[[b]])
      if (b > 1L) srm2lavArgs$lavData <- srmMoments

      ## determine whether to include mean structure
      if (component[b] != "group") {
        srm2lavArgs$meanstructure <- FALSE
      } else if (any(PTb$op == "~1")) {
        srm2lavArgs$meanstructure <- TRUE
      } else if (!is.null(dots$meanstructure)) {
        srm2lavArgs$meanstructure <- dots$meanstructure
      } else srm2lavArgs$meanstructure <- FALSE         #TODO: document default

      ## generate lavMoments
      srmMoments <- eval(as.call(srm2lavArgs))

      ## need to add an ad-hoc mean structure?
      if (!srm2lavArgs$meanstructure && !is.null(srmMoments$sample.mean)) {
        ## sample.mean must be a vector of zeros for this block
        srmMoments$sample.mean[[component]] <- setNames(rep(0, length(ov.names[[b]])),
                                                        nm = ov.names[[b]])
      }


      ## fit SATURATED and BASELINE models
      if (length(blocks) > 1L) {
        ## add block structure
        sat.mod[[length(sat.mod) + 1L]] <- paste("group:", b)
        bas.mod[[length(bas.mod) + 1L]] <- paste("group:", b)
      }
      ## group-level model has no special constraints
      if (fitMore && component[b] == "group") {
        sat.group <- outer(ov.names[[b]], ov.names[[b]], paste, sep = "~~")
        sat.mod <- c(sat.mod, sat.group[lower.tri(sat.group, diag = TRUE)])
        bas.mod <- c(bas.mod, diag(sat.group))
        ## add means?
        #TODO: Safer way to check?
        if (any(srmMoments$sample.mean$group != 0)) {
          sat.mod <- c(sat.mod, paste(ov.names[[b]], "~ 1"))
          bas.mod <- c(bas.mod, paste(ov.names[[b]], "~ 1"))
        }


        ## case-level baseline model needs generalized reciprocity
      } else if (fitMore && component[b] == "case") {
        ## start with all variables
        sat.case <- outer(ov.names[[b]], ov.names[[b]], paste, sep = "~~")
        sat.mod <- c(sat.mod, sat.case[lower.tri(sat.case, diag = TRUE)])

        ## only use diagonal for baseline model ...
        bas.mod <- c(bas.mod, diag(sat.case))
        ## ... but add generalized reciprocity per RR variable
        RRnm <- names(data@varNames$RR)  # loop over all RR variables
        for (v in seq_along(RRnm)) {
          ## if both are modeled, correlate them
          v_ij <- paste(RRnm[v], c("out","in"), sep = "_")
          if (all(v_ij %in% ov.names[[b]])) {
            tmp <- paste(v_ij[[1]], "~~", v_ij[[2]])
            bas.mod <- c(bas.mod, tmp)
          }
        }


        ## dyad-level models need dyadic reciprocity + equality constraints
      } else if (fitMore && component[b] == "dyad") {
        ## loop over all variables, check whether to skip each
        RRnm <- names(data@varNames$RR)
        #TODO: need to concatenate $dyad, equate its covariance with each ij/ji

        ## start with all RR variables
        for (v1 in seq_along(RRnm)) {
          nm1_ij <- paste(RRnm[v1], "ij", sep = "_")
          nm1_ji <- paste(RRnm[v1], "ji", sep = "_")

          for (v2 in v1:length(RRnm)) {
            nm2_ij <- paste(RRnm[v2], "ij", sep = "_")
            nm2_ji <- paste(RRnm[v2], "ji", sep = "_")

            if (v1 == v2) {
              ## equate variances
              if (nm1_ij %in% ov.names[[b]]) {
                tmp <- paste0(nm1_ij, " ~~ var", v1, "*", nm1_ij)
                bas.mod <- c(bas.mod, tmp)
                sat.mod <- c(sat.mod, tmp)
              }
              if (nm1_ji %in% ov.names[[b]]) {
                tmp <- paste0(nm1_ji, " ~~ var", v1, "*", nm1_ji)
                bas.mod <- c(bas.mod, tmp)
                sat.mod <- c(sat.mod, tmp)
              }

              ## dyadic reciprocity
              if (nm1_ij %in% ov.names[[b]] && nm1_ji %in% ov.names[[b]]) {
                tmp <- paste0(nm1_ij, " ~~ ", nm1_ji)
                bas.mod <- c(bas.mod, tmp)
                sat.mod <- c(sat.mod, tmp)
              }

            } else {
              #FIXME?  if (!dots$h1) next

              ## inter-cor
              if (nm1_ij %in% ov.names[[b]] && nm2_ji %in% ov.names[[b]]) {
                tmp <- paste0(nm1_ij, " ~~ inter", v1, v2, "*", nm2_ji)
                sat.mod <- c(sat.mod, tmp)
              }
              if (nm2_ij %in% ov.names[[b]] && nm1_ji %in% ov.names[[b]]) {
                tmp <- paste0(nm2_ij, " ~~ inter", v1, v2, "*", nm1_ji)
                sat.mod <- c(sat.mod, tmp)
              }

              ## intra-cor
              if (nm1_ij %in% ov.names[[b]] && nm2_ij %in% ov.names[[b]]) {
                tmp <- paste0(nm1_ij, " ~~ intra", v1, v2, "*", nm2_ij)
                sat.mod <- c(sat.mod, tmp)
              }
              if (nm1_ji %in% ov.names[[b]] && nm2_ji %in% ov.names[[b]]) {
                tmp <- paste0(nm1_ji, " ~~ intra", v1, v2, "*", nm2_ji)
                sat.mod <- c(sat.mod, tmp)
              }

            }

          } # end v2
        }   # end v1

      } # else, no other components



    } # end b

  } else stop('data= must be an object of class?mvSRM (or a "lavMoments" list)')


  ## assemble lavaan() call
  lavCall <- list(quote(lavaan::lavaan), model = model, data = srmMoments)
  lavCall <- c(lavCall, dots)
  if (!is.null(srmMoments$sample.mean)) lavCall$meanstructure <- TRUE
  ## need to specify custom h1 / baseline model(s) below?
  if (fitMore) {
    ## so they don't get fit to begin with
    lavCall$h1       <- FALSE
    lavCall$baseline <- FALSE
  }
  ## use wrapper function?
  if (is.null(MC$model.type)) {
    ## do nothing: lavCall[[1]] <- lavaan::lavaan # remains the same
  } else if (MC$model.type == "sem") {
    lavCall[[1]] <- quote(lavaan::sem)
  } else if (MC$model.type == "cfa") {
    lavCall[[1]] <- quote(lavaan::cfa)
  }

  ## fit target model
  fit <- eval(as.call(lavCall))
  ## overwrite lavaan's default @call
  fit@call <- MC #FIXME: does this cause problems in lavaan?

  if (fitMore && lavInspect(fit, "converged")) {
    ## fit saturated model regardless
    lavCall[[1]] <- quote(lavaan::lavaan)
    lavCall$model <- sat.mod
    fit@external$h1 <- eval(as.call(lavCall))

    ## also fit custom case/dyad-level baseline model?
    if (dots$baseline) {
      ## fit baseline model
      lavCall$model <- bas.mod
      fit@external$baseline.model <- eval(as.call(lavCall))
      ## also store h1 model in baseline@external
      fit@external$baseline.model@external$h1 <- fit@external$h1
      ## and vice versa (necessary?  takes up space)
      fit@external$h1@external$baseline.model <- fit@external$baseline.model
    }

  }

  fit
}


##' @rdname lavaan.srm
##' @export
cfa.srm <- function(model, data, point = "mean", ...) {
  mc <- match.call(expand.dots = TRUE)
  mc$model.type <- "cfa"
  mc[[1L]] <- quote(lavaan.srm::lavaan.srm)
  eval(mc, parent.frame())
}

##' @rdname lavaan.srm
##' @export
sem.srm <- function(model, data, point = "mean", ...) {
  mc <- match.call(expand.dots = TRUE)
  mc$model.type <- "sem"
  mc[[1L]] <- quote(lavaan.srm::lavaan.srm)
  eval(mc, parent.frame())
}



## -----------------
## Utility functions
## -----------------


## function to
