### Terrence D. Jorgensen
### Last updated: 7 June 2024
### documentation page for the lavaan.srm package


##' The 'lavaan.srm' package.
##'
##' @description
##' This package provides a convenient interface for fitting a multivariate
##' social relations model (SRM; Nestler, 2018) using Markov chain Monte Carlo
##' (MCMC) estimation with the \pkg{rstan} package. This package also provides a
##' convenient interface for fitting a structural equation model (SEM) to the
##' SRM results, thus enabling a two-stage estimator of the SR-SEM (Nestler
##' et al., 2020, 2022) via the \pkg{lavaan} package (Rosseel, 2012).
##'
#TODO: delete this deprecated feature: @docType package

##' @name lavaan.srm-package
##' @aliases lavaan.srm-package
##' @useDynLib lavaan.srm, .registration = TRUE
##' @import methods
##' @import Rcpp
##'
##' @references
##' Nestler, S. (2018). Likelihood estimation of the multivariate social
##' relations model. *Journal of Educational and Behavioral Statistics, 43*(4),
##' 387--406. <https://doi.org/10.3102/1076998617741106>
##'
##' Nestler, S., Lüdtke, O., & Robitzsch, A. (2020). Maximum likelihood
##' estimation of a social relations structural equation model.
##' *Psychometrika, 85*(4), 870--889.
##' <https://doi.org/10.1007/s11336-020-09728-z>
##'
##' Nestler, S., Lüdtke, O., & Robitzsch, A. (2022). Analyzing longitudinal
##' social relations model data using the social relations structural equation
##' model. *Journal of Educational and Behavioral Statistics, 47*(2), 231--260.
##' <https://doi.org/10.3102/10769986211056541>
##'
##' Rosseel, Y. (2012). `lavaan`: An R package for structural equation modeling.
##' *Journal of Statistical Software, 48*(2), 1--36.
##' <https://doi.org/10.18637/jss.v048.i02>
##'
##' Stan Development Team. (2022). `rstan`: The R interface to Stan.
##' R package version 2.21.7. https://mc-stan.org
##'
NULL
