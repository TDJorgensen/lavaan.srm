### Terrence D. Jorgensen
### Last updated: 20 June 2024
### verify the hidden custom @external$h1.model works appropriately
### to yield test stats with the correct df AND scale/shift parameters

library(lavaan.srm)

TESTS <- c("satorra.bentler","scaled.shifted","browne.residual.adf")

checkTests  <- function(obj) {
  cat("UNCORRECTED:\n\n")
  lavaan:::lav_test_print(obj@test, nd = 5)

  cat("\n\nCORRECTED:\n\n")
  print(summary(obj, nd = 5, estimates = FALSE))

  cat("\n\nSTANDARD\n")
  print(lavTestLRT(obj, method = "standard"))

  cat("\n\nSCALED\n")
  sb <- lavTestLRT(obj, test = "satorra.bentler")
  print(sb)
  print(c(scale = attr(sb, "scale")[-1]))

  cat("\n\nSCALED and SHIFTED\n")
  sh <- lavTestLRT(obj, test = "scaled.shifted")
  print(sh)
  print(c(scale = attr(sh, "scale")[-1], shift = attr(sh, "shift")[-1]))

  cat("\n")
  print(lavTestLRT(obj, type = "browne.residual.adf"))
  return(invisible(obj))
}

## -------------------------------
## Example from help page
## -------------------------------

data(data.srm01, package="srm")
## do NOT use the same case-level IDs across groups
data.srm01$egoID   <- paste(data.srm01$Group, data.srm01$Actor  , sep = "_")
data.srm01$alterID <- paste(data.srm01$Group, data.srm01$Partner, sep = "_")

srmOut <- mvsrm(data = data.srm01, rr.vars = paste0("Wert", 1:3),
                IDout = "egoID", IDin = "alterID", IDgroup = "Group",
                ## rstan::sampling() arguments to run this quickly
                chains = 2, iter = 20, seed = 12345,
                cores = ifelse(parallel::detectCores() < 3L, 1L, 2L))


## STAGE 2: Specify and fit an SEM


## specify DYAD-level model
modD <- ' # equal loadings
  f1_ij =~ L1*Wert1_ij + L2*Wert2_ij + L3*Wert3_ij
  f1_ji =~ L1*Wert1_ji + L2*Wert2_ji + L3*Wert3_ji

## equal variances
  f1_ij ~~ phi*f1_ij
  f1_ji ~~ phi*f1_ji

  Wert1_ij ~~ v1*Wert1_ij
  Wert2_ij ~~ v2*Wert2_ij
  Wert3_ij ~~ v3*Wert3_ij

  Wert1_ji ~~ v1*Wert1_ji
  Wert2_ji ~~ v2*Wert2_ji
  Wert3_ji ~~ v3*Wert3_ji

## correlated residuals
  Wert1_ij ~~ Wert1_ji
  Wert2_ij ~~ Wert2_ji
  Wert3_ij ~~ Wert3_ji
'

## fit CFA to DYAD-level components,
## using mvSRM-class object (Stage 1) as input data
fitD <- lavaan.srm(modD, data = srmOut, component = "dyad",
                   ## more arguments passed to lavaan()
                   test = TESTS,
                   std.lv = TRUE,        # fix factor variances == 1
                   auto.cov.lv.x = TRUE) # estimate factor correlation
checkTests(fitD)

summary(fitD, std = TRUE, rsq = TRUE, fit.measures = TRUE,
        ## Ignore standard test, only use Browne's residual-based test.
        ## Also use Browne's residual-based test to calculate fit indices:
        fm.args = list(standard.test = "browne.residual.adf"))



## specify PERSON-level model
modP <- ' ## factor loadings
  f1_out =~ Wert1_out + Wert2_out + Wert3_out
  f1_in  =~ Wert1_in  + Wert2_in  + Wert3_in

## correlated residuals
  Wert1_out ~~ Wert1_in
  Wert2_out ~~ Wert2_in
  Wert3_out ~~ Wert3_in
'

## fit CFA to PERSON-level components,
## using wrapper analogous to lavaan::cfa()
fitP <- cfa.srm(modP, data = srmOut, component = "case", std.lv = TRUE,
                test = TESTS) # multiple test stats
checkTests(fitP)
## fitP is a lavaan-class object, so any lavaan function is available
summary(fitP, std = TRUE, rsq = TRUE, fit.measures = TRUE,
        ## Ignore standard test, only use Browne's residual-based test.
        ## Also use Browne's residual-based test to calculate fit indices:
        fm.args = list(standard.test = "browne.residual.adf"))



## Simultaneously specifying BOTH person & dyad-level models
## requires block-structured syntax, analogous to multilevel SEM.

## specify a model with cross-level measurement equivalence
modPD <- ' group: 1   # eventually use "component: case"

  f1_out =~ L1*Wert1_out + L2*Wert2_out + L3*Wert3_out
  f1_in  =~ L1*Wert1_in  + L2*Wert2_in  + L3*Wert3_in

## correlated residuals
  Wert1_out ~~ Wert1_in
  Wert2_out ~~ Wert2_in
  Wert3_out ~~ Wert3_in

## free factor variances (use dyad-level = 1 as reference "group")
  f1_out ~~ var_perceiver*f1_out
  f1_in  ~~ var_target*f1_in


group: 2   # eventually use "component: dyad"

  f1_ij =~ L1*Wert1_ij + L2*Wert2_ij + L3*Wert3_ij
  f1_ji =~ L1*Wert1_ji + L2*Wert2_ji + L3*Wert3_ji

## equal variances
  f1_ij ~~ var_relationship*f1_ij
  f1_ji ~~ var_relationship*f1_ji

  Wert1_ij ~~ v1*Wert1_ij
  Wert2_ij ~~ v2*Wert2_ij
  Wert3_ij ~~ v3*Wert3_ij

  Wert1_ji ~~ v1*Wert1_ji
  Wert2_ji ~~ v2*Wert2_ji
  Wert3_ji ~~ v3*Wert3_ji

## correlated residuals
  Wert1_ij ~~ Wert1_ji
  Wert2_ij ~~ Wert2_ji
  Wert3_ij ~~ Wert3_ji


## Identification constraint: Sum of factor variances == 1
  var_relationship == 1 - (var_perceiver + var_target)
## Thus, each estimated variance == proportion of total (ICC)
'
## order in component= must match order in model syntax above
fitPD <- lavaan.srm(modPD, data = srmOut, component = c("case","dyad"),
                    auto.var = TRUE, auto.cov.lv.x = TRUE, test = TESTS)
checkTests(fitPD)
summary(fitPD, std = TRUE, rsq = TRUE, fit.measures = TRUE,
        ## Ignore standard test, only use Browne's residual-based test.
        ## Also use Browne's residual-based test to calculate fit indices:
        fm.args = list(standard.test = "browne.residual.adf"))



## -------------------------------
## Example from Aditi's simulation
## -------------------------------

s1ests <- readRDS("s1ests.rds") # in the hidden lavaan.srm/.tests/ directory

mod_combi <- ' group: 1
  Factor_out =~ 1*V1_out + V2_out + V3_out
  Factor_in =~ 1*V1_in + V2_in + V3_in

  Factor_out ~~ Factor_out + Factor_in
  Factor_in ~~ Factor_in

  V1_out ~~ V1_out + V1_in
  V1_in ~~ V1_in
  V2_out ~~ V2_out
  V2_in ~~ V2_in
  V3_out ~~ V3_out + V3_in
  V3_in ~~ V3_in

  group: 2
  Factor_ij =~ 1*V1_ij + FL2*V2_ij + FL3*V3_ij
  Factor_ji =~ 1*V1_ji + FL2*V2_ji + FL3*V3_ji

  Factor_ij ~~ Fvar*Factor_ij + Factor_ji
  Factor_ji ~~ Fvar*Factor_ji

  V1_ij ~~ Ivar1*V1_ij
  V1_ji ~~ Ivar1*V1_ji
  V2_ij ~~ Ivar2*V2_ij + V2_ji
  V2_ji ~~ Ivar2*V2_ji
  V3_ij ~~ Ivar3*V3_ij + V3_ji
  V3_ji ~~ Ivar3*V3_ji
  '

fit_combi <- lavaan.srm(model = mod_combi, data = s1ests,
                        component = c("case", "dyad"), posterior.est = "mean",
                        test = TESTS)
checkTests(fit_combi)

mod_case <- '
  Factor_out =~ 1*V1_out + V2_out + V3_out
  Factor_in =~ 1*V1_in + V2_in + V3_in

  Factor_out ~~ Factor_out + Factor_in
  Factor_in ~~ Factor_in

  V1_out ~~ V1_out + V1_in
  V1_in ~~ V1_in
  V2_out ~~ V2_out
  V2_in ~~ V2_in
  V3_out ~~ V3_out + V3_in
  V3_in ~~ V3_in
  '

fit_case <- lavaan.srm(model = mod_case, data = s1ests,
                       component = "case", posterior.est = "mean",
                       test = TESTS)
checkTests(fit_case)

mod_dyad <- '
  Factor_ij =~ 1*V1_ij + FL2*V2_ij + FL3*V3_ij
  Factor_ji =~ 1*V1_ji + FL2*V2_ji + FL3*V3_ji

  Factor_ij ~~ Fvar*Factor_ij + Factor_ji
  Factor_ji ~~ Fvar*Factor_ji

  V1_ij ~~ Ivar1*V1_ij
  V1_ji ~~ Ivar1*V1_ji
  V2_ij ~~ Ivar2*V2_ij + V2_ji
  V2_ji ~~ Ivar2*V2_ji
  V3_ij ~~ Ivar3*V3_ij + V3_ji
  V3_ji ~~ Ivar3*V3_ji
  '

fit_dyad <- lavaan.srm(model = mod_dyad, data = s1ests, component = "dyad", posterior.est = "mean",
                       test = TESTS)
checkTests(fit_dyad)
