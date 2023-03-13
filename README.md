# Installing `lavaan.srm`

This is an R package whose primary purpose is to facilitate fitting a structural equation model (SEM) to round-robin data.  Installation via the `remotes` package requires a C++ compiler, which can be enabled for 

- Windows OS by installing [Rtools](https://cran.r-project.org/bin/windows/Rtools/)
- Mac OS by installing Xcode / Command Line Tools (open Applications/Utilities/Terminal and run `xcode-select -—install`)

After restarting your computer, open R, install the `remotes` package (if necessary), and use it to install `lavaan.srm`:

```
# install.packages("remotes")
remotes::install_github("TDJorgensen/lavaan.srm")
```


# Estimation algorithm

Estimation is accomplished via a two-stage estimation algorithm:

  1. Fit a multivariate social-relations model (SRM) to round-robin data, in order to estimate covariance matrices for person-level effects (e.g., actor and partner variances, generalized reciprocity) and for dyad-level effects (e.g., relationship variance, dyadic reciprocity).
  2. Fit a SEM using the summary statistics from Stage 1 as data.

Stage 1 is accomplished with Markov chain Monte Carlo estimation via the R package `rstan` (Carpenter et al., 2017; [web site](https://mc-stan.org/)). 

Stage 2 is accomplished with the R package `lavaan` (Rosseel, 2012; [web site](https://lavaan.ugent.be/)).  Uncertainty about Stage-1 estimates is accounted for in Stage 2 using methods developed by Li Cai and colleagues ([2012](https://doi.org/10.3102/1076998612458320), [2019](https://doi.org/10.1080/00273171.2018.1523000)) for multiple imputations of missing data.

This is an efficient alternative to the single-stage maximum-likelihood estimation of the social-relations SEM (SR-SEM) proposed by Nestler et al. ([2020](https://doi.org/10.1007/s11336-020-09728-z)), implemented in the R package `srm` (on [CRAN](https://cran.r-project.org/package=srm) and [GitHub](https://github.com/alexanderrobitzsch/srm)).  The flexible 2-stage estimator in `lavaan.srm` can easily be adapted to allow for numerous practical extensions:

- handling missing data in Stage 1 via data augmentation (analogous to multiple imputation) or case-specific likelihoods (analogous to full-information maximum likelihood, as implemented in `blavaan` when `target = "stan"`).
- introducing a Stage-1 threshold model for binary/ordinal variables
- including level-specific covariates of random effects estimated in Stage 1

These options will be explored, tested, and implemented in future versions of `lavaan.srm`. 



## Stage 1

Round-robin `data=` are passed to `mvsrm()`, along with case-level ID variables indicating the two members of the dyad.  Group IDs can also be provided to distinguish multiple round-robin groups.  The `data=` must be in long format, such that all observed values of a round-robin variable are stored in a single column, with 2 rows for each dyad.  An optional list of round-robin variable names in `data=` to be modeled, but all non-ID variables will be modeled by default.  Descriptions of further arguments and example syntax can be found on the `?mvsrm` help page.

### Extracting SRM results

The `mvsrm()` function fits a multivariate SRM using the Stan software, returning an object of `class?mvSRM` that inherits from `class?stanfit`.  Thus, any `methods(class = "stanfit")` and other `rstan` functions can be applied to the `mvsrm()` output, unless another method exists specifically for `mvSRM-class`.  For example, this package includes `summary()` and `as.matrix()` methods for `mvSRM-class` objects, but adding the argument `as.stanfit=TRUE` allows using the methods defined in the `rstan` package for `stanfit-class` objects.



## Stage 2

The top 4 lines of the table below indicate the `lavaan.srm::lavaan.srm()` function corresponding to the `lavaan::lavaan()` function---namely, the wrappers `cfa.srm()` and `sem.srm()`---all of which return an object of `class?lavaan`.


|   `lavaan` Function    |   `lavaan.srm`   Function |
|:-----------------------|:--------------------------|
|      `lavaan()`        |       `lavaan.srm()`      |
|        `cfa()`         |         `cfa.srm()`       |
|        `sem()`         |         `sem.srm()`       |

The `data=` argument must be an `mvSRM-class` object, as returned by the `mvsrm()` function.  The `component=` argument is required to inform `lavaan.srm()` about which level(s) of analysis to which the `model=` syntax applies.  An optional `posterior.est=` argument can be used to choose a posterior point-estimate of each SRM (co)variance parameter; the default is the posterior mean.  Many other `lavaan()` or `lavOptions()` arguments can be passed to `lavaan.srm()` via `...` (e.g., `std.lv = TRUE`).


### Extracting SR-SEM results

Because `lavaan.srm()` returns a `lavaan-class` object, any relevant `lavaan` functions can be used, as long as the method is relevant for models fitted to summary statistics (i.e., `sample.cov=` rather than raw `data=`).  This is because `lavaan.srm()` internally calculates the relevant `sample.cov=`, `sample.nobs=`, and `NACOV=` arguments to pass to `lavaan()`, so raw/casewise `data=` are not actually analyzed.  Therefore, robust corrections for nonnormality cannot be requested (e.g., `estimator = "MLM"`).  Likewise, many functions from other packages (e.g., `semTools`) are also available for `lavaan.srm()` output, such as `semTools::compRelSEM()`.


### Model tests and fit indices

The `lavaan.srm()` function calls `lavaan()` with some specific `lavOptions()` to provide reasonable defaults.  For example, `lavaan` always provides the "standard" test statistic of data--model fit in its `summary()` output.  But that statistic does **not** provide nominal Type I error rates, so the `lavaan.srm()` function sets `test = "Browne.residual.adf"` to use a statistic that does asymptotically yield nominal Type I error rates.  Thus, in the `summary()` output, one should pay attention only to the test statistic labeled "Browne's residual-based (ADF) test".  Browne's test will be the default returned by `lavTestLRT()`, which is called by `anova()`.

However, `fitMeasures()` still calculates fit indices by default with the standard test.  To choose the better test statistic (using Browne's ADF formula) for calculating fit indices, you can set `fitMeasures(..., fm.args = list(standard.test = "Browne.residual.adf"))`.  Likewise, when calling `summary(..., fit.measures = TRUE)`, one should also add the argument `fm.args = list(standard.test = "Browne.residual.adf")` to `summary()` so that the fit indices are calculated using Browne's test statistic.

Some specialized dyad- and case-level models are automatically fitted to represent the null model (for incremental fit indices) and unrestricted (`h1`) models, the latter of which is necessary for an appropriate dyad-level test of data--model fit.  These specialized models are described in the next sections.


### Specialized "saturated" dyad-level model

There are no implicit equality constraints in the case-level model, so the standard case-level saturated model is appropriate for testing.  But dyad-level components (i.e., with `_ij` and `_ji` suffix) are "symmetric" for indistinguishable dyads, so their variances should be constrained to equality for each round-robin variable.  Likewise, between-variable covariances should be constrained to equality for within-person / intrapersonal components (i.e., `_ij` with `_ij`) and for between-person / interpersonal components (i.e., `_ij` with `_ji`).  See Nestler ([2018](https://doi.org/10.3102/1076998617741106), p. 390) for details.

The standard saturated model places no such constraints on the implicitly equal (co)variances.  Although their point estimates would be unaffected, their associated 
standard errors (*SE*s) would be larger than necessary because every such equal parameter would be estimated twice rather than once.  If the hypothesized model does not specify these equality constraints, then its *SE*s would likewise be too large.

The same concern arises regarding the fitted model's test statistic.  Specifying the implicit equality constraints in the fitted model would appropriately decrease its *SE*s, as well as increase its degrees of freedom (*df*).  However, the test statistic would be unchanged because the parameter estimates themselves would be unaffected by the equality constraints.  The same holds for the saturated model, whose $df=0$ would yield the same perfect fit as an otherwise unconstrained model that imposes the same dyad-level equality constraints (i.e., whose $df>0$).  Therefore, the fitted and saturated models should be synchronized---that is, the test statistic of a dyad-level SEM with equality constraints appropriate for indistinguishable dyads should be calculated by comparing it to an unrestricted (but not fully saturated) model with analogous equality constraints.  In other words, the test statistic reflects misspecification due only to the hypothesized (not implicit) constraints, so the $df$ of its test statistic should only reflect the number of hypothesized constraints. 

A synchronized `h1` model is stored in the `lavaan.srm()` output whenever modeling dyad-level data, unless it is turned off by setting `lavaan.srm(..., h1 = FALSE)`.  It is stored in the slot `@external$h1`, where it is detected and used by `lavTestLRT()` and `fitMeasures()`.  



### Specialized "null" models for case-level and dyad-level

Likewise, the null model used to calculate incremental fit indices (e.g., CFI and TLI) should also be synchronized with the fitted and saturated models, as described above.  By default, `lavaan::fitMeasures()` uses the "independence model" when the `baseline.model=` argument is not specified, which sets all covariances to zero and estimates all variances. Instead, the same implicit equality constraints should be imposed on variances of each round-robin variable when dyads are indistinguishable.  

Because the `_ij` and `_ji` components are measured on the same variable, their residuals should be correlated in any fitted model.  Thus, the dyadic reciprocity in the null model should also be estimated, as no reasonable null hypothesis would assume dyadic reciprocity is zero.  There are no implicit equality constraints in the case-level model, but it would also be unreasonable to assume generalized reciprocity is zero, so the null model should also estimate each generalized reciprocity.

A synchronized null model is stored in the `lavaan.srm()` output whenever modeling case- or dyad-level data, unless it is turned off by setting `lavaan.srm(..., baseline = FALSE)`.  It is stored in the slot `@external$baseline.model`, where it is detected and used by `fitMeasures()`. When modeling dyad-level data, the synchronized `h1` model is also stored in the slot `@external$baseline.model@external$h1`, so the null model's fit statistic has the appropriate *df*.  Note that users can always specify and fit their own custom null model to be passed to `fitMeasures()`, which might need to reflect further constraints (e.g., consistent with evaluating measurement invariance; see Widamin & Thompson, [2003](https://psycnet.apa.org/doi/10.1037/1082-989X.8.1.16)).


