# `lavaan.srm`

This is an R package whose primary purpose is to facilitate fitting a structural equation model (SEM) to round-robin data.  Estimation is accomplished via a two-stage estimation algorithm:

  1. Fit a multivariate social-relations model (SRM) to round-robin data, in order to estimate covariance matrices for person-level effects (e.g., actor and partner variances, generalized reciprocity) and for dyad-level effects (e.g., relationship variance, dyadic reciprocity).
  2. Fit a SEM using the summary statistics from Step 1 as data.

Step 1 is accomplished with Markov chain Monte Carlo estimation via the R package `rstan` (Carpenter et al., 2017; [web site](https://mc-stan.org/)). 

Step 2 is accomplished with the R package `lavaan` (Rosseel, 2012; [web site](https://lavaan.ugent.be/)).  Uncertainty about Stage-1 estimates is accounted for in Stage 2 using the multiple-imputation methods developed by Li Cai and colleagues ([2012](https://doi.org/10.3102/1076998612458320), [2019](https://doi.org/10.1080/00273171.2018.1523000)).




## Connecting `lavaan.mi` to `lavaan`

The top 4 lines of the table below indicate the `lavaan.mi::lavaan.mi()` function corresponds to the `lavaan::lavaan()` function, as do the wrappers `cfa.mi()`, `sem.mi()`, and `growth.mi()`, all of which return an object of class `lavaan.mi`.


|   `lavaan` Function    |   `lavaan.mi`    Function |
|:-----------------------|:--------------------------|
|      `lavaan()`        |       `lavaan.srm()`      |
|        `cfa()`         |         `cfa.srm()`       |
|        `sem()`         |         `sem.srm()`       |
| `parameterEstimates()` |        `summary()`        |
|    `lavTestLRT()`      |     `lavTestLRT.mi()`     |
|    `lavTestWald()`     |     `lavTestWald.mi()`    |
|    `lavTestScore()`    |     `lavTestScore.mi()`   |
|    `modIndices()`      |     `modIndices.mi()`     |

The `data=` argument must be a `list` of `data.frame`s rather than a single `data.frame`, and the specified `lavaan::model.syntax` will be applied to each `data.frame` in the `list`.


## Inheritance from `lavaanList`

The `lavaan.mi` class inherits from class `lavaanList` (from the `lavaan` package), extending it with a few additional and customized slots (see `class?lavaanList` and `class?lavaan.mi` for details).


## Methods for `lavaan` and `lavaan.mi` objects

Many remaining methods written for `lavaan`-class objects are also available for `lavaan.mi`-class objects (also see `class?lavaan.mi` for a list and descriptions of arguments):

- `show()`
- `fitMeasures()`
- `anova()`
- `nobs()`
- `coef()`
- `vcov()`
- `fitted()`
- `residuals()`


## Connecting `lavaan.mi` to `semTools`

Many of the functions in the `semTools` package continue to be available for `lavaan.mi`-class objects (e.g., `semTools::reliability()`), and functions such as `semTools::parcelAllocation()` and `semTools::plausibleValues()` are enhanced by allowing their output to be analyzed using `lavaan.mi()`.


