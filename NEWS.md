---
output:
  html_document: default
  pdf_document: default
---
`MatchIt` News and Updates
======

# MatchIt 4.1.0

* Coarsened exact matching (i.e., `matchit()` with `method = "cem"`) has been completely rewritten and no longer involves the `cem` package, eliminating some spurious warning messages and fixing some bugs. All the same arguments can still be used, so old code will run, though some results will differ slightly. Additional options are available for matching and performance has improved. See `?method_cem` for details on the differences between the implementation in the current version of `MatchIt` and that in `cem` and older versions of `MatchIt`. In general, these changes make coarsened exact matching function as one would expect it to, circumventing some peculiarities and bugs in the `cem` package.

* Variable ratio matching is now compatible with `method = "optimal"` in the same way it is with `method = "nearest"`, i.e., by using the `min.controls` and `max.controls` arguments.

* With `method = "full"` and `method = "optimal"`, the maximum problem size has been set to unlimited, so that larger datasets can be used with these methods without error. They may take a long time to run, though.

* Processing improvements with `method = "optimal"` due to rewriting some functions in `Rcpp`.

* Using `method = "optimal"` runs more smoothly when combining it with exact matching through the `exact` argument.

* When using `ratio` different from 1 with `method = "nearest"` and `method = "optimal"` and with exact matching, errors and warnings about the number of units that will be matched are clearer. Certain `ratio`s that would produce errors now only produce warnings.

* Fixed a bug when no argument was supplied to `data` in `matchit()`.

* Improvements to vignettes and documentation.

# MatchIt 4.0.1

* Restored `cem` functionality after it had been taken down and re-uploaded.

* Added `pkgdown` website.

* Computing matching weights after matching with replacement is faster due to programming in `Rcpp`.

* Fixed issues with `Rcpp` code that required C++11. C++11 has been added to SystemRequirements in DESCRIPTION, and `MatchIt` now requires R version 3.1.0 or later.

# MatchIt 4.0.0

## General Fixes and New Features

* `match.data()`, which is used to create matched datasets, has a few new arguments. The `data` argument can be supplied with a dataset that will have the matching weights and subclasses added. If not supplied, `match.data()` will try to figure out the appropriate dataset like it did in the past. The `drop.unmatched` argument controls whether unmatched units are dropped from the output. The default is `TRUE`, consistent with past behavior. Warnings are now more informative.

* `get_matches()`, which seems to have been rarely used since it performed a similar function to `match.data()`, has been revamped. It creates a dataset with one row per unit per matched pair. If a unit is part of two separate pairs (e.g., as a result of matching with replacement), it will get two rows in the output dataset. The goal here was to be able to implement standard error estimators that rely both on repeated use of the same unit and subclass/pair membership, e.g., Austin & Cafri (2020). Otherwise, it functions similarly to `match.data()`. *NOTE: the changes to `get_matches()` are breaking changes! Legacy code will not work with the new syntax!*

* `print.matchit()` has completely changed and now prints information about the matching type and specifications. `summary.matchit()` contains all the information that was in the old `print` method.

* A new function, `add_s.weights()`, adds sampling weights to `matchit` objects for use in balance checking and effect estimation. Sampling weights can also be directly supplied to `matchit()` through the new `s.weights` argument. A new vignette describing how to using `MatchIt` with sampling weights is available at `vignette("sampling-weights")`.

* The included dataset, `lalonde`, now uses a `race` variable instead of separate `black` and `hispan` variables. This makes it easier to see how character variables are treated by `MatchIt` functions.

* Added extensive documentation for every function, matching method, and distance specification. Documentation no longer links to `gking.harvard.edu/matchit` as it now stands alone.

## `matchit()`
* An argument to `data` is no longer required if the variables in `formula` are present in the environment.

* When missing values are present in the dataset but not in the treatment or matching variables, the error that used to appear no longer does.

* The `exact` argument can be supplied either as a character vector of names of variables in `data` or as a one-sided formula. A full cross of all included variables will be used to create bins within which matching will take place.

* The `mahvars` argument can also be supplied either as a character vector of names of variables in `data` or as a one-sided formula. Mahalanobis distance matching will occur on the variables in the formula, processed by `model.matrix()`. Use this when performing Mahalanobis distance matching on some variables within a caliper defined by the propensity scores estimated from the variables in the main `formula` using the argument to `distance`. For regular Mahalanobis distance matching (without a propensity score caliper), supply the variables in the main `formula` and set `distance = "mahalanobis"`.

* The `caliper` argument can now be specified as a numeric vector with a caliper for each variable named in it. This means you can separately impose calipers on individual variables as well as or instead of the propensity score. For example, to require that units within pairs must be no more than .2 standard deviations of `X1` away from each other, one could specify `caliper = c(X1 = .2)`. A new option `std.caliper` allows the choice of whether the caliper is in standard deviation units or not, and one value per entry in `caliper` can be supplied. An unnamed entry to `caliper` applies the caliper to the propensity score and the default of `std.caliper` is `FALSE`, so this doesn't change the behavior of old code. These options only apply to the methods that accept calipers, namely `"nearest"`, `"genetic"`, and `"full"`.

* A new `estimand` argument can be supplied to specify the target estimand of the analysis. For all methods, the ATT and ATC are available with the ATT as the default, consistent with prior behavior. For some methods, the ATE is additionally available. Note that setting the estimand doesn't actually mean that estimand is being targeted; if calipers, common support, or other restrictions are applied, the target population will shift from that requested. `estimand` just triggers the choice of which level of the treatment is focal and what formula should be used to compute weights from subclasses.

* In methods that accept it, `m.order` can be set to "`data`", which matches in the order the data appear. With `distance = "mahalanobis"`, `m.order` can be "`random`" or "`data`", with "`data`" as the default. Otherwise, `m.order` can be `"largest"`, `"smallest"`, `"random"`, or `"data"`, with `"largest"` as the default (consistent with prior behavior).

* The output to `matchit()` has changed slightly; the component `X` is now a data frame, the result of a call to `model.frame()` with the formula provided. If `exact` or `mahvars` are specified, their variables are included as well, if not already present. It is included for all methods and is the same for all methods. In the past, it was the result of a call to `model.matrix()` and was only included for some methods.

* When key arguments are supplied to methods that don't accept them, a warning will be thrown.

* `method` can be set to `NULL` to not perform matching but create a `matchit` object, possibly with a propensity score estimated using `distance` or with a common support restriction using `discard`, for the purpose of supplying to `summary.matchit()` to assess balance prior to matching.

### `method = "nearest"`

* Matching is much faster due to re-programming with `Rcpp`.

* With `method = "nearest"`, a `subclass` component containing pair membership is now included in the output when `replace = FALSE` (the default), as it has been with optimal and full matching.

* When using `method = "nearest"` with `distance = "mahalanobis"`, factor variables can now be included in the main `formula`. The design matrix no longer has to be full rank because a generalized inverse is used to compute the Mahalanobis distance.

* Unless `m.order = "random"`, results will be identical across runs. Previously, several random choices would occur to break ties. Ties are broken based on the order of the data; shuffling the order of the data may therefore yield different matches.

* When using `method = "nearest"` with a caliper specified, the nearest control unit will be matched to the treated unit if one is available. Previously, a random control unit within the caliper would be selected. This eliminates the need for the `calclosest` argument, which has been removed.

* Variable ratio extremal matching as described by Ming & Rosenbaum (2000) can be implemented using the new `min.controls` and `max.controls` arguments.

* Added ability to display a progress bar during matching, which can be activated by setting `verbose = TRUE`.

### `method = "optimal"` and `method = "full"`

* Fixed bug in `method = "optimal"`, which produced results that did not match `optmatch`. Now they do.

* Added support for optimal and full Mahalanobis distance matching by setting `method = "mahalanobis"` with `method = "optimal"` and `method = "full"`. Previously, both methods would perform a random match if `method` was set to `"mahalanobis"`. Now they use the native support in `optmatch::pairmatch()` and `optmatch::fullmatch()` for Mahalanobis distance matching.

* Added support for exact matching with `method = "optimal"` and `method = "full"`. As with `method = "nearest"`, the names of the variables for which exact matches are required should be supplied to the `exact` argument. This relies on `optmatch::exactMatch()`.

* The warning that used to occur about the order of the match not guaranteed to be the same as the original data no longer occurs.

* For `method = "full"`, the `estimand` argument can be set to `"ATT"`, `"ATC"`, or `"ATE"` to compute matching weights that correspond to the given estimand. See `?matchit` for details on how weights are computed for each `estimand`.

### `method = "genetic"`

* Fixed a bug with `method = "genetic"` that caused an error with some `ratio` greater than 1.

* The default of `replace` in `method = "genetic"` is now `FALSE`, as it is with `method = "nearest"`.

* When `verbose = FALSE`, the default, no output is printed with `method = "genetic"`. With `verbose = TRUE`, the printed output of `Matching::GenMatch()` with `print.level = 2` is displayed.

* The `exact` argument now correctly functions with `method = "genetic"`. Previously, it would have to be specified in accordance with its use in `Matching::GenMatch()`.

* Different ways to match on variables are now allowed with `method = "genetic"`, similar to how they are with `method = "nearest"`. If `distance = "mahalanobis"`, no propensity score will be computed, and genetic matching will be performed just on the variables supplied to `formula`. If `mahvars` is specified, genetic matching will be performed on the variables supplied to `mahvars`, but balance will be optimized on all covariates supplied to `formula`. Otherwise, genetic matching will be performed on the variables supplied to `formula` and the propensity score. Previously, `mahvars` was ignored. Balance is now always optimized on the variables included in `formula` and never on the propensity score, whereas in the past the propensity score was always included in the balance optimization.

* The `caliper` argument now works as it does with `method = "nearest"` and other methods rather than needing to be supplied in a way that `Matching::Match()` would accept.

* A `subclass` component is now included in the output when `replace = FALSE` (the default), as it has been with optimal and full matching.

### `method = "cem"` and `method = "exact"`

* With `method = "cem"`, the `k2k` argument is now recognized. Previously it was ignored unless an argument to `k2k.method` was supplied.

* The `estimand` argument can be set to `"ATT"`, `"ATC"`, or `"ATE"` to compute matching weights that correspond to the given estimand. Previously only ATT weights were computed. See `?matchit` for details on how weights are computed for each `estimand`.

### `method = "subclass"`

* Performance improvements.

* A new argument, `min.n`, can be supplied, which controls the minimum size a treatment group can be in each subclass. When any estimated subclass doesn't have enough members from a treatment group, units from other subclasses are pulled to fill it so that every subclass will have at least `min.n` units from each treatment group. This uses the same mechanism as is used in `WeightIt`. The default `min.n` is 1 to ensure there are at least one treated and control unit in each subclass.

* Rather than producing warnings and just using the default number of subclasses (6), when an inappropriate argument is supplied to `subclass`, an error will occur.

* The new `subclass` argument to `summary()` can be used to control whether subclass balance statistics are computed; it can be `TRUE` (display balance for all subclasses), `FALSE` (display balance for no subclasses), or a vector of subclass indices on which to assess balance. The default is `FALSE`.

* With `summary()`, balance aggregating across subclasses is now computed using subclass weights instead of by combining the subclass-specific balance statistics.

* The `sub.by` argument has been replaced with `estimand`, which can be set to `"ATT"`, `"ATC"`, or `"ATE"` to replace the `sub.by` inputs of `"treat"`, `"control"`, and `"all"`, respectively. Previously, weights for `sub.by` that wasn't `"treat"` were incorrect; they are now correctly computed for all inputs to `estimand`.

### `distance`

* The allowable options to `distance` have changed slightly. The input should be either `"mahalanobis"` for Mahalanobis distance matching (without a propensity score caliper), a numeric vector of distance values (i.e., values whose absolute pairwise differences form the distances), or one of the allowable options. The new allowable values include `"glm"` for propensity scores estimated with `glm()`, `"gam"` for propensity scores estimated with `mgcv::gam()`, `"rpart"` for propensity scores estimated with `rpart::rpart()`, `"nnet"` for propensity scores estimated with `nnet::nnet()`, `"cbps"` for propensity scores estimated with `CBPS::CBPS()`, or `bart` for propensity scores estimated with `dbarts::bart2()`. To specify a link (e.g., for probit regression), specify an argument to the new `link` parameter. For linear versions of the propensity score, specify `link` as `"linear.{link}"`. For example, for linear probit regression propensity scores, one should specify `distance = "glm", link = "linear.probit"`. The default `distance` is `"glm"` and the default link is `"logit"`, so these can be omitted if either is desired. Not all methods accept a `link`, and for those that don't, it will be ignored. If an old-style `distance` is supplied, it will be converted to an appropriate specification with a warning (except for `distance = "logit"`, which will be converted without a warning).

* Added `"cbps"` as option for `distance`. This estimates propensity scores using the covariate balancing propensity score (CBPS) algorithm as implemented in the `CBPS` package. Set `link = "linear"` to use a linear version of the CBPS.

* Added `"bart"` as an option for `distance`. This estimates propensity scores using Bayesian Additive Regression Trees (BART) as implemented in the `dbarts` package.

* Added `"randomforest"` as an option for `distance`. This estimates propensity scores using random forests as implemented in the `randomForest` package.

* Bugs in `distance = "rpart"` have been fixed.

## `summary.matchit()`

* When `interactions = TRUE`, interactions are no longer computed with the distance measure or between dummy variables of the same factor. Variable names are cleaned up and easier to read.

* The argument to `addlvariables` can be specified as a data frame or matrix of covariates, a formula with the additional covariates (and transformations) on the right side, or a character vector containing the names of the additional covariates. For the latter two, if the variables named do not exist in the `X` component of the `matchit` output object or in the environment, an argument to `data` can be supplied to `summary()` that contains these variables.

* The output for `summary()` is now the same for all methods (except subclassification). Previously there were different methods for a few different types of matching. 

* The eCDF median (and QQ median) statistics have been replaced with the variance ratio, which is better studied and part of several sets of published recommendations. The eCDF and QQ median statistics provide little information above and beyond the corresponding mean statistics. The variance ratio uses the variances weighted by the matching weights.

* The eCDF and QQ statistics have been adjusted. Both now use the weights that were computed as part of the matching. The eCDF and QQ statistics for binary variables are set to the difference in group proportions. The standard deviation of the control group has been removed from the output.

* The default for `standardize` is now `TRUE`, so that standardized mean differences and eCDF statistics will be displayed by default.

* A new column for the average absolute pair difference for each covariate is included in the output. The values indicate how far treated and control units within pairs are from each other. An additional argument to `summary.matchit()`, `pair.dist`, controls whether this value is computed. It can take a long time for some matching methods and could be omitted to speed up computation.

* Balance prior to matching can now be suppressed by setting `un = FALSE`.

* Percent balance improvement can now be suppressed by setting `improvement = FALSE`. When `un = FALSE`, `improvement` is automatically set to `FALSE`.

## `plot.matchit()`

* Plots now use weighted summaries when weights are present, removing the need for the `num.draws` argument.

* Added a new plot type, `"ecdf"`, which creates empirical CDF plots before and after matching.

* The appearance of some plots has improved (e.g., text is appropriately centered, axes are more clearly labeled). For eQQ plots with binary variables or variables that take on only a few values, the plots look more like clusters than snakes.

* The argument to `type` can be abbreviated (e.g., `"j"` for jitter).

* Fixed a bug that caused all plots generated after using `plot(., type = "hist")` to be small.

* When specifying an argument to `which.xs` to control for which variables balance is displayed graphically, the input should be the name of the original variable rather than the version that appears in the `summary()` output. In particular, if a factor variable was supplied to `matchit()`, it should be referred to by its name rather than the names of its split dummies. This makes it easier to view balance on factor variables without having to know or type the names of all their levels.

* eQQ plots can now be used with all matching methods. Previously, attempting `plot()` after `method = "exact"` would fail.

## `plot.summary.matchit()`

* The summary plot has been completely redesigned. It is now a Love plot made using `graphics::dotchart()`. A few options are available for ordering the variables, presenting absolute or raw standardized mean differences, and placing threshold lines on the plots. For a more sophisticated interface, see `cobalt::love.plot()`, which natively supports `matchit` objects and uses `ggplot2` as its engine.
