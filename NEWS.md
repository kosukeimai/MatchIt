`MatchIt` News and Updates
======

# MatchIt (development version)

## General Fixes

* When missing values are present in the dataset but not in the treatment or matching variables, the error that used to appear no longer does.

* An argument to `data` is no longer required if the variables in `formula` are present in the environment.

* The `exact` argument can be supplied either as a character vector of names of variables in `data` or as a one-sided formula. A full cross of all included variables will be used to create bins within which matching will take place.

* The `mahvars` argument can also be supplied either as a character vector of names of variables in `data` or as a one-sided formula. Mahalanobis distance matching will occur on the variables in the formula, processed by `model.matrix()`. Use this when performing Mahalanobis distance matching on some variables within a caliper defined by the propensity scores estimated from the variables in the main `formula` using the argument to `distance`. For regular Mahalanobis distance matching, supply the variables in the main `formula` and set `distance = "mahalanobis"`.

* The allowable options to `distance` have changed slightly. The input should be either `"mahalanobis"` for Mahalanobis distance matching (without a caliper), a numeric vector of distance values (i.e., values whose absolute pairwise differences form the distances), or one of the allowable options. The new allowable values include `"glm"` for propensity scores estimated with `glm()`, `"gam"` for propensity scores estimated with `mgcv::gam()`, `"rpart"` for propensity scores estimated with `rpart::rpart()`, `"nnet"` for propensity scores estimated with `nnet::nnet()`, or `"cbps"` for propensity scores estimated with `CBPS::CBPS()`. To specify a link (e.g., for probit regression), specify an argument to the new `link` parameter. For linear versions of the propensity score, specify `link` as `"linear.{link}"`. For example, for linear probit regression propensity scores, one should specify `distance = "glm", link = "linear.probit"`. The default `distance` is `"glm"` and the default link is `"logit"`, so these can be omitted if either is desired. Not all methods accept a `link`, and for these, it will be ignored. If an old-style `distance` is supplied, it will be converted to an appropriate specification with a warning (except for `distance = "logit"`, which will be converted without a warning).

* The output to `matchit()` has changed slightly; the component `X` is now a data frame, the result of a call to `get_all_vars()` with the formula provided. If `exact` or `mahvars` are specified, their variables are included as well, if not already present. It is included for all methods and is the same for all methods. In the past, it was the result of a call to `model.matrix()` and was only included for some methods.

* Added `"cbps"` as option for `distance`. This estimates propensity scores using the covariate balancing propensity score (CBPS) algorithm as implemented in the `CBPS` package. Set `link = "linear"` to use a linear version of the CBPS.

* Added `"bart"` as an option for `distance`. This estimates propensity scores using Bayesian Additive Regression Trees (BART) as implemented in the `dbarts` package.

* In methods that accept it, `m.order` can be set to "`data`", which matches in the order the data appear. With `distance = "mahalanobis"`, `m.order` can be "`random`" or "`data`", with "`data`" as the default. Otherwise, `m.order` can be `"largest"`, `"smallest"`, `"random"`, or `"data"`, with `"largest"` as the default (consistent with prior behavior).

* `match.data()`, which is used to create matched datasets, has a few new arguments. The `data` argument can be supplied with a dataset that will have the matching weights and subclasses added. If not supplied, `match.data()` will try to figure out the appropriate dataset like it did in the past. The `drop.unmatched` aregument controls whether unmatched units are dropped from the output. The default is `TRUE`, consistent with prior behavior. Warnings are now more informative.

* The included dataset, `lalonde`, now uses a `race` variable instead of separate `black` and `hispan` variables. This makes it easier to see how character variables are treated by `MatchIt` functions.

* Bugs in `distance = "rpart"` have been fixed.

* When arguments are supplied to methods that don't accept them, a warning will be thrown.

## `method = "nearest"`

* With `method = "nearest"`, a `subclass` component is now included in the output when `replace = FALSE` (the default), as it has been with optimal and full matching.

* Performance improvements for nearest neighbor matching.

* When using `method = "nearest"` with `distance = "mahalanobis"`, factor variables can now be included in the main `formula`. The design matrix no longer has to be full rank because a generalized inverse is used to compute the Mahalanobis distance.

* When using Mahalanobis distance matching, unless `m.order = "random"`, results will be identical across runs. Previously, several random choices would occur to break ties. Ties are broken based on the order of the data; shuffling the order of the data may therefore yield different matches.

* When using `method = "nearest"` with a caliper specified, the nearest control unit will be matched to the treated unit if one is available. Previously, a random control unit within the caliper would be selected. This eliminates the need for the `calclosest` argument, which has been removed.

## `method = "optimal"` and `method = "full"`

* Fixed bug in `method = "optimal"`, which produced results that did not match `optmatch`. Now they do.

* Added support for optimal and full Mahalanobis distance matching by setting `method = "mahalanobis"` with `method = "optimal"` and `method = "full"`. Previously, both methods would perform a random match if `method` was set to `"mahalanobis"`. Now they use the native support in `optmatch::pairmatch()` and `optmatch::fullmatch()` for Mahalanobis distance matching.

* Added support for exact matching with `method = "optimal"` and `method = "full"`. As with `method = "nearest"`, the names of the variables for which exact matches are required should be supplied to the `exact` argument. This relies on `optmatch::exactMatch()`.

* The warning that used to occur when using `method = "optimal"` and `method = "full"` about the order of the match not guaranteed to be the same as the original data no longer occurs.

* The `match.matrix` component of the `matchit()` output object is now produced with `method = "full"`. Note that the matrix can be much wider than with other methods because all control units that are not discarded are matched.

## `method = "genetic"`

* Fixed a bug with `method = "genetic"` that caused an error with some `ratio` greater than 1.

* The default of `replace` in `method = "genetic"` is now `FALSE`, as it is with `method = "nearest"`.

* When `verbose = FALSE`, the default, no output is printed with `method = "genetic"`. With `verbose = TRUE`, the printed output of `Matching::GenMatch()` with `print.level = 2` is displayed.

* The `exact` argument now correctly functions with `method = "genetic"`. Previously, it would have to be specified in accordance with its use in `Matching::GenMatch()`.

* Different ways to match on variables are now allowed with `method = "genetic"`, similar to how they are with `method = "nearest"`. If `distance = "mahalanobis"`, no distance measure will be computed, and genetic matching will be performed just on the variables supplied to `formula`. Otherwise, if `mahvars` are not specified, genetic matching will be performed on the variables supplied to `formula` and the distance measure. If `mahvars` are specified, genetic matching will be performed on the variables supplied to `mahvars` and the distance measure, but balance will be optimized on all covariates supplied to `formula`. Previously, `mahvars` was ignored. Balance is now always optimized on the variables included in `formula` and never on the distance measure, whereas in the past the distance measure was always included in the balance optimization.

* The `caliper` argument now works as it does with `method = "nearest"`, so that it applies only to the distance measure and not the other matching covariates, as it did in the past. It is ignored with `distance = "mahalanobis"`.

* With `method = "genetic"`, a `subclass` component is now included in the output when `replace = FALSE` (the default), as it has been with optimal and full matching.

## `method = "cem"`

* With `method = "cem"`, the `k2k` argument is now recognized. Previously it was ignored unless an argument to `k2k.method` was supplied.

## `method = "subclass"`

* Performance improvements.

* When any estimated subclass doesn't have any members from a treatment group, units from other subclasses are pulled to fill it so that every subclass will have at least one unit from each treatment group. This uses the same mechanism as is used in `WeightIt`.

* Rather than producing warnings and using the default number of subclasses (6), when an inappropriate argument is supplied to `subclass`, an error will occur.

* The new `subclass` argument to `summary()` can be used to control whether subclass balance statistics are computed; it can be `TRUE` (display balance for all subclasses), `FALSE` (display balance for no subclasses), or a vector of subclass indices on which to assess balance. The default is `FALSE`.

* With `summary()`, balance aggregating across subclasses is now computed using subclass weights instead of by combining the subclass-specific balance statistics.

## `summary.matchit()`

* In `summary.matchit()`, when `interactions = TRUE`, interactions are no longer computed with the distance measure or between dummies variables of the same factor. Variable names are cleaned up and easier to read.

* In `summary.matchit()`, the argument to `addlvariables` can be specified as a data frame or matrix of covariates, a formula with the additional covariates (and transformations) on the right side, or a character vector containing the names of the additional covariates. For the latter two, if the variables named do not exist in the `X` component of the `matchit` output object or in the environment, an argument to `data` can be supplied to `summary()` that contains these variables.

* The output for `summary()` is now the same for all methods. Previously there were different methods for a few different types of matching. In addition, the QQ and eCDF statistics have been adjusted. Both now use the weights that were computed as part of the matching. The QQ and eCDF statistics for binary variables are set to the difference in group proportions. The standard deviation of the control group has been removed from the output.

## `plot.matchit()`

* Plots now use weighted summaries when weights are present, removing the need for the `num.draws` argument.

* The appearance of some plots has improved (e.g., text is appropriately centered, axes are more clearly labelled). For QQ plots with binary variables or variables that take on only a few values, the plots look more like clusters than snakes.

* The argument to `type` can be abbreviated (e.g., `"j"` for jitter).

* Fixed a bug that caused all plots generated after using `plot(., type = "hist")` to be small.

* When specifying an argument to `which.xs` to control for which variables balance is displayed graphically, the input should be the name of the original variable rather than the version that appears in the `summary()` output. In particular, if a factor variable was supplied to `matchit()`, it should be referred to by its name rather than the names of its split dummies. This makes it easier to view balance on factor variables without having to know or type the names of all their levels.

* QQ plots can now be used with all matching methods. Previously, attempting `plot()` after `method = "exact"` would fail.

## `plot.matchit.summary()`

* The summary plot has been completley redesigned. It is now a Love plot made using `graphics::dotchart()`. A few options are available for ordering the variables, presenting absolute or raw standardized mean differences, and placing threshold lines on the plots. For a more sophisticated interface, see `cobalt::love.plot()`, which natively supports `matchit` objects and uses `ggplot2` as its engine.