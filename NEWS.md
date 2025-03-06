---
output:
  html_document: default
  pdf_document: default
---
`MatchIt` News and Updates
======

# MatchIt 4.7.1

* Updates for CRAN compatibility.

* Typo fixes in documentation and vignettes, thanks to @iagogv3.

# MatchIt 4.7.0

* For nearest neighbor matching, optimal full matching, and genetic matching, calipers can now be negative, which forces paired units to be further away from each other on the given variables.

* When using `method = "optimal"` to do 1:1 matching on a propensity score, a new preprocessing algorithm is used to speed up the matching process. This algorithm ensures the resulting match is as good as a match that would have been found without the preprocessing step while shrinking the size of the matching problem. Typically, the same matched set will be selected, but units may be paired differently (this is because the general optimization problem often has multiple solutions with the same value of the objective function; this algorithm adds an additional constraint to select among fewer such solutions). This algorithm is described in [Sävje (2020)](https://doi.org/10.1214/19-STS739) and is implemented in C++ to be fast.

* Fixed a bug when matching with a nonzero `ratio` where subclass membership was incorrectly calculated. Thanks to Simon Loewe (@simon-lowe) for originally pointing it out. (#207, #208)

* `match.data()` has been renamed to `match_data()`, but `match.data()` will remain as an alias for backward compatibility.

* Fixed a bug with printing.

* Documentation fixes.

# MatchIt 4.6.0

Most improvements are related to performance. Some of these dramatically improve speeds for large datasets. Most come from improvements to `Rcpp` code.

* When using `method = "nearest"`, `m.order` can now be set to `"farthest"` to prioritize hard-to-match treated units. Note this **does not** implement "far matching" but simply changes the order in which the closest matches are selected.

* Speed improvements to `method = "nearest"`, especially when matching on a propensity score.

* Speed improvements to `summary()` when `pair.dist = TRUE` and a `match.matrix` component is not included in the output (e.g., for `method = "full"` or `method = "quick"`).

* Speed improvements to `method = "subclass"` with `min.n` greater than 0.

* A new `normalize` argument has been added to `matchit()`. When set to `TRUE` (the default, which used to be the only option), the nonzero weights in each treatment group are rescaled to have an average of 1. When `FALSE`, the weights generated directly by the matching are returned instead.

* When using `method = "nearest"` with `m.order = "closest"`, the full distance matrix is no longer computed, which increases support for larger samples. This uses an adaptation of an algorithm described by [Rassen et al. (2012)](https://doi.org/10.1002/pds.3263).

* When using `method = "nearest"` with `verbose = TRUE`, the progress bar now displays an estimate of how much time remains.

* When using `method = "nearest"` with `m.order = "closest"` and `ratio` greater than 1, all eligible units will receive their first match before any receive their second, etc. Previously, the closest pairs would be matched regardless of whether other units had been matched. This ensures consistency with other `m.order` arguments.

* Speed and memory improvements to `method = "cem"` with many covariates and a large sample size. Previous versions used a Cartesian expansion of all levels of factor variables, which could easily explode.

* When using `method = "cem"` with `k2k = TRUE`, `m.order` can be set to select the matching order. Allowable options include `"data"` (the default), `"closest"`, `"farthest"`, and `"random"`. `"closest"` is recommended, but `"data"` is the default for now to remain consistent with previous versions.

* Documentation updates.

* Fixed a bug when using `method = "optimal"` or `method = "full"` with `discard` specified and `data` given as a tibble (`tbl_df` object). (#185)

* Fixed a bug when using `method = "cardinality"` with a single covariate. (#194)

# MatchIt 4.5.5

* When using `method = "cardinality"`, a new solver, HiGHS, can be requested by setting `solver = "highs"`, which relies on the `highs` package. This is much faster and more reliable than GLPK and is free and easy to install as a regular R package with no additional requirements.

* Fixed a bug when using `method = "optimal"` with `discard` and `exact` specified. Thanks to @NikNakk for the issue and fix. (#171)

# MatchIt 4.5.4

* With `method = "nearest"`, `m.order` can now be set to `"closest"` to request that the closest potential pairs are matched first. This can be used whether a propensity score is used or not.

* Fixed bugs when `distance = NULL` and no covariates are specified in `matchit()`.

* Changed "empirical cumulative density function" to "empirical cumulative distribution function" in documentation. (#166)

* Fixed a bug where calipers would not work properly on some systems. Thanks to Bill Dunlap for the solution. (#163)

* Fixed a bug when `.` was present in formulas. Thanks to @dmolitor. (#167)

* Fixed a bug when nearest neighbor matching for the ATC with `distance` supplied as a numeric distance matrix.

# MatchIt 4.5.3

* Error messages have been improved using `chk` and `rlang`, which are now dependencies.

* Fixed a bug when using `method = "nearest"` with `replace = TRUE` and `ratio` greater than 1. Thanks to Julia Kretschmann. (#159)

* Fixed a bug when using `method = "nearest"` with `exact` and `ratio` greater than 1. Thanks to Sarah Conner.

* Fixed a bug that would occur due to numerical imprecision in `plot.matchit()`. Thanks to @hkmztrk. (#158)

* Fixed bugs when using `method = "cem"` where a covariate was to be omitted from coarsening. Thanks to @jfhelmer. (#160)

* Fixed some typos in the vignettes. Thanks to @fBedecarrats. (#156)

* Updated vignettes to use `marginaleffects` v0.11.0 syntax.

# MatchIt 4.5.2

* Fixed a bug when using `method = "quick"` with `exact` specified. Thanks to @m-marquis. (#149)

* Improved performance and fixed some bugs when using `exact` in cases where some strata contain units from only one treatment group. Thanks to @m-marquis and others for pointing these out. (#151)

# MatchIt 4.5.1

* Nearest neighbor matching now uses a much faster algorithm (up to 6x times faster) when `distance` is a propensity score and `mahvars` is not specified. Differences in sort order might cause results to differ from previous versions if there are units with identical propensity scores.

* Template matching has been renamed profile matching in all documentation.

* After cardinality or profile matching using `method = "cardinality"` with `ratio` set to a whole number, it is possible to perform optimal Mahalanobis distance matching in the matched sample by supplying the desired matching variables to `mahvars`. Previously, the user had to run a separate pairing step.

* Fixed some typos in the vignettes.

* Fixed a bug where character variables would be flagged as non-finite. Thanks to @isfraser. (#138)

* Added alt text to images in README and vignettes. (#134)

# MatchIt 4.5.0

* Generalized full matching, as described by [Sävje, Higgins, and Sekhon (2021)](https://doi.org/10.1017/pan.2020.32), can now be implemented by setting `method = "quick"` in `matchit()`. It is a dramatically faster alternative to optimal full matching that can support much larger datasets and otherwise has similar balancing performance. See `?method_quick` and `vignette("matching-methods")` for more information. This functionality relies on the `quickmatch` package.

* The package structure has been updated, include with the use of Roxygen for documentation. This should not affect use, but the source code will look different from that of previous versions.

* When `method = "subclass"` and `min.n = 0` (which is not the default), any units not placed into a subclass are now considered "unmatched" and given weights of 0. Previously they were left in.

* When `method = "genetic"`, the default `distance.tolerance` is now 0. In previous versions, this argument was ignored; now it is not.

* For `plot.matchit()`, the `which.xs` argument can be specified as a one-sided formula. A new `data` argument is allowed if the variables in that formula are not among the original covariates.

* When a factor variable is supplied to `plot.matchit()` with `type = "density"`, the plot now displays all factor levels in the same plot instead of in separate plots for each level, similar to `cobalt::bal.plot()`.

* The "Estimating Effects" vignette (`vignette("estimating-effects")`) has been rewritten to be much shorter (and hopefully clearer) and to use the `marginaleffects` package, which is now a Suggested package. The new vignette focuses on using g-computation to estimate treatment effects using a single workflow with slight modifications for different situations.

* The error message when covariates have missing or non-finite values is now clearer, identifying which variables are afflicted. This fixes a bug mentioned in #115.

* Fixed a bug when using `matchit()` with `method = "cem"`, `k2k = TRUE`, and `k2k.method = NULL`. Thanks to Florian B. Mayr.

* Fixed a bug when using `method = "optimal"` and `method = "full"` with `exact` and `antiexact` specified, wherein a warning would occur about the `drop` argument in subsetting.

* Fixed a bug where `antiexact` would not work correctly with `method = "nearest"`. Thanks to @gli-1. (#119)

* Fixed typos in the documentation and vignettes.

* Calculating pair distances in `summary()` with `pair.dist = TRUE` is now faster.

* Improved printing of balance results when no covariates are supplied.

* Updates to the Estimating Effects vignette that dramatically increase the speed of the cluster bootstrap for average marginal effects after matching. Thanks to Yohei Hashimoto for pointing out the inefficiency.

* Updates to the Assessing Balance vignette to fix errors

* All vignettes and help files are better protected against Suggested packages not available on CRAN.

# MatchIt 4.4.0

* `optmatch` has returned to CRAN, now with an open-source license! A new `solver` argument can be passed to `matchit()` with `method = "full"` and `method = "optimal"` to control the solver used to perform the optimization used in the matching. Note that using the default (open source) solver LEMON may yield results different from those obtained prior to `optmatch` 0.10.0. For reproducibility questions, please contact the `optmatch` maintainers.

* New functions have been added to compute the Euclidean distance (`euclidean_dist()`), scaled Euclidean distance (`scaled_euclidean_dist()`), Mahalanobis distance (`mahalanobis_dist()`), and robust Mahalanobis distance (`robust_mahalanobis_dist()`). They produce distance matrices that can be supplied to the `distance` argument of `matchit()`, but see below.

* New distance options are available for `matchit()` based on the distance functions above: `"robust_mahalanobis"`, `"euclidean"`, and `"scaled_euclidean"`, which complement `"mahalanobis"`. Similar to `"mahalanobis"`, these do not involve estimating a propensity score but rather operate on the covariates directly. These can be used for nearest neighbor matching, optimal matching, full matching, and coarsened exact matching with `k2k = TRUE`.

* The Mahalanobis distance is now computed using the pooled within-group covariance matrix (computed by treatment group-mean centering each covariate before computing the covariance in the full sample), in line with how it is computed in `optmatch` and recommended by Rubin (1980) among others. This will cause results to differ between this version and prior versions of `MatchIt` that used the Mahalanobis distance computed ignoring group membership.

* Added the `unit.id` argument to `matchit()` with `method = "nearest"`, which defines unit IDs so that if a control observation with a given unit ID has been matched to a treated unit, no other control units with the same ID can be used as future matches, ensuring each unit ID is used no more than once. This is useful when, e.g., multiple rows correspond to the same control firm but you only want each control firm to be matched once, in which case firm ID would be supplied to `unit.id`. See [here](https://stats.stackexchange.com/questions/349784/propensity-matching-and-analysis-of-resultant-data-on-a-data-set-with-repeated-m) for an example use case.

* In `summary.matchit()`, `improvement` is now set to `FALSE` by default to hide the percentage improvement in balance. Set to `TRUE` to recover prior behavior.

* Added clearer errors when required packages are missing for certain `distance` methods.

* Fixed a bug when using `matchit()` with `method = "nearest"`, `ratio` greater than 1, and `reuse.max` specified. The bug allowed a previously matched control unit to be matched to the same treatment unit, thereby essentially ignoring the `ratio` argument. It now works as intended.

* Fixed a bug in `matchit()` with `method = "nearest"` when `distance` was supplied as a matrix and `Inf` values were present.

* Fixed a bug when using exact matching that caused an infinite loop when variable levels contained commas. Thanks to @bking124. (#111)

* Fixed a bug introduced by `optmatch` version 0.10.3.

* Documentation updates.

* Updated the logo, thanks to [Ben Stillerman](https://stillben.com).

# MatchIt 4.3.4

* `optmatch` has been removed from CRAN. Instructions on installing it are in `?method_optimal` and `?method_full`.

* When `s.weights` are supplied with `distance = "randomforest"`, the weights are supplied to `randomForest::randomForest()`.

* Improved conditional use of packages, especially `optmatch`. This may mean that certain examples fail to run in the vignettes.

# MatchIt 4.3.3

* Fixed a bug where `rbind.matchdata()` would produce datasets twice their expected length. Thanks to @sconti555. (#98)

# MatchIt 4.3.2

* Fixed a bug where the `q.cut` component of the `matchit` object when `method = "subclass"` was not included. Now it is. Thanks to @aldencabajar. (#92)

* The `nn` and `qn` components of the `matchit` object have been removed. They are now computed by `summary.matchit()` and included in the `summary.matchit` object.

* Removed the code to disable compiler checks to satisfy CRAN requirements.

# MatchIt 4.3.1

* Added the `reuse.max` argument to `matchit()` with `method = "nearest"`. This controls the maximum number of times each control unit can be used as a match. Setting `reuse.max = 1` is equivalent to matching without replacement (i.e., like setting `replace = FALSE`), and setting `reuse.max = Inf` is equivalent to matching with replacement with no restriction on the reuse of controls (i.e., like setting `replace = TRUE`). Values in between restrict how many times each control unit can be used as a match. Higher values will tend to improve balance but decrease precision.

* Mahalanobis distance matching with `method = "nearest"` is now a bit faster. 

* Fixed a bug where `method = "full"` would fail when some exact matching strata contained exactly one treated unit and exactly one control unit. (#88)

* Fixed a bug introduced in 4.3.0 where the inclusion of character variables would cause the error `"Non-finite values are not allowed in the covariates."` Thanks to Moaath Mustafa.

* Documentation updates.

# MatchIt 4.3.0

* Cardinality and template matching can now be used by setting `method = "cardinality"` in `matchit()`. These methods use mixed integer programming to directly select a matched subsample without pairing or stratifying units that satisfied user-supplied balance constraints. Their results can be dramatically improved when using the Gurobi optimizer. See `?method_cardinality` and `vignette("matching-methods")` for more information.

* Added `"lasso"`, `"ridge"`, and `"elasticnet"` as options for `distance`. These estimate propensity scores using lasso, ridge, or elastic net regression, respectively, as implemented in the `glmnet` package.

* Added `"gbm"` as an option for `distance`. This estimates propensity scores using generalized boosted models as implemented in the `gbm` package. This implementation differs from that in `twang` by using cross-validation or out-of-bag error to choose the tuning parameter as opposed to balance.

* A new argument, `include.obj`, has been added to `matchit()`. When `TRUE`, the intermediate matching object created internally will be included in the output in the `obj` component. See the individual methods pages for information on what is included in each output.  This is ignored for some methods.

* Density plots can now be requested using `plot.matchit()` by setting `type = "density"`. These display the density of each covariate in the treatment groups before and after matching and are similar to the plots created by `cobalt::bal.plot()`. Density plots can be easier to interpret than eCDF plots. `vignette("assessing-balance")` has been updated with this addition.

* A clearer error is now produced when the treatment variable is omitted from the `formula` argument to `matchit()`.

* Improvements in how `match.data()` finds the original dataset. It's still always safer to supply an argument to `data`, but now `match.data()` will look in the environment of the `matchit` formula, then the calling environment of `match.data()`, then the `model` component of the `matchit` object. A clearer error message is now printed when a valid dataset cannot be found in these places.

* Fixed a bug that would occur when using `summary.matchit()` with just one covariate.

* When `verbose = TRUE` and a propensity score is estimated (i.e., using the `distance` argument), a message saying so will be displayed.

* Fixed a bug in `print.matchit()` where it would indicate that the propensity score was used in a caliper if any caliper was specified, even if not on the propensity score. Now, it will only indicate that the propensity score was used in a caliper if it actually was.

* Fixed a bug in `plot.matchit()` that would occur when a level of a factor had no values.

* Speed improvements for `method = "full"` with `exact` specified. These changes can make current results differ slightly from past results when the `tol` value is high. It is recommended to always use a low value of `tol`.

* Typo fixes in documentation and vignettes.

* Fixed a bug where supplying a "GAM" string to the `distance` argument (i.e., using the syntax prior to version 4.0.0) would ignore the link supplied.

* When an incompatible argument is supplied to `matchit()` (e.g., `reestimate` with `distance = "mahalanobis"`), an error or warning will only be produced when that argument has been set to a value other than its default (e.g., so setting `reestimate = FALSE` will no longer throw an error). This fixes an issue brought up by Vu Ng when using `MatchThem`.

* A clearer error is produced when non-finite values are present in the covariates.

# MatchIt 4.2.0

* `distance` can now be supplied as a distance matrix containing pairwise distances with nearest neighbor, optimal, and full matching. This means users can create a distance matrix outside `MatchIt` (e.g., using `optmatch::match_on()` or `dist()`) and `matchit()` will use those distances in the matching. See `?distance` for details.

* Added `rbind.matchdata()` method for `matchdata` and `getmatches` objects (the output of `match.data()` and `get_matches()`, respectively) to avoid subclass conflicts when combining matched samples after matching within subgroups.

* Added a section in `vignette("estimating-effects")` on moderation analysis with matching, making use of the new `rbind()` method.

* Added `antiexact` argument to perform anti-exact matching, i.e., matching that ensures treated and control units have different values of certain variables. See [here](https://stackoverflow.com/questions/66526115/propensity-score-matching-with-panel-data) and [here](https://stackoverflow.com/questions/61120201/avoiding-duplicates-from-propensity-score-matching?rq=1) for examples where this feature was requested and might be useful. Anti-exact matching works with nearest neighbor, optimal, full, and genetic matching. The argument to `antiexact` should be similar to an argument to `exact`: either a string or a one-sided `formula` containing the names of the anti-exact matching variables.

* Slight speed improvements for nearest neighbor matching, especially with `exact` specified.

* With `method = "nearest"`, `verbose = TRUE`, and `exact` specified, separate messages and progress bars will be shown for each subgroup of the `exact` variable(s).

* A spurious warning that would appear when using a large `ratio` with `replace = TRUE` and `method = "nearest"` no longer appears.

* Fixed a bug when trying to supply `distance` as a labeled numeric vector (e.g., resulting from `haven`).

* Fixed some typos in the documentation and vignettes.

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
