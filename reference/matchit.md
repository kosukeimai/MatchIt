# Matching for Causal Inference

`matchit()` is the main function of *MatchIt* and performs pairing,
subset selection, and subclassification with the aim of creating
treatment and control groups balanced on included covariates. *MatchIt*
implements the suggestions of Ho, Imai, King, and Stuart (2007) for
improving parametric statistical models by preprocessing data with
nonparametric matching methods.

This page documents the overall use of `matchit()`, but for specifics of
how `matchit()` works with individual matching methods, see the
individual pages linked in the Details section below.

## Usage

``` r
matchit(
  formula,
  data = NULL,
  method = "nearest",
  distance = "glm",
  link = "logit",
  distance.options = list(),
  estimand = "ATT",
  exact = NULL,
  mahvars = NULL,
  antiexact = NULL,
  discard = "none",
  reestimate = FALSE,
  s.weights = NULL,
  replace = FALSE,
  m.order = NULL,
  caliper = NULL,
  std.caliper = TRUE,
  ratio = 1,
  verbose = FALSE,
  include.obj = FALSE,
  normalize = TRUE,
  ...
)
```

## Arguments

- formula:

  a two-sided [`formula`](https://rdrr.io/r/stats/formula.html) object
  containing the treatment and covariates to be used in creating the
  distance measure used in the matching. This formula will be supplied
  to the functions that estimate the distance measure. The formula
  should be specified as `A ~ X1 + X2 + ...` where `A` represents the
  treatment variable and `X1` and `X2` are covariates.

- data:

  a data frame containing the variables named in `formula` and possible
  other arguments. If not found in `data`, the variables will be sought
  in the environment.

- method:

  the matching method to be used. The allowed methods are
  [`"nearest"`](https://kosukeimai.github.io/MatchIt/reference/method_nearest.md)
  for nearest neighbor matching (on the propensity score by default),
  [`"optimal"`](https://kosukeimai.github.io/MatchIt/reference/method_optimal.md)
  for optimal pair matching,
  [`"full"`](https://kosukeimai.github.io/MatchIt/reference/method_full.md)
  for optimal full matching,
  [`"quick"`](https://kosukeimai.github.io/MatchIt/reference/method_quick.md)
  for generalized (quick) full matching,
  [`"genetic"`](https://kosukeimai.github.io/MatchIt/reference/method_genetic.md)
  for genetic matching,
  [`"cem"`](https://kosukeimai.github.io/MatchIt/reference/method_cem.md)
  for coarsened exact matching,
  [`"exact"`](https://kosukeimai.github.io/MatchIt/reference/method_exact.md)
  for exact matching,
  [`"cardinality"`](https://kosukeimai.github.io/MatchIt/reference/method_cardinality.md)
  for cardinality and profile matching, and
  [`"subclass"`](https://kosukeimai.github.io/MatchIt/reference/method_subclass.md)
  for subclassification. When set to `NULL`, no matching will occur, but
  propensity score estimation and common support restrictions will still
  occur if requested. See the linked pages for each method for more
  details on what these methods do, how the arguments below are used by
  each on, and what additional arguments are allowed.

- distance:

  the distance measure to be used. Can be either the name of a method of
  estimating propensity scores (e.g., `"glm"`), the name of a method of
  computing a distance matrix from the covariates (e.g.,
  `"mahalanobis"`), a vector of already-computed distance measures, or a
  matrix of pairwise distances. See
  [`distance`](https://kosukeimai.github.io/MatchIt/reference/distance.md)
  for allowable options. The default is `"glm"` for propensity scores
  estimated with logistic regression using
  [`glm()`](https://rdrr.io/r/stats/glm.html). Ignored for some methods;
  see individual methods pages for information on whether and how the
  distance measure is used.

- link:

  when `distance` is specified as a string, an additional argument
  controlling the link function used in estimating the distance measure.
  Allowable options depend on the specific `distance` value specified.
  See
  [`distance`](https://kosukeimai.github.io/MatchIt/reference/distance.md)
  for allowable options with each option. The default is `"logit"`,
  which, along with `distance = "glm"`, identifies the default measure
  as logistic regression propensity scores.

- distance.options:

  a named list containing additional arguments supplied to the function
  that estimates the distance measure as determined by the argument to
  `distance`. See
  [`distance`](https://kosukeimai.github.io/MatchIt/reference/distance.md)
  for an example of its use.

- estimand:

  a string containing the name of the target estimand desired. Can be
  one of `"ATT"`, `"ATC"`, or `"ATE"`. Default is `"ATT"`. See Details
  and the individual methods pages for information on how this argument
  is used.

- exact:

  for methods that allow it, for which variables exact matching should
  take place. Can be specified as a string containing the names of
  variables in `data` to be used or a one-sided formula with the desired
  variables on the right-hand side (e.g., `~ X3 + X4`). See the
  individual methods pages for information on whether and how this
  argument is used.

- mahvars:

  for methods that allow it, on which variables Mahalanobis distance
  matching should take place when `distance` corresponds to propensity
  scores. Usually used to perform Mahalanobis distance matching within
  propensity score calipers, where the propensity scores are computed
  using `formula` and `distance`. Can be specified as a string
  containing the names of variables in `data` to be used or a one-sided
  formula with the desired variables on the right-hand side (e.g.,
  `~ X3 + X4`). See the individual methods pages for information on
  whether and how this argument is used.

- antiexact:

  for methods that allow it, for which variables anti-exact matching
  should take place. Anti-exact matching ensures paired individuals do
  not have the same value of the anti-exact matching variable(s). Can be
  specified as a string containing the names of variables in `data` to
  be used or a one-sided formula with the desired variables on the
  right-hand side (e.g., `~ X3 + X4`). See the individual methods pages
  for information on whether and how this argument is used.

- discard:

  a string containing a method for discarding units outside a region of
  common support. When a propensity score is estimated or supplied to
  `distance` as a vector, the options are `"none"`, `"treated"`,
  `"control"`, or `"both"`. For `"none"`, no units are discarded for
  common support. Otherwise, units whose propensity scores fall outside
  the corresponding region are discarded. Can also be a `logical` vector
  where `TRUE` indicates the unit is to be discarded. Default is
  `"none"` for no common support restriction. See Details.

- reestimate:

  if `discard` is not `"none"` and propensity scores are estimated,
  whether to re-estimate the propensity scores in the remaining sample.
  Default is `FALSE` to use the propensity scores estimated in the
  original sample.

- s.weights:

  an optional numeric vector of sampling weights to be incorporated into
  propensity score models and balance statistics. Can also be specified
  as a string containing the name of variable in `data` to be used or a
  one-sided formula with the variable on the right-hand side (e.g.,
  `~ SW`). Not all propensity score models accept sampling weights; see
  [`distance`](https://kosukeimai.github.io/MatchIt/reference/distance.md)
  for information on which do and do not, and see
  [`vignette("sampling-weights")`](https://kosukeimai.github.io/MatchIt/articles/sampling-weights.md)
  for details on how to use sampling weights in a matching analysis.

- replace:

  for methods that allow it, whether matching should be done with
  replacement (`TRUE`), where control units are allowed to be matched to
  several treated units, or without replacement (`FALSE`), where control
  units can only be matched to one treated unit each. See the individual
  methods pages for information on whether and how this argument is
  used. Default is `FALSE` for matching without replacement.

- m.order:

  for methods that allow it, the order that the matching takes place.
  Allowable options depend on the matching method. The default of `NULL`
  corresponds to `"largest"` when a propensity score is estimated or
  supplied as a vector and `"data"` otherwise.

- caliper:

  for methods that allow it, the width(s) of the caliper(s) to use in
  matching. Should be a numeric vector with each value named according
  to the variable to which the caliper applies. To apply to the distance
  measure, the value should be unnamed. See the individual methods pages
  for information on whether and how this argument is used. Positive
  values require the distance between paired units to be no larger than
  the supplied caliper; negative values require the distance between
  paired units to be larger than the absolute value value of the
  supplied caliper. The default is `NULL` for no caliper.

- std.caliper:

  `logical`; when a caliper is specified, whether the the caliper is in
  standard deviation units (`TRUE`) or raw units (`FALSE`). Can either
  be of length 1, applying to all calipers, or of length equal to the
  length of `caliper`. Default is `TRUE`.

- ratio:

  for methods that allow it, how many control units should be matched to
  each treated unit in k:1 matching. Should be a single integer value.
  See the individual methods pages for information on whether and how
  this argument is used. The default is 1 for 1:1 matching.

- verbose:

  `logical`; whether information about the matching process should be
  printed to the console. What is printed depends on the matching
  method. Default is `FALSE` for no printing other than warnings.

- include.obj:

  `logical`; whether to include any objects created in the matching
  process in the output, i.e., by the functions from other packages
  `matchit()` calls. What is included depends on the matching method.
  Default is `FALSE`.

- normalize:

  `logical`; whether to rescale the nonzero weights in each treatment
  group to have an average of 1. Default is `TRUE`. See "How Matching
  Weights Are Computed" below for more details.

- ...:

  additional arguments passed to the functions used in the matching
  process. See the individual methods pages for information on what
  additional arguments are allowed for each method.

## Value

When `method` is something other than `"subclass"`, a `matchit` object
with the following components:

- match.matrix:

  a matrix containing the matches. The row names correspond to the
  treated units and the values in each row are the names (or indices) of
  the control units matched to each treated unit. When treated units are
  matched to different numbers of control units (e.g., with variable
  ratio matching or matching with a caliper), empty spaces will be
  filled with `NA`. Not included when `method` is `"full"`, `"cem"`
  (unless `k2k = TRUE`), `"exact"`, `"quick"`, or `"cardinality"`
  (unless `mahvars` is supplied and `ratio` is an integer).

- subclass:

  a factor containing matching pair/stratum membership for each unit.
  Unmatched units will have a value of `NA`. Not included when
  `replace = TRUE` or when `method = "cardinality"` unless `mahvars` is
  supplied and `ratio` is an integer.

- weights:

  a numeric vector of estimated matching weights. Unmatched and
  discarded units will have a weight of zero.

- model:

  the fit object of the model used to estimate propensity scores when
  `distance` is specified as a method of estimating propensity scores.
  When `reestimate = TRUE`, this is the model estimated after discarding
  units.

- X:

  a data frame of covariates mentioned in `formula`, `exact`, `mahvars`,
  `caliper`, and `antiexact`.

- call:

  the `matchit()` call.

- info:

  information on the matching method and distance measures used.

- estimand:

  the argument supplied to `estimand`.

- formula:

  the `formula` supplied.

- treat:

  a vector of treatment status converted to zeros (0) and ones (1) if
  not already in that format.

- distance:

  a vector of distance values (i.e., propensity scores) when `distance`
  is supplied as a method of estimating propensity scores or a numeric
  vector.

- discarded:

  a logical vector denoting whether each observation was discarded
  (`TRUE`) or not (`FALSE`) by the argument to `discard`.

- s.weights:

  the vector of sampling weights supplied to the `s.weights` argument,
  if any.

- exact:

  a one-sided formula containing the variables, if any, supplied to
  `exact`.

- mahvars:

  a one-sided formula containing the variables, if any, supplied to
  `mahvars`.

- obj:

  when `include.obj = TRUE`, an object containing the intermediate
  results of the matching procedure. See the individual methods pages
  for what this component will contain.

When `method = "subclass"`, a `matchit.subclass` object with the same
components as above except that `match.matrix` is excluded and one
additional component, `q.cut`, is included, containing a vector of the
distance measure cutpoints used to define the subclasses. See
[`method_subclass`](https://kosukeimai.github.io/MatchIt/reference/method_subclass.md)
for details.

## Details

Details for the various matching methods can be found at the following
help pages:

- [`method_nearest`](https://kosukeimai.github.io/MatchIt/reference/method_nearest.md)
  for nearest neighbor matching

- [`method_optimal`](https://kosukeimai.github.io/MatchIt/reference/method_optimal.md)
  for optimal pair matching

- [`method_full`](https://kosukeimai.github.io/MatchIt/reference/method_full.md)
  for optimal full matching

- [`method_quick`](https://kosukeimai.github.io/MatchIt/reference/method_quick.md)
  for generalized (quick) full matching

- [`method_genetic`](https://kosukeimai.github.io/MatchIt/reference/method_genetic.md)
  for genetic matching

- [`method_cem`](https://kosukeimai.github.io/MatchIt/reference/method_cem.md)
  for coarsened exact matching

- [`method_exact`](https://kosukeimai.github.io/MatchIt/reference/method_exact.md)
  for exact matching

- [`method_cardinality`](https://kosukeimai.github.io/MatchIt/reference/method_cardinality.md)
  for cardinality and profile matching

- [`method_subclass`](https://kosukeimai.github.io/MatchIt/reference/method_subclass.md)
  for subclassification

The pages contain information on what the method does, which of the
arguments above are allowed with them and how they are interpreted, and
what additional arguments can be supplied to further tune the method.
Note that the default method with no arguments supplied other than
`formula` and `data` is 1:1 nearest neighbor matching without
replacement on a propensity score estimated using a logistic regression
of the treatment on the covariates. This is not the same default offered
by other matching programs, such as those in *Matching*, `teffects` in
Stata, or `PROC PSMATCH` in SAS, so care should be taken if trying to
replicate the results of those programs.

When `method = NULL`, no matching will occur, but any propensity score
estimation and common support restriction will. This can be a simple way
to estimate the propensity score for use in future matching
specifications without having to re-estimate it each time. The
`matchit()` output with no matching can be supplied to
[`summary()`](https://rdrr.io/r/base/summary.html) to examine balance
prior to matching on any of the included covariates and on the
propensity score if specified. All arguments other than `distance`,
`discard`, and `reestimate` will be ignored.

See
[`distance`](https://kosukeimai.github.io/MatchIt/reference/distance.md)
for details on the several ways to specify the `distance`, `link`, and
`distance.options` arguments to estimate propensity scores and create
distance measures.

When the treatment variable is not a `0/1` variable, it will be coerced
to one and returned as such in the `matchit()` output (see section
Value, below). The following rules are used: 1) if `0` is one of the
values, it will be considered the control and the other value the
treated; 2) otherwise, if the variable is a factor, `levels(treat)[1]`
will be considered control and the other value the treated; 3)
otherwise, `sort(unique(treat))[1]` will be considered control and the
other value the treated. It is safest to ensure the treatment variable
is a `0/1` variable.

The `discard` option implements a common support restriction. It can
only be used when a distance measure is an estimated propensity score or
supplied as a vector and is ignored for some matching methods. When
specified as `"treated"`, treated units whose distance measure is
outside the range of distance measures of the control units will be
discarded. When specified as `"control"`, control units whose distance
measure is outside the range of distance measures of the treated units
will be discarded. When specified as `"both"`, treated and control units
whose distance measure is outside the intersection of the range of
distance measures of the treated units and the range of distance
measures of the control units will be discarded. When
`reestimate = TRUE` and `distance` corresponds to a propensity
score-estimating function, the propensity scores are re-estimated in the
remaining units prior to being used for matching or calipers.

Caution should be used when interpreting effects estimated with various
values of `estimand`. Setting `estimand = "ATT"` doesn't necessarily
mean the average treatment effect in the treated is being estimated; it
just means that for matching methods, treated units will be untouched
and given weights of 1 and control units will be matched to them (and
the opposite for `estimand = "ATC"`). If a caliper is supplied or
treated units are removed for common support or some other reason (e.g.,
lacking matches when using exact matching), the actual estimand targeted
is not the ATT but the treatment effect in the matched sample. The
argument to `estimand` simply triggers which units are matched to which,
and for stratification-based methods (exact matching, CEM, full
matching, and subclassification), determines the formula used to compute
the stratification weights.

### How Matching Weights Are Computed

Matching weights are computed in one of two ways depending on whether
matching was done with replacement or not.

#### Matching without replacement and subclassification

For matching *without* replacement (except for cardinality matching),
including subclassification, each unit is assigned to a subclass, which
represents the pair they are a part of (in the case of k:1 matching) or
the stratum they belong to (in the case of exact matching, coarsened
exact matching, full matching, or subclassification). The formula for
computing the weights depends on the argument supplied to `estimand`. A
new "stratum propensity score" (\\p^s_i\\) is computed for each unit
\\i\\ as \\p^s_i = \frac{1}{n_s}\sum\_{j: s_j =s_i}{I(A_j=1)}\\ where
\\n_s\\ is the size of subclass \\s\\ and \\I(A_j=1)\\ is 1 if unit
\\j\\ is treated and 0 otherwise. That is, the stratum propensity score
for stratum \\s\\ is the proportion of units in stratum \\s\\ that are
in the treated group, and all units in stratum \\s\\ are assigned that
stratum propensity score. This is distinct from the propensity score
used for matching, if any. Weights are then computed using the standard
formulas for inverse probability weights with the stratum propensity
score inserted:

- for the ATT, weights are 1 for the treated units and
  \\\frac{p^s}{1-p^s}\\ for the control units

- for the ATC, weights are \\\frac{1-p^s}{p^s}\\ for the treated units
  and 1 for the control units

- for the ATE, weights are \\\frac{1}{p^s}\\ for the treated units and
  \\\frac{1}{1-p^s}\\ for the control units.

For cardinality matching, all matched units receive a weight of 1.

#### Matching with replacement

For matching *with* replacement, units are not assigned to unique
strata. For the ATT, each treated unit gets a weight of 1. Each control
unit is weighted as the sum of the inverse of the number of control
units matched to the same treated unit across its matches. For example,
if a control unit was matched to a treated unit that had two other
control units matched to it, and that same control was matched to a
treated unit that had one other control unit matched to it, the control
unit in question would get a weight of \\1/3 + 1/2 = 5/6\\. For the ATC,
the same is true with the treated and control labels switched. The
weights are computed using the `match.matrix` component of the
`matchit()` output object.

#### Normalized weights

When `normalize = TRUE` (the default), in each treatment group, weights
are divided by the mean of the nonzero weights in that treatment group
to make the weights sum to the number of units in that treatment group
(i.e., to have an average of 1).

#### Sampling weights

If sampling weights are included through the `s.weights` argument, they
will be included in the `matchit()` output object but not incorporated
into the matching weights.
[`match_data()`](https://kosukeimai.github.io/MatchIt/reference/match_data.md),
which extracts the matched set from a `matchit` object, combines the
matching weights and sampling weights.

## References

Ho, D. E., Imai, K., King, G., & Stuart, E. A. (2007). Matching as
Nonparametric Preprocessing for Reducing Model Dependence in Parametric
Causal Inference. *Political Analysis*, 15(3), 199–236.
[doi:10.1093/pan/mpl013](https://doi.org/10.1093/pan/mpl013)

Ho, D. E., Imai, K., King, G., & Stuart, E. A. (2011). MatchIt:
Nonparametric Preprocessing for Parametric Causal Inference. *Journal of
Statistical Software*, 42(8).
[doi:10.18637/jss.v042.i08](https://doi.org/10.18637/jss.v042.i08)

## See also

[`summary.matchit()`](https://kosukeimai.github.io/MatchIt/reference/summary.matchit.md)
for balance assessment after matching,
[`plot.matchit()`](https://kosukeimai.github.io/MatchIt/reference/plot.matchit.md)
for plots of covariate balance and propensity score overlap after
matching.

- [`vignette("MatchIt")`](https://kosukeimai.github.io/MatchIt/articles/MatchIt.md)
  for an introduction to matching with *MatchIt*

- [`vignette("matching-methods")`](https://kosukeimai.github.io/MatchIt/articles/matching-methods.md)
  for descriptions of the variety of matching methods and options
  available

- [`vignette("assessing-balance")`](https://kosukeimai.github.io/MatchIt/articles/assessing-balance.md)
  for information on assessing the quality of a matching specification

- [`vignette("estimating-effects")`](https://kosukeimai.github.io/MatchIt/articles/estimating-effects.md)
  for instructions on how to estimate treatment effects after matching

- [`vignette("sampling-weights")`](https://kosukeimai.github.io/MatchIt/articles/sampling-weights.md)
  for a guide to using *MatchIt* with sampling weights.

## Author

Daniel Ho, Kosuke Imai, Gary King, and Elizabeth Stuart wrote the
original package. Starting with version 4.0.0, Noah Greifer is the
primary maintainer and developer.

## Examples

``` r
data("lalonde")

# Default: 1:1 NN PS matching w/o replacement

m.out1 <- matchit(treat ~ age + educ + race + nodegree +
                    married + re74 + re75,
                  data = lalonde)
m.out1
#> A `matchit` object
#>  - method: 1:1 nearest neighbor matching without replacement
#>  - distance: Propensity score
#>              - estimated with logistic regression
#>  - number of obs.: 614 (original), 370 (matched)
#>  - target estimand: ATT
#>  - covariates: age, educ, race, nodegree, married, re74, re75
summary(m.out1)
#> 
#> Call:
#> matchit(formula = treat ~ age + educ + race + nodegree + married + 
#>     re74 + re75, data = lalonde)
#> 
#> Summary of Balance for All Data:
#>            Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean
#> distance          0.5774        0.1822          1.7941     0.9211    0.3774
#> age              25.8162       28.0303         -0.3094     0.4400    0.0813
#> educ             10.3459       10.2354          0.0550     0.4959    0.0347
#> raceblack         0.8432        0.2028          1.7615          .    0.6404
#> racehispan        0.0595        0.1422         -0.3498          .    0.0827
#> racewhite         0.0973        0.6550         -1.8819          .    0.5577
#> nodegree          0.7081        0.5967          0.2450          .    0.1114
#> married           0.1892        0.5128         -0.8263          .    0.3236
#> re74           2095.5737     5619.2365         -0.7211     0.5181    0.2248
#> re75           1532.0553     2466.4844         -0.2903     0.9563    0.1342
#>            eCDF Max
#> distance     0.6444
#> age          0.1577
#> educ         0.1114
#> raceblack    0.6404
#> racehispan   0.0827
#> racewhite    0.5577
#> nodegree     0.1114
#> married      0.3236
#> re74         0.4470
#> re75         0.2876
#> 
#> Summary of Balance for Matched Data:
#>            Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean
#> distance          0.5774        0.3629          0.9739     0.7566    0.1321
#> age              25.8162       25.3027          0.0718     0.4568    0.0847
#> educ             10.3459       10.6054         -0.1290     0.5721    0.0239
#> raceblack         0.8432        0.4703          1.0259          .    0.3730
#> racehispan        0.0595        0.2162         -0.6629          .    0.1568
#> racewhite         0.0973        0.3135         -0.7296          .    0.2162
#> nodegree          0.7081        0.6378          0.1546          .    0.0703
#> married           0.1892        0.2108         -0.0552          .    0.0216
#> re74           2095.5737     2342.1076         -0.0505     1.3289    0.0469
#> re75           1532.0553     1614.7451         -0.0257     1.4956    0.0452
#>            eCDF Max Std. Pair Dist.
#> distance     0.4216          0.9740
#> age          0.2541          1.3938
#> educ         0.0757          1.2474
#> raceblack    0.3730          1.0259
#> racehispan   0.1568          1.0743
#> racewhite    0.2162          0.8390
#> nodegree     0.0703          1.0106
#> married      0.0216          0.8281
#> re74         0.2757          0.7965
#> re75         0.2054          0.7381
#> 
#> Sample Sizes:
#>           Control Treated
#> All           429     185
#> Matched       185     185
#> Unmatched     244       0
#> Discarded       0       0
#> 

# 1:1 NN Mahalanobis distance matching w/ replacement and
# exact matching on married and race

m.out2 <- matchit(treat ~ age + educ + race + nodegree +
                    married + re74 + re75,
                  data = lalonde,
                  distance = "mahalanobis",
                  replace = TRUE,
                  exact = ~ married + race)
m.out2
#> A `matchit` object
#>  - method: 1:1 nearest neighbor matching with replacement
#>  - distance: Mahalanobis - number of obs.: 614 (original), 265 (matched)
#>  - target estimand: ATT
#>  - covariates: age, educ, race, nodegree, married, re74, re75
summary(m.out2, un = TRUE)
#> 
#> Call:
#> matchit(formula = treat ~ age + educ + race + nodegree + married + 
#>     re74 + re75, data = lalonde, distance = "mahalanobis", exact = ~married + 
#>     race, replace = TRUE)
#> 
#> Summary of Balance for All Data:
#>            Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean
#> age              25.8162       28.0303         -0.3094     0.4400    0.0813
#> educ             10.3459       10.2354          0.0550     0.4959    0.0347
#> raceblack         0.8432        0.2028          1.7615          .    0.6404
#> racehispan        0.0595        0.1422         -0.3498          .    0.0827
#> racewhite         0.0973        0.6550         -1.8819          .    0.5577
#> nodegree          0.7081        0.5967          0.2450          .    0.1114
#> married           0.1892        0.5128         -0.8263          .    0.3236
#> re74           2095.5737     5619.2365         -0.7211     0.5181    0.2248
#> re75           1532.0553     2466.4844         -0.2903     0.9563    0.1342
#>            eCDF Max
#> age          0.1577
#> educ         0.1114
#> raceblack    0.6404
#> racehispan   0.0827
#> racewhite    0.5577
#> nodegree     0.1114
#> married      0.3236
#> re74         0.4470
#> re75         0.2876
#> 
#> Summary of Balance for Matched Data:
#>            Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean
#> age              25.8162       25.6162          0.0280     0.6513    0.0466
#> educ             10.3459       10.3946         -0.0242     1.1564    0.0065
#> raceblack         0.8432        0.8432          0.0000          .    0.0000
#> racehispan        0.0595        0.0595          0.0000          .    0.0000
#> racewhite         0.0973        0.0973          0.0000          .    0.0000
#> nodegree          0.7081        0.7135         -0.0119          .    0.0054
#> married           0.1892        0.1892          0.0000          .    0.0000
#> re74           2095.5737     1861.6424          0.0479     1.4978    0.0286
#> re75           1532.0553     1091.6516          0.1368     2.0335    0.0347
#>            eCDF Max Std. Pair Dist.
#> age          0.1784          0.4918
#> educ         0.0324          0.2070
#> raceblack    0.0000          0.0000
#> racehispan   0.0000          0.0000
#> racewhite    0.0000          0.0000
#> nodegree     0.0054          0.0119
#> married      0.0000          0.0000
#> re74         0.1784          0.2606
#> re75         0.0811          0.2445
#> 
#> Sample Sizes:
#>               Control Treated
#> All            429.       185
#> Matched (ESS)   34.89     185
#> Matched         80.       185
#> Unmatched      349.         0
#> Discarded        0.         0
#> 

# 2:1 NN Mahalanobis distance matching within caliper defined
# by a probit pregression PS

m.out3 <- matchit(treat ~ age + educ + race + nodegree +
                    married + re74 + re75,
                  data = lalonde,
                  distance = "glm",
                  link = "probit",
                  mahvars = ~ age + educ + re74 + re75,
                  caliper = .1,
                  ratio = 2)
m.out3
#> A `matchit` object
#>  - method: 2:1 nearest neighbor matching without replacement
#>  - distance: Mahalanobis [matching]
#>              Propensity score [caliper]
#> 
#>              - estimated with probit regression
#>  - caliper: <distance> (0.029)
#>  - number of obs.: 614 (original), 257 (matched)
#>  - target estimand: ATT
#>  - covariates: age, educ, race, nodegree, married, re74, re75
summary(m.out3, un = TRUE)
#> 
#> Call:
#> matchit(formula = treat ~ age + educ + race + nodegree + married + 
#>     re74 + re75, data = lalonde, distance = "glm", link = "probit", 
#>     mahvars = ~age + educ + re74 + re75, caliper = 0.1, ratio = 2)
#> 
#> Summary of Balance for All Data:
#>            Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean
#> distance          0.5773        0.1817          1.8276     0.8777    0.3774
#> age              25.8162       28.0303         -0.3094     0.4400    0.0813
#> educ             10.3459       10.2354          0.0550     0.4959    0.0347
#> raceblack         0.8432        0.2028          1.7615          .    0.6404
#> racehispan        0.0595        0.1422         -0.3498          .    0.0827
#> racewhite         0.0973        0.6550         -1.8819          .    0.5577
#> nodegree          0.7081        0.5967          0.2450          .    0.1114
#> married           0.1892        0.5128         -0.8263          .    0.3236
#> re74           2095.5737     5619.2365         -0.7211     0.5181    0.2248
#> re75           1532.0553     2466.4844         -0.2903     0.9563    0.1342
#>            eCDF Max
#> distance     0.6413
#> age          0.1577
#> educ         0.1114
#> raceblack    0.6404
#> racehispan   0.0827
#> racewhite    0.5577
#> nodegree     0.1114
#> married      0.3236
#> re74         0.4470
#> re75         0.2876
#> 
#> Summary of Balance for Matched Data:
#>            Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean
#> distance          0.5113        0.4932          0.0835     1.0778    0.0253
#> age              26.0721       24.9459          0.1574     0.4310    0.0928
#> educ             10.4144       10.3514          0.0314     0.6279    0.0171
#> raceblack         0.7387        0.7252          0.0372          .    0.0135
#> racehispan        0.0991        0.0946          0.0190          .    0.0045
#> racewhite         0.1622        0.1802         -0.0608          .    0.0180
#> nodegree          0.6667        0.6396          0.0594          .    0.0270
#> married           0.1892        0.2297         -0.1035          .    0.0405
#> re74           3016.7936     2280.7013          0.1506     1.8738    0.0569
#> re75           2023.1731     1525.9838          0.1544     2.0215    0.0434
#>            eCDF Max Std. Pair Dist.
#> distance     0.1441          0.0860
#> age          0.3198          0.9487
#> educ         0.0586          0.7324
#> raceblack    0.0135          0.0565
#> racehispan   0.0045          0.4924
#> racewhite    0.0180          0.3236
#> nodegree     0.0270          0.5725
#> married      0.0405          0.4722
#> re74         0.2117          0.5502
#> re75         0.1081          0.5885
#> 
#> Sample Sizes:
#>               Control Treated
#> All            429.       185
#> Matched (ESS)  131.78     111
#> Matched        146.       111
#> Unmatched      283.        74
#> Discarded        0.         0
#> 

# Optimal full PS matching for the ATE within calipers on
# PS, age, and educ
m.out4 <- matchit(treat ~ age + educ + race + nodegree +
                    married + re74 + re75,
                  data = lalonde,
                  method = "full",
                  estimand = "ATE",
                  caliper = c(.1, age = 2, educ = 1),
                  std.caliper = c(TRUE, FALSE, FALSE))
m.out4
#> A `matchit` object
#>  - method: Optimal full matching
#>  - distance: Propensity score [caliper]
#> 
#>              - estimated with logistic regression
#>  - caliper: <distance> (0.029), age (2), educ (1)
#>  - number of obs.: 614 (original), 314 (matched)
#>  - target estimand: ATE
#>  - covariates: age, educ, race, nodegree, married, re74, re75
summary(m.out4, un = TRUE)
#> 
#> Call:
#> matchit(formula = treat ~ age + educ + race + nodegree + married + 
#>     re74 + re75, data = lalonde, method = "full", estimand = "ATE", 
#>     caliper = c(0.1, age = 2, educ = 1), std.caliper = c(TRUE, 
#>         FALSE, FALSE))
#> 
#> Summary of Balance for All Data:
#>            Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean
#> distance          0.5774        0.1822          1.7569     0.9211    0.3774
#> age              25.8162       28.0303         -0.2419     0.4400    0.0813
#> educ             10.3459       10.2354          0.0448     0.4959    0.0347
#> raceblack         0.8432        0.2028          1.6708          .    0.6404
#> racehispan        0.0595        0.1422         -0.2774          .    0.0827
#> racewhite         0.0973        0.6550         -1.4080          .    0.5577
#> nodegree          0.7081        0.5967          0.2355          .    0.1114
#> married           0.1892        0.5128         -0.7208          .    0.3236
#> re74           2095.5737     5619.2365         -0.5958     0.5181    0.2248
#> re75           1532.0553     2466.4844         -0.2870     0.9563    0.1342
#>            eCDF Max
#> distance     0.6444
#> age          0.1577
#> educ         0.1114
#> raceblack    0.6404
#> racehispan   0.0827
#> racewhite    0.5577
#> nodegree     0.1114
#> married      0.3236
#> re74         0.4470
#> re75         0.2876
#> 
#> Summary of Balance for Matched Data:
#>            Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean
#> distance          0.3515        0.3473          0.0185     1.0175    0.0146
#> age              22.4064       22.0133          0.0430     0.8438    0.0156
#> educ             10.7009       10.6440          0.0230     1.0814    0.0107
#> raceblack         0.4618        0.4586          0.0083          .    0.0032
#> racehispan        0.1497        0.1083          0.1387          .    0.0414
#> racewhite         0.3885        0.4331         -0.1125          .    0.0446
#> nodegree          0.5792        0.5824         -0.0068          .    0.0032
#> married           0.2894        0.2677          0.0482          .    0.0217
#> re74           2679.9337     2952.3545         -0.0461     1.0852    0.0401
#> re75           1430.9478     1745.5573         -0.0966     1.2043    0.0535
#>            eCDF Max Std. Pair Dist.
#> distance     0.0625          0.0449
#> age          0.1078          0.1260
#> educ         0.0594          0.2067
#> raceblack    0.0032          0.0325
#> racehispan   0.0414          0.5148
#> racewhite    0.0446          0.3771
#> nodegree     0.0032          0.2018
#> married      0.0217          0.5360
#> re74         0.1923          0.5470
#> re75         0.1493          0.6744
#> 
#> Sample Sizes:
#>               Control Treated
#> All             429.   185.  
#> Matched (ESS)   133.1   39.37
#> Matched         203.   111.  
#> Unmatched       226.    74.  
#> Discarded         0.     0.  
#> 
# Subclassification on a logistic PS with 10 subclasses after
# discarding controls outside common support of PS

s.out1 <- matchit(treat ~ age + educ + race + nodegree +
                    married + re74 + re75,
                  data = lalonde,
                  method = "subclass",
                  distance = "glm",
                  discard = "control",
                  subclass = 10)
s.out1
#> A `matchit` object
#>  - method: Subclassification (10 subclasses)
#>  - distance: Propensity score [common support]
#> 
#>              - estimated with logistic regression
#>  - common support: control units dropped
#>  - number of obs.: 614 (original), 557 (matched)
#>  - target estimand: ATT
#>  - covariates: age, educ, race, nodegree, married, re74, re75
summary(s.out1, un = TRUE)
#> 
#> Call:
#> matchit(formula = treat ~ age + educ + race + nodegree + married + 
#>     re74 + re75, data = lalonde, method = "subclass", distance = "glm", 
#>     discard = "control", subclass = 10)
#> 
#> Summary of Balance for All Data:
#>            Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean
#> distance          0.5774        0.1822          1.7941     0.9211    0.3774
#> age              25.8162       28.0303         -0.3094     0.4400    0.0813
#> educ             10.3459       10.2354          0.0550     0.4959    0.0347
#> raceblack         0.8432        0.2028          1.7615          .    0.6404
#> racehispan        0.0595        0.1422         -0.3498          .    0.0827
#> racewhite         0.0973        0.6550         -1.8819          .    0.5577
#> nodegree          0.7081        0.5967          0.2450          .    0.1114
#> married           0.1892        0.5128         -0.8263          .    0.3236
#> re74           2095.5737     5619.2365         -0.7211     0.5181    0.2248
#> re75           1532.0553     2466.4844         -0.2903     0.9563    0.1342
#>            eCDF Max
#> distance     0.6444
#> age          0.1577
#> educ         0.1114
#> raceblack    0.6404
#> racehispan   0.0827
#> racewhite    0.5577
#> nodegree     0.1114
#> married      0.3236
#> re74         0.4470
#> re75         0.2876
#> 
#> Summary of Balance Across Subclasses
#>            Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean
#> distance          0.5774        0.5710          0.0293     0.9338    0.0158
#> age              25.8162       25.3714          0.0622     0.4577    0.0866
#> educ             10.3459       10.4094         -0.0316     0.6894    0.0150
#> raceblack         0.8432        0.8262          0.0469          .    0.0170
#> racehispan        0.0595        0.0676         -0.0343          .    0.0081
#> racewhite         0.0973        0.1062         -0.0302          .    0.0089
#> nodegree          0.7081        0.6782          0.0658          .    0.0299
#> married           0.1892        0.1785          0.0274          .    0.0107
#> re74           2095.5737     2232.5096         -0.0280     1.3102    0.0449
#> re75           1532.0553     1643.4179         -0.0346     1.4216    0.0472
#>            eCDF Max
#> distance     0.0541
#> age          0.3043
#> educ         0.0425
#> raceblack    0.0170
#> racehispan   0.0081
#> racewhite    0.0089
#> nodegree     0.0299
#> married      0.0107
#> re74         0.2731
#> re75         0.1841
#> 
#> Sample Sizes:
#>               Control Treated
#> All            429.       185
#> Matched (ESS)   72.59     185
#> Matched        372.       185
#> Unmatched        0.         0
#> Discarded       57.         0
#> 
```
