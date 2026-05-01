# Optimal Pair Matching

In
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md),
setting `method = "optimal"` performs optimal pair matching. The
matching is optimal in the sense that that sum of the absolute pairwise
distances in the matched sample is as small as possible. The method
functionally relies on
[`optmatch::fullmatch()`](https://rdrr.io/pkg/optmatch/man/fullmatch.html).

Advantages of optimal pair matching include that the matching order is
not required to be specified and it is less likely that extreme
within-pair distances will be large, unlike with nearest neighbor
matching. Generally, however, as a subset selection method, optimal pair
matching tends to perform similarly to nearest neighbor matching in that
similar subsets of units will be selected to be matched.

This page details the allowable arguments with `method = "optmatch"`.
See
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
for an explanation of what each argument means in a general context and
how it can be specified.

Below is how
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
is used for optimal pair matching:


    matchit(formula,
            data = NULL,
            method = "optimal",
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
            ratio = 1,
            min.controls = NULL,
            max.controls = NULL,
            verbose = FALSE,
            ...) 

## Arguments

- formula:

  a two-sided [formula](https://rdrr.io/r/stats/formula.html) object
  containing the treatment and covariates to be used in creating the
  distance measure used in the matching. This formula will be supplied
  to the functions that estimate the distance measure.

- data:

  a data frame containing the variables named in `formula`. If not found
  in `data`, the variables will be sought in the environment.

- method:

  set here to `"optimal"`.

- distance:

  the distance measure to be used. See
  [`distance`](https://kosukeimai.github.io/MatchIt/reference/distance.md)
  for allowable options. Can be supplied as a distance matrix.

- link:

  when `distance` is specified as a method of estimating propensity
  scores, an additional argument controlling the link function used in
  estimating the distance measure. See
  [`distance`](https://kosukeimai.github.io/MatchIt/reference/distance.md)
  for allowable options with each option.

- distance.options:

  a named list containing additional arguments supplied to the function
  that estimates the distance measure as determined by the argument to
  `distance`.

- estimand:

  a string containing the desired estimand. Allowable options include
  `"ATT"` and `"ATC"`. See Details.

- exact:

  for which variables exact matching should take place.

- mahvars:

  for which variables Mahalanobis distance matching should take place
  when `distance` corresponds to a propensity score (e.g., for caliper
  matching or to discard units for common support). If specified, the
  distance measure will not be used in matching.

- antiexact:

  for which variables anti-exact matching should take place. Anti-exact
  matching is processed using
  [`optmatch::antiExactMatch()`](https://rdrr.io/pkg/optmatch/man/antiExactMatch.html).

- discard:

  a string containing a method for discarding units outside a region of
  common support. Only allowed when `distance` is not `"mahalanobis"`
  and not a matrix.

- reestimate:

  if `discard` is not `"none"`, whether to re-estimate the propensity
  score in the remaining sample prior to matching.

- s.weights:

  the variable containing sampling weights to be incorporated into
  propensity score models and balance statistics.

- ratio:

  how many control units should be matched to each treated unit for k:1
  matching. For variable ratio matching, see section "Variable Ratio
  Matching" in Details below.

- min.controls, max.controls:

  for variable ratio matching, the minimum and maximum number of
  controls units to be matched to each treated unit. See section
  "Variable Ratio Matching" in Details below.

- verbose:

  `logical`; whether information about the matching process should be
  printed to the console. What is printed depends on the matching
  method. Default is `FALSE` for no printing other than warnings.

- ...:

  additional arguments passed to
  [`optmatch::fullmatch()`](https://rdrr.io/pkg/optmatch/man/fullmatch.html).
  Allowed arguments include `tol` and `solver`. See the
  [`optmatch::fullmatch()`](https://rdrr.io/pkg/optmatch/man/fullmatch.html)
  documentation for details. In general, `tol` should be set to a low
  number (e.g., `1e-7`) to get a more precise solution (default is
  `1e-3`).

  The arguments `replace`, `caliper`, and `m.order` are ignored with a
  warning.

## Details

### Mahalanobis Distance Matching

Mahalanobis distance matching can be done one of two ways:

1.  If no propensity score needs to be estimated, `distance` should be
    set to `"mahalanobis"`, and Mahalanobis distance matching will occur
    using all the variables in `formula`. Arguments to `discard` and
    `mahvars` will be ignored. For example, to perform simple
    Mahalanobis distance matching, the following could be run:


        matchit(treat ~ X1 + X2, method = "nearest",
                distance = "mahalanobis") 

    With this code, the Mahalanobis distance is computed using `X1` and
    `X2`, and matching occurs on this distance. The `distance` component
    of the
    [`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
    output will be empty.

2.  If a propensity score needs to be estimated for common support with
    `discard`, `distance` should be whatever method is used to estimate
    the propensity score or a vector of distance measures, i.e., it
    should not be `"mahalanobis"`. Use `mahvars` to specify the
    variables used to create the Mahalanobis distance. For example, to
    perform Mahalanobis after discarding units outside the common
    support of the propensity score in both groups, the following could
    be run:


        matchit(treat ~ X1 + X2 + X3, method = "nearest",
                distance = "glm", discard = "both",
                mahvars = ~ X1 + X2) 

    With this code, `X1`, `X2`, and `X3` are used to estimate the
    propensity score (using the `"glm"` method, which by default is
    logistic regression), which is used to identify the common support.
    The actual matching occurs on the Mahalanobis distance computed only
    using `X1` and `X2`, which are supplied to `mahvars`. The estimated
    propensity scores will be included in the `distance` component of
    the
    [`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
    output.

### Estimand

The `estimand` argument controls whether control units are selected to
be matched with treated units (`estimand = "ATT"`) or treated units are
selected to be matched with control units (`estimand = "ATC"`). The
"focal" group (e.g., the treated units for the ATT) is typically made to
be the smaller treatment group, and a warning will be thrown if it is
not set that. Setting `estimand = "ATC"` is equivalent to swapping all
treated and control labels for the treatment variable. When
`estimand = "ATC"`, the `match.matrix` component of the output will have
the names of the control units as the rownames and be filled with the
names of the matched treated units (opposite to when
`estimand = "ATT"`). Note that the argument supplied to `estimand`
doesn't necessarily correspond to the estimand actually targeted; it is
merely a switch to trigger which treatment group is considered "focal".

### Variable Ratio Matching

[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
can perform variable ratio matching, which involves matching a different
number of control units to each treated unit. When `ratio > 1`, rather
than requiring all treated units to receive `ratio` matches, the
arguments to `max.controls` and `min.controls` can be specified to
control the maximum and minimum number of matches each treated unit can
have. `ratio` controls how many total control units will be matched:
`n1 * ratio` control units will be matched, where `n1` is the number of
treated units, yielding the same total number of matched controls as
fixed ratio matching does.

Variable ratio matching can be used with any `distance` specification.
`ratio` does not have to be an integer but must be greater than 1 and
less than `n0/n1`, where `n0` and `n1` are the number of control and
treated units, respectively. Setting `ratio = n0/n1` performs a
restricted form of full matching where all control units are matched. If
`min.controls` is not specified, it is set to 1 by default.
`min.controls` must be less than `ratio`, and `max.controls` must be
greater than `ratio`. See the Examples section of
[`method_nearest()`](https://kosukeimai.github.io/MatchIt/reference/method_nearest.md)
for an example of their use, which is the same as it is with optimal
matching.

## Note

Optimal pair matching is a restricted form of optimal full matching
where the number of treated units in each subclass is equal to 1,
whereas in unrestricted full matching, multiple treated units can be
assigned to the same subclass.
[`optmatch::pairmatch()`](https://rdrr.io/pkg/optmatch/man/pairmatch.html)
is simply a wrapper for
[`optmatch::fullmatch()`](https://rdrr.io/pkg/optmatch/man/fullmatch.html),
which performs optimal full matching and is the workhorse for
[`method_full`](https://kosukeimai.github.io/MatchIt/reference/method_full.md).
In the same way,
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
uses
[`optmatch::fullmatch()`](https://rdrr.io/pkg/optmatch/man/fullmatch.html)
under the hood, imposing the restrictions that make optimal full
matching function like optimal pair matching (which is simply to set
`min.controls >= 1` and to pass `ratio` to the `mean.controls`
argument). This distinction is not important for regular use but may be
of interest to those examining the source code.

The option `"optmatch_max_problem_size"` is automatically set to `Inf`
during the matching process, different from its default in *optmatch*.
This enables matching problems of any size to be run, but may also let
huge, infeasible problems get through and potentially take a long time
or crash R. See
[`optmatch::setMaxProblemSize()`](https://rdrr.io/pkg/optmatch/man/setMaxProblemSize.html)
for more details.

A preprocessing algorithm describe by Sävje (2020;
[doi:10.1214/19-STS739](https://doi.org/10.1214/19-STS739) ) is used to
improve the speed of the matching when 1:1 matching on a propensity
score. It does so by adding an additional constraint that guarantees a
solution as optimal as the solution that would have been found without
the constraint, and that constraint often dramatically reduces the size
of the matching problem at no cost. However, this may introduce
differences between the results obtained by *MatchIt* and by *optmatch*,
though such differences will shrink when smaller values of `tol` are
used.

## Outputs

All outputs described in
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
are returned with `method = "optimal"`. When `include.obj = TRUE` in the
call to
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md),
the output of the call to
[`optmatch::fullmatch()`](https://rdrr.io/pkg/optmatch/man/fullmatch.html)
will be included in the output. When `exact` is specified, this will be
a list of such objects, one for each stratum of the `exact` variables.

## References

In a manuscript, be sure to cite the following paper if using
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
with `method = "optimal"`:

Hansen, B. B., & Klopfer, S. O. (2006). Optimal Full Matching and
Related Designs via Network Flows. Journal of Computational and
Graphical Statistics, 15(3), 609–627.
[doi:10.1198/106186006X137047](https://doi.org/10.1198/106186006X137047)

For example, a sentence might read:

*Optimal pair matching was performed using the MatchIt package (Ho,
Imai, King, & Stuart, 2011) in R, which calls functions from the
optmatch package (Hansen & Klopfer, 2006).*

## See also

[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
for a detailed explanation of the inputs and outputs of a call to
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md).

[`optmatch::fullmatch()`](https://rdrr.io/pkg/optmatch/man/fullmatch.html),
which is the workhorse.

[`method_full`](https://kosukeimai.github.io/MatchIt/reference/method_full.md)
for optimal full matching, of which optimal pair matching is a special
case, and which relies on similar machinery.

## Examples

``` r
data("lalonde")

#1:1 optimal PS matching with exact matching on race
m.out1 <- matchit(treat ~ age + educ + race +
                    nodegree + married + re74 + re75,
                  data = lalonde,
                  method = "optimal",
                  exact = ~race)
#> Warning: Fewer control units than treated units in some `exact` strata; not all
#> treated units will get a match.
m.out1
#> A `matchit` object
#>  - method: 1:1 optimal pair matching
#>  - distance: Propensity score
#>              - estimated with logistic regression
#>  - number of obs.: 614 (original), 232 (matched)
#>  - target estimand: ATT
#>  - covariates: age, educ, race, nodegree, married, re74, re75
summary(m.out1)
#> 
#> Call:
#> matchit(formula = treat ~ age + educ + race + nodegree + married + 
#>     re74 + re75, data = lalonde, method = "optimal", exact = ~race)
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
#> distance          0.4972        0.4901          0.0325     1.0026    0.0056
#> age              25.7845       25.9569         -0.0241     0.4984    0.0711
#> educ             10.0948       10.3966         -0.1501     0.6260    0.0286
#> raceblack         0.7500        0.7500          0.0000          .    0.0000
#> racehispan        0.0948        0.0948          0.0000          .    0.0000
#> racewhite         0.1552        0.1552          0.0000          .    0.0000
#> nodegree          0.6810        0.6121          0.1517          .    0.0690
#> married           0.2414        0.2672         -0.0660          .    0.0259
#> re74           2862.5644     3046.0687         -0.0376     1.4317    0.0569
#> re75           1785.7677     1873.2021         -0.0272     1.6241    0.0506
#>            eCDF Max Std. Pair Dist.
#> distance     0.0776          0.0436
#> age          0.2069          1.2024
#> educ         0.0776          1.1705
#> raceblack    0.0000          0.0000
#> racehispan   0.0000          0.0000
#> racewhite    0.0000          0.0000
#> nodegree     0.0690          0.9102
#> married      0.0259          0.5503
#> re74         0.2759          0.7299
#> re75         0.1379          0.8456
#> 
#> Sample Sizes:
#>           Control Treated
#> All           429     185
#> Matched       116     116
#> Unmatched     313      69
#> Discarded       0       0
#> 

#2:1 optimal matching on the scaled Euclidean distance
m.out2 <- matchit(treat ~ age + educ + race +
                    nodegree + married + re74 + re75,
                  data = lalonde,
                  method = "optimal",
                  ratio = 2,
                  distance = "scaled_euclidean")
m.out2
#> A `matchit` object
#>  - method: 2:1 optimal pair matching
#>  - distance: Scaled Euclidean - number of obs.: 614 (original), 555 (matched)
#>  - target estimand: ATT
#>  - covariates: age, educ, race, nodegree, married, re74, re75
summary(m.out2, un = FALSE)
#> 
#> Call:
#> matchit(formula = treat ~ age + educ + race + nodegree + married + 
#>     re74 + re75, data = lalonde, method = "optimal", distance = "scaled_euclidean", 
#>     ratio = 2)
#> 
#> Summary of Balance for Matched Data:
#>            Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean
#> age              25.8162       26.4757         -0.0922     0.4947    0.0658
#> educ             10.3459       10.2730          0.0363     0.5810    0.0286
#> raceblack         0.8432        0.2351          1.6726          .    0.6081
#> racehispan        0.0595        0.1216         -0.2629          .    0.0622
#> racewhite         0.0973        0.6432         -1.8422          .    0.5459
#> nodegree          0.7081        0.6243          0.1843          .    0.0838
#> married           0.1892        0.4432         -0.6487          .    0.2541
#> re74           2095.5737     4120.1633         -0.4143     0.7972    0.1577
#> re75           1532.0553     2107.1613         -0.1786     1.0608    0.0935
#>            eCDF Max Std. Pair Dist.
#> age          0.1405          0.5477
#> educ         0.0865          0.3750
#> raceblack    0.6081          1.6726
#> racehispan   0.0622          0.2629
#> racewhite    0.5459          1.8422
#> nodegree     0.0838          0.1843
#> married      0.2541          0.6487
#> re74         0.4081          0.4806
#> re75         0.2595          0.3195
#> 
#> Sample Sizes:
#>           Control Treated
#> All           429     185
#> Matched       370     185
#> Unmatched      59       0
#> Discarded       0       0
#> 
```
