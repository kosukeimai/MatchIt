# Nearest Neighbor Matching

In
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md),
setting `method = "nearest"` performs greedy nearest neighbor matching.
A distance is computed between each treated unit and each control unit,
and, one by one, each treated unit is assigned a control unit as a
match. The matching is "greedy" in the sense that there is no action
taken to optimize an overall criterion; each match is selected without
considering the other matches that may occur subsequently.

This page details the allowable arguments with `method = "nearest"`. See
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
for an explanation of what each argument means in a general context and
how it can be specified.

Below is how
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
is used for nearest neighbor matching:

    matchit(formula,
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
            replace = TRUE,
            m.order = NULL,
            caliper = NULL,
            ratio = 1,
            min.controls = NULL,
            max.controls = NULL,
            verbose = FALSE,
            ...) 

## Arguments

- formula:

  a two-sided [formula](https://rdrr.io/r/stats/formula.html) object
  containing the treatment and covariates to be used in creating the
  distance measure used in the matching.

- data:

  a data frame containing the variables named in `formula`. If not found
  in `data`, the variables will be sought in the environment.

- method:

  set here to `"nearest"`.

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

  for which variables exact matching should take place; two units with
  different values of an exact matching variable will not be paired.

- mahvars:

  for which variables Mahalanobis distance matching should take place
  when `distance` corresponds to a propensity score (e.g., for caliper
  matching or to discard units for common support). If specified, the
  distance measure will not be used in matching.

- antiexact:

  for which variables anti-exact matching should take place; two units
  with the same value of an anti-exact matching variable will not be
  paired.

- discard:

  a string containing a method for discarding units outside a region of
  common support. Only allowed when `distance` corresponds to a
  propensity score.

- reestimate:

  if `discard` is not `"none"`, whether to re-estimate the propensity
  score in the remaining sample prior to matching.

- s.weights:

  the variable containing sampling weights to be incorporated into
  propensity score models and balance statistics.

- replace:

  whether matching should be done with replacement (i.e., whether
  control units can be used as matches multiple times). See also the
  `reuse.max` argument below. Default is `FALSE` for matching without
  replacement.

- m.order:

  the order that the matching takes place. Allowable options include
  `"largest"`, where matching takes place in descending order of
  distance measures; `"smallest"`, where matching takes place in
  ascending order of distance measures; `"closest"`, where matching
  takes place in ascending order of the smallest distance between units;
  `"farthest"`, where matching takes place in descending order of the
  smallest distance between units; `"random"`, where matching takes
  place in a random order; and `"data"` where matching takes place based
  on the order of units in the data. When `m.order = "random"`, results
  may differ across different runs of the same code unless a seed is set
  and specified with [`set.seed()`](https://rdrr.io/r/base/Random.html).
  The default of `NULL` corresponds to `"largest"` when a propensity
  score is estimated or supplied as a vector and `"data"` otherwise. See
  Details for more information.

- caliper:

  the width(s) of the caliper(s) used for caliper matching. Two units
  with a difference on a caliper variable larger than the caliper will
  not be paired. See Details and Examples.

- std.caliper:

  `logical`; when calipers are specified, whether they are in standard
  deviation units (`TRUE`) or raw units (`FALSE`).

- ratio:

  how many control units should be matched to each treated unit for k:1
  matching. For variable ratio matching, see section "Variable Ratio
  Matching" in Details below. When `ratio` is greater than 1, all
  treated units will be attempted to be matched with a control unit
  before any treated unit is matched with a second control unit, etc.
  This reduces the possibility that control units will be used up before
  some treated units receive any matches.

- min.controls, max.controls:

  for variable ratio matching, the minimum and maximum number of
  controls units to be matched to each treated unit. See section
  "Variable Ratio Matching" in Details below.

- verbose:

  `logical`; whether information about the matching process should be
  printed to the console. When `TRUE`, a progress bar implemented using
  *RcppProgress* will be displayed along with an estimate of the time
  remaining.

- ...:

  additional arguments that control the matching specification:

  `reuse.max`

  :   `numeric`; the maximum number of times each control can be used as
      a match. Setting `reuse.max = 1` corresponds to matching without
      replacement (i.e., `replace = FALSE`), and setting
      `reuse.max = Inf` corresponds to traditional matching with
      replacement (i.e., `replace = TRUE`) with no limit on the number
      of times each control unit can be matched. Other values restrict
      the number of times each control can be matched when matching with
      replacement. `replace` is ignored when `reuse.max` is specified.

  `unit.id`

  :   one or more variables containing a unit ID for each observation,
      i.e., in case multiple observations correspond to the same unit.
      Once a control observation has been matched, no other observation
      with the same unit ID can be used as matches. This ensures each
      control unit is used only once even if it has multiple
      observations associated with it. Omitting this argument is the
      same as giving each observation a unique ID.

## Details

### Mahalanobis Distance Matching

Mahalanobis distance matching can be done one of two ways:

1.  If no propensity score needs to be estimated, `distance` should be
    set to `"mahalanobis"`, and Mahalanobis distance matching will occur
    using all the variables in `formula`. Arguments to `discard` and
    `mahvars` will be ignored, and a caliper can only be placed on named
    variables. For example, to perform simple Mahalanobis distance
    matching, the following could be run:

        matchit(treat ~ X1 + X2, method = "nearest",
                distance = "mahalanobis") 

    With this code, the Mahalanobis distance is computed using `X1` and
    `X2`, and matching occurs on this distance. The `distance` component
    of the
    [`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
    output will be empty.

2.  If a propensity score needs to be estimated for any reason, e.g.,
    for common support with `discard` or for creating a caliper,
    `distance` should be whatever method is used to estimate the
    propensity score or a vector of distance measures. Use `mahvars` to
    specify the variables used to create the Mahalanobis distance. For
    example, to perform Mahalanobis within a propensity score caliper,
    the following could be run:

        matchit(treat ~ X1 + X2 + X3, method = "nearest",
                distance = "glm", caliper = .25,
                mahvars = ~ X1 + X2) 

    With this code, `X1`, `X2`, and `X3` are used to estimate the
    propensity score (using the `"glm"` method, which by default is
    logistic regression), which is used to create a matching caliper.
    The actual matching occurs on the Mahalanobis distance computed only
    using `X1` and `X2`, which are supplied to `mahvars`. Units whose
    propensity score difference is larger than the caliper will not be
    paired, and some treated units may therefore not receive a match.
    The estimated propensity scores will be included in the `distance`
    component of the
    [`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
    output. See Examples.

### Estimand

The `estimand` argument controls whether control units are selected to
be matched with treated units (`estimand = "ATT"`) or treated units are
selected to be matched with control units (`estimand = "ATC"`). The
"focal" group (e.g., the treated units for the ATT) is typically made to
be the smaller treatment group, and a warning will be thrown if it is
not set that way unless `replace = TRUE`. Setting `estimand = "ATC"` is
equivalent to swapping all treated and control labels for the treatment
variable. When `estimand = "ATC"`, the default `m.order` is
`"smallest"`, and the `match.matrix` component of the output will have
the names of the control units as the rownames and be filled with the
names of the matched treated units (opposite to when
`estimand = "ATT"`). Note that the argument supplied to `estimand`
doesn't necessarily correspond to the estimand actually targeted; it is
merely a switch to trigger which treatment group is considered "focal".

### Variable Ratio Matching

[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
can perform variable ratio "extremal" matching as described by Ming and
Rosenbaum (2000;
[doi:10.1111/j.0006-341X.2000.00118.x](https://doi.org/10.1111/j.0006-341X.2000.00118.x)
). This method tends to result in better balance than fixed ratio
matching at the expense of some precision. When `ratio > 1`, rather than
requiring all treated units to receive `ratio` matches, each treated
unit is assigned a value that corresponds to the number of control units
they will be matched to. These values are controlled by the arguments
`min.controls` and `max.controls`, which correspond to \\\alpha\\ and
\\\beta\\, respectively, in Ming and Rosenbaum (2000), and trigger
variable ratio matching to occur. Some treated units will receive
`min.controls` matches and others will receive `max.controls` matches
(and one unit may have an intermediate number of matches); how many
units are assigned each number of matches is determined by the algorithm
described in Ming and Rosenbaum (2000, p119). `ratio` controls how many
total control units will be matched: `n1 * ratio` control units will be
matched, where `n1` is the number of treated units, yielding the same
total number of matched controls as fixed ratio matching does.

Variable ratio matching cannot be used with Mahalanobis distance
matching or when `distance` is supplied as a matrix. The calculations of
the numbers of control units each treated unit will be matched to occurs
without consideration of `caliper` or `discard`. `ratio` does not have
to be an integer but must be greater than 1 and less than `n0/n1`, where
`n0` and `n1` are the number of control and treated units, respectively.
Setting `ratio = n0/n1` performs a crude form of full matching where all
control units are matched. If `min.controls` is not specified, it is set
to 1 by default. `min.controls` must be less than `ratio`, and
`max.controls` must be greater than `ratio`. See Examples below for an
example of their use.

### Using `m.order = "closest"` or `"farthest"`

`m.order` can be set to `"closest"` or `"farthest"`, which work
regardless of how the distance measure is specified. This matches in
order of the distance between units. First, all the closest match is
found for all treated units and the pairwise distances computed; when
`m.order = "closest"` the pair with the smallest of the distances is
matched first, and when `m.order = "farthest"`, the pair with the
largest of the distances is matched first. Then, the pair with the
second smallest (or largest) is matched second. If the matched control
is ineligible (i.e., because it has already been used in a prior match),
a new match is found for the treated unit, the new pair's distance is
re-computed, and the pairs are re-ordered by distance.

Using `m.order = "closest"` ensures that the best possible matches are
given priority, and in that sense should perform similarly to
`m.order = "smallest"`. It can be used to ensure the best matches,
especially when matching with a caliper. Using `m.order = "farthest"`
ensures that the hardest units to match are given their best chance to
find a close match, and in that sense should perform similarly to
`m.order = "largest"`. It can be used to reduce the possibility of
extreme imbalance when there are hard-to-match units competing for
controls. Note that `m.order = "farthest"` **does not** implement "far
matching" (i.e., finding the farthest control unit from each treated
unit); it defines the order in which the closest matches are selected.

### Reproducibility

Nearest neighbor matching involves a random component only when
`m.order = "random"` (or when the propensity is estimated using a method
with randomness; see
[`distance`](https://kosukeimai.github.io/MatchIt/reference/distance.md)
for details), so a seed must be set in that case using
[`set.seed()`](https://rdrr.io/r/base/Random.html) to ensure
reproducibility. Otherwise, it is purely deterministic, and any ties are
broken based on the order in which the data appear.

## Outputs

All outputs described in
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
are returned with `method = "nearest"`. When `replace = TRUE`, the
`subclass` component is omitted. `include.obj` is ignored.

## References

In a manuscript, you don't need to cite another package when using
`method = "nearest"` because the matching is performed completely within
*MatchIt*. For example, a sentence might read:

*Nearest neighbor matching was performed using the MatchIt package (Ho,
Imai, King, & Stuart, 2011) in R.*

## See also

[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
for a detailed explanation of the inputs and outputs of a call to
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md).

[`method_optimal()`](https://kosukeimai.github.io/MatchIt/reference/method_optimal.md)
for optimal pair matching, which is similar to nearest neighbor matching
without replacement except that an overall distance criterion is
minimized (i.e., as an alternative to specifying `m.order`).

## Examples

``` r
data("lalonde")

# 1:1 greedy NN matching on the PS
m.out1 <- matchit(treat ~ age + educ + race + nodegree +
                    married + re74 + re75,
                  data = lalonde,
                  method = "nearest")
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
#>     re74 + re75, data = lalonde, method = "nearest")
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

# 3:1 NN Mahalanobis distance matching with
# replacement within a PS caliper
m.out2 <- matchit(treat ~ age + educ + race + nodegree +
                    married + re74 + re75,
                  data = lalonde,
                  method = "nearest",
                  replace = TRUE,
                  mahvars = ~ age + educ + re74 + re75,
                  ratio = 3,
                  caliper = .02)
m.out2
#> A `matchit` object
#>  - method: 3:1 nearest neighbor matching with replacement
#>  - distance: Mahalanobis [matching]
#>              Propensity score [caliper]
#> 
#>              - estimated with logistic regression
#>  - caliper: <distance> (0.006)
#>  - number of obs.: 614 (original), 300 (matched)
#>  - target estimand: ATT
#>  - covariates: age, educ, race, nodegree, married, re74, re75
summary(m.out2, un = FALSE)
#> 
#> Call:
#> matchit(formula = treat ~ age + educ + race + nodegree + married + 
#>     re74 + re75, data = lalonde, method = "nearest", mahvars = ~age + 
#>     educ + re74 + re75, replace = TRUE, caliper = 0.02, ratio = 3)
#> 
#> Summary of Balance for Matched Data:
#>            Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean
#> distance          0.5660        0.5660          0.0003     0.9880    0.0037
#> age              25.5000       24.0010          0.2095     0.5178    0.0835
#> educ             10.3086       10.3930         -0.0420     0.7356    0.0115
#> raceblack         0.8210        0.8169          0.0113          .    0.0041
#> racehispan        0.0679        0.0576          0.0435          .    0.0103
#> racewhite         0.1111        0.1255         -0.0486          .    0.0144
#> nodegree          0.7037        0.6852          0.0407          .    0.0185
#> married           0.1914        0.1132          0.1996          .    0.0782
#> re74           1831.9348     2119.9268         -0.0589     1.2783    0.0528
#> re75           1400.4354     1386.7669          0.0042     2.2735    0.0623
#>            eCDF Max Std. Pair Dist.
#> distance     0.0432          0.0128
#> age          0.3508          0.9787
#> educ         0.0381          0.7836
#> raceblack    0.0041          0.0358
#> racehispan   0.0103          0.2643
#> racewhite    0.0144          0.1845
#> nodegree     0.0185          0.5556
#> married      0.0782          0.4588
#> re74         0.2984          0.4634
#> re75         0.2428          0.5274
#> 
#> Sample Sizes:
#>               Control Treated
#> All            429.       185
#> Matched (ESS)   64.23     162
#> Matched        138.       162
#> Unmatched      291.        23
#> Discarded        0.         0
#> 

# 1:1 NN Mahalanobis distance matching within calipers
# on re74 and re75 and exact matching on married and race
m.out3 <- matchit(treat ~ age + educ + re74 + re75,
                  data = lalonde,
                  method = "nearest",
                  distance = "mahalanobis",
                  exact = ~ married + race,
                  caliper = c(re74 = .2, re75 = .15))
#> Warning: Fewer control units than treated units in some `exact` strata; not all
#> treated units will get a match.
m.out3
#> A `matchit` object
#>  - method: 1:1 nearest neighbor matching without replacement
#>  - distance: Mahalanobis - caliper: re74 (1295.593), re75 (494.352)
#>  - number of obs.: 614 (original), 166 (matched)
#>  - target estimand: ATT
#>  - covariates: age, educ, re74, re75, married, race
summary(m.out3, un = FALSE)
#> 
#> Call:
#> matchit(formula = treat ~ age + educ + re74 + re75, data = lalonde, 
#>     method = "nearest", distance = "mahalanobis", exact = ~married + 
#>         race, caliper = c(re74 = 0.2, re75 = 0.15))
#> 
#> Summary of Balance for Matched Data:
#>            Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean
#> age              24.7711       24.6024          0.0236     0.4804    0.0753
#> educ             10.3976       10.0482          0.1738     0.4610    0.0361
#> re74            705.7326      864.0131         -0.0324     0.9038    0.0270
#> re75            592.6384      640.9355         -0.0150     1.0194    0.0150
#> married           0.1687        0.1687          0.0000          .    0.0000
#> raceblack         0.7349        0.7349          0.0000          .    0.0000
#> racehispan        0.0723        0.0723          0.0000          .    0.0000
#> racewhite         0.1928        0.1928          0.0000          .    0.0000
#>            eCDF Max Std. Pair Dist.
#> age          0.2410          0.7544
#> educ         0.0964          0.7131
#> re74         0.2410          0.0537
#> re75         0.1084          0.0358
#> married      0.0000          0.0000
#> raceblack    0.0000          0.0000
#> racehispan   0.0000          0.0000
#> racewhite    0.0000          0.0000
#> 
#> Sample Sizes:
#>           Control Treated
#> All           429     185
#> Matched        83      83
#> Unmatched     346     102
#> Discarded       0       0
#> 

# 2:1 variable ratio NN matching on the PS
m.out4 <- matchit(treat ~ age + educ + race + nodegree +
                    married + re74 + re75,
                  data = lalonde,
                  method = "nearest",
                  ratio = 2,
                  min.controls = 1,
                  max.controls = 12)
m.out4
#> A `matchit` object
#>  - method: Variable ratio 2:1 nearest neighbor matching without replacement
#>  - distance: Propensity score
#>              - estimated with logistic regression
#>  - number of obs.: 614 (original), 555 (matched)
#>  - target estimand: ATT
#>  - covariates: age, educ, race, nodegree, married, re74, re75
summary(m.out4, un = FALSE)
#> 
#> Call:
#> matchit(formula = treat ~ age + educ + race + nodegree + married + 
#>     re74 + re75, data = lalonde, method = "nearest", ratio = 2, 
#>     min.controls = 1, max.controls = 12)
#> 
#> Summary of Balance for Matched Data:
#>            Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean
#> distance          0.5774        0.3608          0.9835     0.7425    0.1414
#> age              25.8162       25.4459          0.0518     0.4560    0.0830
#> educ             10.3459       10.4887         -0.0710     0.5654    0.0236
#> raceblack         0.8432        0.4703          1.0259          .    0.3730
#> racehispan        0.0595        0.2110         -0.6408          .    0.1515
#> racewhite         0.0973        0.3187         -0.7472          .    0.2214
#> nodegree          0.7081        0.6517          0.1240          .    0.0564
#> married           0.1892        0.2218         -0.0833          .    0.0326
#> re74           2095.5737     2614.8070         -0.1063     1.1615    0.0611
#> re75           1532.0553     1714.9169         -0.0568     1.3486    0.0498
#>            eCDF Max Std. Pair Dist.
#> distance     0.4216          0.5600
#> age          0.2457          1.4044
#> educ         0.0659          1.3187
#> raceblack    0.3730          0.5129
#> racehispan   0.1515          1.0629
#> racewhite    0.2214          0.8390
#> nodegree     0.0564          1.0760
#> married      0.0326          1.0420
#> re74         0.2912          0.9210
#> re75         0.2020          0.8587
#> 
#> Sample Sizes:
#>               Control Treated
#> All               429     185
#> Matched (ESS)     202     185
#> Matched           370     185
#> Unmatched          59       0
#> Discarded           0       0
#> 

# Some units received 1 match and some received 12
table(table(m.out4$subclass[m.out4$treat == 0]))
#> 
#>   1  10  12 
#> 168   1  16 
```
