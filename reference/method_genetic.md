# Genetic Matching

In
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md),
setting `method = "genetic"` performs genetic matching. Genetic matching
is a form of nearest neighbor matching where distances are computed as
the generalized Mahalanobis distance, which is a generalization of the
Mahalanobis distance with a scaling factor for each covariate that
represents the importance of that covariate to the distance. A genetic
algorithm is used to select the scaling factors. The scaling factors are
chosen as those which maximize a criterion related to covariate balance,
which can be chosen, but which by default is the smallest p-value in
covariate balance tests among the covariates. This method relies on and
is a wrapper for
[`Matching::GenMatch()`](https://rdrr.io/pkg/Matching/man/GenMatch.html)
and [`Matching::Match()`](https://rdrr.io/pkg/Matching/man/Match.html),
which use
[`rgenoud::genoud()`](https://rdrr.io/pkg/rgenoud/man/genoud.html) to
perform the optimization using the genetic algorithm.

This page details the allowable arguments with `method = "genetic"`. See
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
for an explanation of what each argument means in a general context and
how it can be specified.

Below is how
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
is used for genetic matching:


    matchit(formula,
            data = NULL,
            method = "genetic",
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
            ratio = 1,
            verbose = FALSE,
            ...) 

## Arguments

- formula:

  a two-sided [formula](https://rdrr.io/r/stats/formula.html) object
  containing the treatment and covariates to be used in creating the
  distance measure used in the matching. This formula will be supplied
  to the functions that estimate the distance measure and is used to
  determine the covariates whose balance is to be optimized.

- data:

  a data frame containing the variables named in `formula`. If not found
  in `data`, the variables will be sought in the environment.

- method:

  set here to `"genetic"`.

- distance:

  the distance measure to be used. See
  [`distance`](https://kosukeimai.github.io/MatchIt/reference/distance.md)
  for allowable options. When set to a method of estimating propensity
  scores or a numeric vector of distance values, the distance measure is
  included with the covariates in `formula` to be supplied to the
  generalized Mahalanobis distance matrix unless `mahvars` is specified.
  Otherwise, only the covariates in `formula` are supplied to the
  generalized Mahalanobis distance matrix to have their scaling factors
  chosen. `distance` *cannot* be supplied as a distance matrix.
  Supplying any method of computing a distance matrix (e.g.,
  `"mahalanobis"`) has the same effect of omitting propensity score but
  does not affect how the distance between units is computed otherwise.

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

  when a distance corresponds to a propensity score (e.g., for caliper
  matching or to discard units for common support), which covariates
  should be supplied to the generalized Mahalanobis distance matrix for
  matching. If unspecified, all variables in `formula` will be supplied
  to the distance matrix. Use `mahvars` to only supply a subset. Even if
  `mahvars` is specified, balance will be optimized on all covariates in
  `formula`. See Details.

- antiexact:

  for which variables anti-exact matching should take place. Anti-exact
  matching is processed using the `restrict` argument to
  [`Matching::GenMatch()`](https://rdrr.io/pkg/Matching/man/GenMatch.html)
  and
  [`Matching::Match()`](https://rdrr.io/pkg/Matching/man/Match.html).

- discard:

  a string containing a method for discarding units outside a region of
  common support. Only allowed when `distance` corresponds to a
  propensity score.

- reestimate:

  if `discard` is not `"none"`, whether to re-estimate the propensity
  score in the remaining sample prior to matching.

- s.weights:

  the variable containing sampling weights to be incorporated into
  propensity score models and balance statistics. These are also
  supplied to `GenMatch()` for use in computing the balance t-test
  p-values in the process of matching.

- replace:

  whether matching should be done with replacement.

- m.order:

  the order that the matching takes place. Allowable options include
  `"largest"`, where matching takes place in descending order of
  distance measures; `"smallest"`, where matching takes place in
  ascending order of distance measures; `"random"`, where matching takes
  place in a random order; and `"data"` where matching takes place based
  on the order of units in the data. When `m.order = "random"`, results
  may differ across different runs of the same code unless a seed is set
  and specified with [`set.seed()`](https://rdrr.io/r/base/Random.html).
  The default of `NULL` corresponds to `"largest"` when a propensity
  score is estimated or supplied as a vector and `"data"` otherwise.

- caliper:

  the width(s) of the caliper(s) used for caliper matching. See Details
  and Examples.

- std.caliper:

  `logical`; when calipers are specified, whether they are in standard
  deviation units (`TRUE`) or raw units (`FALSE`).

- ratio:

  how many control units should be matched to each treated unit for k:1
  matching. Should be a single integer value.

- verbose:

  `logical`; whether information about the matching process should be
  printed to the console. When `TRUE`, output from `GenMatch()` with
  `print.level = 2` will be displayed. Default is `FALSE` for no
  printing other than warnings.

- ...:

  additional arguments passed to
  [`Matching::GenMatch()`](https://rdrr.io/pkg/Matching/man/GenMatch.html).
  Potentially useful options include `pop.size`, `max.generations`, and
  `fit.func`. If `pop.size` is not specified, a warning from *Matching*
  will be thrown reminding you to change it. Note that the `ties` and
  `CommonSupport` arguments are set to `FALSE` and cannot be changed. If
  `distance.tolerance` is not specified, it is set to 0, whereas the
  default in *Matching* is 1e-5.

## Details

In genetic matching, covariates play three roles: 1) as the variables on
which balance is optimized, 2) as the variables in the generalized
Mahalanobis distance between units, and 3) in estimating the propensity
score. Variables supplied to `formula` are always used for role (1), as
the variables on which balance is optimized. When `distance` corresponds
to a propensity score, the covariates are also used to estimate the
propensity score (unless it is supplied). When `mahvars` is specified,
the named variables will form the covariates that go into the distance
matrix. Otherwise, the variables in `formula` along with the propensity
score will go into the distance matrix. This leads to three ways to use
`distance` and `mahvars` to perform the matching:

1.  When `distance` corresponds to a propensity score and `mahvars` *is
    not* specified, the covariates in `formula` along with the
    propensity score are used to form the generalized Mahalanobis
    distance matrix. This is the default and most typical use of
    `method = "genetic"` in
    [`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md).

2.  When `distance` corresponds to a propensity score and `mahvars` *is*
    specified, the covariates in `mahvars` are used to form the
    generalized Mahalanobis distance matrix. The covariates in `formula`
    are used to estimate the propensity score and have their balance
    optimized by the genetic algorithm. The propensity score is not
    included in the generalized Mahalanobis distance matrix.

3.  When `distance` is a method of computing a distance matrix
    (e.g.,`"mahalanobis"`), no propensity score is estimated, and the
    covariates in `formula` are used to form the generalized Mahalanobis
    distance matrix. Which specific method is supplied has no bearing on
    how the distance matrix is computed; it simply serves as a signal to
    omit estimation of a propensity score.

When a caliper is specified, any variables mentioned in `caliper`,
possibly including the propensity score, will be added to the matching
variables used to form the generalized Mahalanobis distance matrix. This
is because *Matching* doesn't allow for the separation of caliper
variables and matching variables in genetic matching.

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
Note that while `GenMatch()` and `Match()` support the ATE as an
estimand,
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
only supports the ATT and ATC for genetic matching.

### Reproducibility

Genetic matching involves a random component, so a seed must be set
using [`set.seed()`](https://rdrr.io/r/base/Random.html) to ensure
reproducibility. When `cluster` is used for parallel processing, the
seed must be compatible with parallel processing (e.g., by setting
`kind = "L'Ecuyer-CMRG"`).

## Outputs

All outputs described in
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
are returned with `method = "genetic"`. When `replace = TRUE`, the
`subclass` component is omitted. When `include.obj = TRUE` in the call
to
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md),
the output of the call to
[`Matching::GenMatch()`](https://rdrr.io/pkg/Matching/man/GenMatch.html)
will be included in the output.

## References

In a manuscript, be sure to cite the following papers if using
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
with `method = "genetic"`:

Diamond, A., & Sekhon, J. S. (2013). Genetic matching for estimating
causal effects: A general multivariate matching method for achieving
balance in observational studies. Review of Economics and Statistics,
95(3), 932–945.
[doi:10.1162/REST_a_00318](https://doi.org/10.1162/REST_a_00318)

Sekhon, J. S. (2011). Multivariate and Propensity Score Matching
Software with Automated Balance Optimization: The Matching package for
R. Journal of Statistical Software, 42(1), 1–52.
[doi:10.18637/jss.v042.i07](https://doi.org/10.18637/jss.v042.i07)

For example, a sentence might read:

*Genetic matching was performed using the MatchIt package (Ho, Imai,
King, & Stuart, 2011) in R, which calls functions from the Matching
package (Diamond & Sekhon, 2013; Sekhon, 2011).*

## See also

[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
for a detailed explanation of the inputs and outputs of a call to
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md).

[`Matching::GenMatch()`](https://rdrr.io/pkg/Matching/man/GenMatch.html)
and [`Matching::Match()`](https://rdrr.io/pkg/Matching/man/Match.html),
which do the work.

## Examples

``` r
data("lalonde")

# 1:1 genetic matching with PS as a covariate
m.out1 <- matchit(treat ~ age + educ + race + nodegree +
                    married + re74 + re75,
                  data = lalonde,
                  method = "genetic",
                  pop.size = 10) #use much larger pop.size
m.out1
#> A `matchit` object
#>  - method: 1:1 genetic matching without replacement
#>  - distance: Propensity score
#>              - estimated with logistic regression
#>  - number of obs.: 614 (original), 370 (matched)
#>  - target estimand: ATT
#>  - covariates: age, educ, race, nodegree, married, re74, re75
summary(m.out1)
#> 
#> Call:
#> matchit(formula = treat ~ age + educ + race + nodegree + married + 
#>     re74 + re75, data = lalonde, method = "genetic", pop.size = 10)
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
#> distance          0.5774        0.3510          1.0280     0.6880    0.1748
#> age              25.8162       25.9135         -0.0136     0.4954    0.0670
#> educ             10.3459       10.0865          0.1290     0.5967    0.0245
#> raceblack         0.8432        0.4703          1.0259          .    0.3730
#> racehispan        0.0595        0.2811         -0.9372          .    0.2216
#> racewhite         0.0973        0.2486         -0.5107          .    0.1514
#> nodegree          0.7081        0.6486          0.1308          .    0.0595
#> married           0.1892        0.3405         -0.3864          .    0.1514
#> re74           2095.5737     3393.9055         -0.2657     0.7718    0.0949
#> re75           1532.0553     2011.2096         -0.1488     1.0304    0.0713
#>            eCDF Max Std. Pair Dist.
#> distance     0.4216          1.0326
#> age          0.1946          1.1664
#> educ         0.0703          0.4463
#> raceblack    0.3730          1.0259
#> racehispan   0.2216          1.3943
#> racewhite    0.1514          0.5107
#> nodegree     0.0595          0.1308
#> married      0.1514          0.6625
#> re74         0.2973          0.3498
#> re75         0.2054          0.4251
#> 
#> Sample Sizes:
#>           Control Treated
#> All           429     185
#> Matched       185     185
#> Unmatched     244       0
#> Discarded       0       0
#> 

# 2:1 genetic matching with replacement without PS
m.out2 <- matchit(treat ~ age + educ + race + nodegree +
                    married + re74 + re75,
                  data = lalonde,
                  method = "genetic",
                  replace = TRUE,
                  ratio = 2,
                  distance = "mahalanobis",
                  pop.size = 10) #use much larger pop.size
m.out2
#> A `matchit` object
#>  - method: 2:1 genetic matching with replacement
#>  - distance: Mahalanobis - number of obs.: 614 (original), 302 (matched)
#>  - target estimand: ATT
#>  - covariates: age, educ, race, nodegree, married, re74, re75
summary(m.out2, un = FALSE)
#> 
#> Call:
#> matchit(formula = treat ~ age + educ + race + nodegree + married + 
#>     re74 + re75, data = lalonde, method = "genetic", distance = "mahalanobis", 
#>     replace = TRUE, ratio = 2, pop.size = 10)
#> 
#> Summary of Balance for Matched Data:
#>            Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean
#> age              25.8162       25.0595          0.1058     0.7234    0.0412
#> educ             10.3459       10.3189          0.0134     1.0229    0.0074
#> raceblack         0.8432        0.8378          0.0149          .    0.0054
#> racehispan        0.0595        0.0595          0.0000          .    0.0000
#> racewhite         0.0973        0.1027         -0.0182          .    0.0054
#> nodegree          0.7081        0.7081          0.0000          .    0.0000
#> married           0.1892        0.1622          0.0690          .    0.0270
#> re74           2095.5737     2071.5698          0.0049     1.3109    0.0287
#> re75           1532.0553     1104.6774          0.1328     1.7164    0.0382
#>            eCDF Max Std. Pair Dist.
#> age          0.1703          0.4397
#> educ         0.0324          0.2420
#> raceblack    0.0054          0.0149
#> racehispan   0.0000          0.0000
#> racewhite    0.0054          0.0182
#> nodegree     0.0000          0.0000
#> married      0.0270          0.1104
#> re74         0.1757          0.2469
#> re75         0.0865          0.2999
#> 
#> Sample Sizes:
#>               Control Treated
#> All            429.       185
#> Matched (ESS)   40.74     185
#> Matched        117.       185
#> Unmatched      312.         0
#> Discarded        0.         0
#> 

# 1:1 genetic matching on just age, educ, re74, and re75
# within calipers on PS and educ; other variables are
# used to estimate PS
m.out3 <- matchit(treat ~ age + educ + race + nodegree +
                    married + re74 + re75,
                  data = lalonde,
                  method = "genetic",
                  mahvars = ~ age + educ + re74 + re75,
                  caliper = c(.05, educ = 2),
                  std.caliper = c(TRUE, FALSE),
                  pop.size = 10) #use much larger pop.size
m.out3
#> A `matchit` object
#>  - method: 1:1 genetic matching without replacement
#>  - distance: Mahalanobis [matching]
#>              Propensity score [caliper]
#> 
#>              - estimated with logistic regression
#>  - caliper: <distance> (0.015), educ (2)
#>  - number of obs.: 614 (original), 206 (matched)
#>  - target estimand: ATT
#>  - covariates: age, educ, race, nodegree, married, re74, re75
summary(m.out3, un = FALSE)
#> 
#> Call:
#> matchit(formula = treat ~ age + educ + race + nodegree + married + 
#>     re74 + re75, data = lalonde, method = "genetic", mahvars = ~age + 
#>     educ + re74 + re75, caliper = c(0.05, educ = 2), std.caliper = c(TRUE, 
#>     FALSE), pop.size = 10)
#> 
#> Summary of Balance for Matched Data:
#>            Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean
#> distance          0.4937        0.4870          0.0305     1.0220    0.0105
#> age              25.5340       24.9612          0.0801     0.4989    0.0784
#> educ             10.1553       10.2913         -0.0676     0.9739    0.0133
#> raceblack         0.7184        0.6990          0.0534          .    0.0194
#> racehispan        0.1068        0.0971          0.0411          .    0.0097
#> racewhite         0.1748        0.2039         -0.0983          .    0.0291
#> nodegree          0.7087        0.6699          0.0854          .    0.0388
#> married           0.2233        0.2039          0.0496          .    0.0194
#> re74           2681.1038     2015.1915          0.1363     2.2258    0.0576
#> re75           1891.5406     1535.6796          0.1105     1.9425    0.0316
#>            eCDF Max Std. Pair Dist.
#> distance     0.0777          0.0390
#> age          0.2718          0.9186
#> educ         0.0680          0.3573
#> raceblack    0.0194          0.1068
#> racehispan   0.0097          0.3695
#> racewhite    0.0291          0.1638
#> nodegree     0.0388          0.3844
#> married      0.0194          0.5454
#> re74         0.2524          0.5029
#> re75         0.1068          0.5474
#> 
#> Sample Sizes:
#>           Control Treated
#> All           429     185
#> Matched       103     103
#> Unmatched     326      82
#> Discarded       0       0
#> 
```
