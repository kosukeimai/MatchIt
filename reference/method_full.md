# Optimal Full Matching

In
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md),
setting `method = "full"` performs optimal full matching, which is a
form of subclassification wherein all units, both treatment and control
(i.e., the "full" sample), are assigned to a subclass and receive at
least one match. The matching is optimal in the sense that that sum of
the absolute distances between the treated and control units in each
subclass is as small as possible. The method relies on and is a wrapper
for
[`optmatch::fullmatch()`](https://rdrr.io/pkg/optmatch/man/fullmatch.html).

Advantages of optimal full matching include that the matching order is
not required to be specified, units do not need to be discarded, and it
is less likely that extreme within-subclass distances will be large,
unlike with standard subclassification. The primary output of full
matching is a set of matching weights that can be applied to the matched
sample; in this way, full matching can be seen as a robust alternative
to propensity score weighting, robust in the sense that the propensity
score model does not need to be correct to estimate the treatment effect
without bias. Note: with large samples, the optimization may fail or run
very slowly; one can try using
[`method = "quick"`](https://kosukeimai.github.io/MatchIt/reference/method_quick.md)
instead, which also performs full matching but can be much faster.

This page details the allowable arguments with `method = "full"`. See
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
for an explanation of what each argument means in a general context and
how it can be specified.

Below is how
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
is used for optimal full matching:


    matchit(formula,
            data = NULL,
            method = "full",
            distance = "glm",
            link = "logit",
            distance.options = list(),
            estimand = "ATT",
            exact = NULL,
            mahvars = NULL,
            anitexact = NULL,
            discard = "none",
            reestimate = FALSE,
            s.weights = NULL,
            caliper = NULL,
            std.caliper = TRUE,
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

  set here to `"full"`.

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
  `"ATT"`, `"ATC"`, and `"ATE"`. The estimand controls how the weights
  are computed; see the Computing Weights section at
  [`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
  for details.

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
  common support. Only allowed when `distance` corresponds to a
  propensity score.

- reestimate:

  if `discard` is not `"none"`, whether to re-estimate the propensity
  score in the remaining sample prior to matching.

- s.weights:

  the variable containing sampling weights to be incorporated into
  propensity score models and balance statistics.

- caliper:

  the width(s) of the caliper(s) used for caliper matching. Calipers are
  processed by
  [`optmatch::caliper()`](https://rdrr.io/pkg/optmatch/man/caliper-methods.html).
  Positive and negative calipers are allowed. See Notes and Examples.

- std.caliper:

  `logical`; when calipers are specified, whether they are in standard
  deviation units (`TRUE`) or raw units (`FALSE`).

- verbose:

  `logical`; whether information about the matching process should be
  printed to the console.

- ...:

  additional arguments passed to
  [`optmatch::fullmatch()`](https://rdrr.io/pkg/optmatch/man/fullmatch.html).
  Allowed arguments include `min.controls`, `max.controls`,
  `omit.fraction`, `mean.controls`, `tol`, and `solver`. See the
  [`optmatch::fullmatch()`](https://rdrr.io/pkg/optmatch/man/fullmatch.html)
  documentation for details. In general, `tol` should be set to a low
  number (e.g., `1e-7`) to get a more precise solution.

  The arguments `replace`, `m.order`, and `ratio` are ignored with a
  warning.

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
    propensity score or a vector of distance measures, i.e., it should
    not be `"mahalanobis"`. Use `mahvars` to specify the variables used
    to create the Mahalanobis distance. For example, to perform
    Mahalanobis within a propensity score caliper, the following could
    be run:


        matchit(treat ~ X1 + X2 + X3, method = "nearest",
                distance =  "glm", caliper = .25,
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

## Note

Calipers can only be used when `min.controls` is left at its default.

The option `"optmatch_max_problem_size"` is automatically set to `Inf`
during the matching process, different from its default in *optmatch*.
This enables matching problems of any size to be run, but may also let
huge, infeasible problems get through and potentially take a long time
or crash R. See
[`optmatch::setMaxProblemSize()`](https://rdrr.io/pkg/optmatch/man/setMaxProblemSize.html)
for more details.

## Outputs

All outputs described in
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
are returned with `method = "full"` except for `match.matrix`. This is
because matching strata are not indexed by treated units as they are in
some other forms of matching. When `include.obj = TRUE` in the call to
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md),
the output of the call to
[`optmatch::fullmatch()`](https://rdrr.io/pkg/optmatch/man/fullmatch.html)
will be included in the output. When `exact` is specified, this will be
a list of such objects, one for each stratum of the `exact` variables.

## References

In a manuscript, be sure to cite the following paper if using
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
with `method = "full"`:

Hansen, B. B., & Klopfer, S. O. (2006). Optimal Full Matching and
Related Designs via Network Flows. *Journal of Computational and
Graphical Statistics*, 15(3), 609–627.
[doi:10.1198/106186006X137047](https://doi.org/10.1198/106186006X137047)

For example, a sentence might read:

*Optimal full matching was performed using the MatchIt package (Ho,
Imai, King, & Stuart, 2011) in R, which calls functions from the
optmatch package (Hansen & Klopfer, 2006).*

Theory is also developed in the following article:

Hansen, B. B. (2004). Full Matching in an Observational Study of
Coaching for the SAT. Journal of the American Statistical Association,
99(467), 609–618.
[doi:10.1198/016214504000000647](https://doi.org/10.1198/016214504000000647)

## See also

[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
for a detailed explanation of the inputs and outputs of a call to
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md).

[`optmatch::fullmatch()`](https://rdrr.io/pkg/optmatch/man/fullmatch.html),
which is the workhorse.

[`method_optimal`](https://kosukeimai.github.io/MatchIt/reference/method_optimal.md)
for optimal pair matching, which is a special case of optimal full
matching, and which relies on similar machinery. Results from
`method = "optimal"` can be replicated with `method = "full"` by setting
`min.controls`, `max.controls`, and `mean.controls` to the desired
`ratio`.

[`method_quick`](https://kosukeimai.github.io/MatchIt/reference/method_quick.md)
for fast generalized quick matching, which is very similar to optimal
full matching but can be dramatically faster at the expense of
optimality and is less customizable.

## Examples

``` r
data("lalonde")

# Optimal full PS matching
m.out1 <- matchit(treat ~ age + educ + race + nodegree +
                    married + re74 + re75,
                  data = lalonde,
                  method = "full")
m.out1
#> A `matchit` object
#>  - method: Optimal full matching
#>  - distance: Propensity score
#>              - estimated with logistic regression
#>  - number of obs.: 614 (original), 614 (matched)
#>  - target estimand: ATT
#>  - covariates: age, educ, race, nodegree, married, re74, re75
summary(m.out1)
#> 
#> Call:
#> matchit(formula = treat ~ age + educ + race + nodegree + married + 
#>     re74 + re75, data = lalonde, method = "full")
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
#> distance          0.5774        0.5762          0.0054     0.9930    0.0041
#> age              25.8162       24.8095          0.1407     0.4976    0.0795
#> educ             10.3459       10.3452          0.0004     0.5830    0.0206
#> raceblack         0.8432        0.8347          0.0236          .    0.0086
#> racehispan        0.0595        0.0657         -0.0266          .    0.0063
#> racewhite         0.0973        0.0996         -0.0078          .    0.0023
#> nodegree          0.7081        0.7056          0.0056          .    0.0025
#> married           0.1892        0.1368          0.1338          .    0.0524
#> re74           2095.5737     2363.4473         -0.0548     1.1080    0.0424
#> re75           1532.0553     1632.4020         -0.0312     1.8588    0.0704
#>            eCDF Max Std. Pair Dist.
#> distance     0.0486          0.0192
#> age          0.3131          1.3111
#> educ         0.0548          1.2390
#> raceblack    0.0086          0.0324
#> racehispan   0.0063          0.5400
#> racewhite    0.0023          0.3911
#> nodegree     0.0025          0.9593
#> married      0.0524          0.4715
#> re74         0.2492          0.8654
#> re75         0.2366          0.8099
#> 
#> Sample Sizes:
#>               Control Treated
#> All            429.       185
#> Matched (ESS)   52.11     185
#> Matched        429.       185
#> Unmatched        0.         0
#> Discarded        0.         0
#> 

# Optimal full Mahalanobis distance matching within a PS caliper
m.out2 <- matchit(treat ~ age + educ + race + nodegree +
                    married + re74 + re75,
                  data = lalonde,
                  method = "full",
                  caliper = .01,
                  mahvars = ~ age + educ + re74 + re75)
m.out2
#> A `matchit` object
#>  - method: Optimal full matching
#>  - distance: Mahalanobis [matching]
#>              Propensity score [caliper]
#> 
#>              - estimated with logistic regression
#>  - caliper: <distance> (0.003)
#>  - number of obs.: 614 (original), 349 (matched)
#>  - target estimand: ATT
#>  - covariates: age, educ, race, nodegree, married, re74, re75
summary(m.out2, un = FALSE)
#> 
#> Call:
#> matchit(formula = treat ~ age + educ + race + nodegree + married + 
#>     re74 + re75, data = lalonde, method = "full", mahvars = ~age + 
#>     educ + re74 + re75, caliper = 0.01)
#> 
#> Summary of Balance for Matched Data:
#>            Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean
#> distance          0.5564        0.5566         -0.0012     0.9885    0.0028
#> age              25.1940       24.4645          0.1020     0.5323    0.0681
#> educ             10.3060       10.5529         -0.1228     0.6813    0.0193
#> raceblack         0.7985        0.8022         -0.0103          .    0.0037
#> racehispan        0.0672        0.0596          0.0321          .    0.0076
#> racewhite         0.1343        0.1382         -0.0131          .    0.0039
#> nodegree          0.6716        0.6742         -0.0057          .    0.0026
#> married           0.1642        0.1189          0.1157          .    0.0453
#> re74           1504.8003     2669.4440         -0.2383     0.4985    0.0889
#> re75           1242.9898     1479.5812         -0.0735     1.8382    0.0762
#>            eCDF Max Std. Pair Dist.
#> distance     0.0373          0.0063
#> age          0.2709          1.1149
#> educ         0.0812          1.1082
#> raceblack    0.0037          0.0105
#> racehispan   0.0076          0.5788
#> racewhite    0.0039          0.4490
#> nodegree     0.0026          0.8782
#> married      0.0453          0.4951
#> re74         0.3362          0.6849
#> re75         0.2368          0.6460
#> 
#> Sample Sizes:
#>               Control Treated
#> All            429.       185
#> Matched (ESS)   61.97     134
#> Matched        215.       134
#> Unmatched      214.        51
#> Discarded        0.         0
#> 

# Optimal full Mahalanobis distance matching within calipers
# of 500 on re74 and re75
m.out3 <- matchit(treat ~ age + educ + re74 + re75,
                  data = lalonde,
                  distance = "mahalanobis",
                  method = "full",
                  caliper = c(re74 = 500,
                              re75 = 500),
                  std.caliper = FALSE)
m.out3
#> A `matchit` object
#>  - method: Optimal full matching
#>  - distance: Mahalanobis - caliper: re74 (500), re75 (500)
#>  - number of obs.: 614 (original), 391 (matched)
#>  - target estimand: ATT
#>  - covariates: age, educ, re74, re75
summary(m.out3,
        addlvariables = ~race + nodegree + married,
        data = lalonde,
        un = FALSE)
#> 
#> Call:
#> matchit(formula = treat ~ age + educ + re74 + re75, data = lalonde, 
#>     method = "full", distance = "mahalanobis", caliper = c(re74 = 500, 
#>         re75 = 500), std.caliper = FALSE)
#> 
#> Summary of Balance for Matched Data:
#>            Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean
#> age              25.8438       26.6668         -0.1150     0.6314    0.0391
#> educ             10.2562       10.2998         -0.0217     0.5406    0.0266
#> re74            902.3844      927.5659         -0.0052     0.9952    0.0089
#> re75            661.1024      674.8033         -0.0043     1.0111    0.0090
#> raceblack         0.8500        0.3065          1.4948          .    0.5435
#> racehispan        0.0437        0.1258         -0.3469          .    0.0820
#> racewhite         0.1062        0.5677         -1.5570          .    0.4614
#> nodegree          0.7312        0.6768          0.1197          .    0.0544
#> married           0.1625        0.4717         -0.7895          .    0.3092
#>            eCDF Max Std. Pair Dist.
#> age          0.0868          0.5037
#> educ         0.0623          0.5850
#> re74         0.1957          0.0249
#> re75         0.0871          0.0421
#> raceblack    0.5435          1.7611
#> racehispan   0.0820          0.7605
#> racewhite    0.4614          1.8935
#> nodegree     0.0544          0.3086
#> married      0.3092          1.1113
#> 
#> Sample Sizes:
#>               Control Treated
#> All            429.       185
#> Matched (ESS)   54.18     160
#> Matched        231.       160
#> Unmatched      198.        25
#> Discarded        0.         0
#> 
```
