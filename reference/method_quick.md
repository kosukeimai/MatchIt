# Fast Generalized Full Matching

In
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md),
setting `method = "quick"` performs generalized full matching, which is
a form of subclassification wherein all units, both treatment and
control (i.e., the "full" sample), are assigned to a subclass and
receive at least one match. It uses an algorithm that is extremely fast
compared to optimal full matching, which is why it is labeled as
"quick", at the expense of true optimality. The method is described in
Sävje, Higgins, & Sekhon (2021). The method relies on and is a wrapper
for
[`quickmatch::quickmatch()`](https://rdrr.io/pkg/quickmatch/man/quickmatch.html).

Advantages of generalized full matching include that the matching order
is not required to be specified, units do not need to be discarded, and
it is less likely that extreme within-subclass distances will be large,
unlike with standard subclassification. The primary output of
generalized full matching is a set of matching weights that can be
applied to the matched sample; in this way, generalized full matching
can be seen as a robust alternative to propensity score weighting,
robust in the sense that the propensity score model does not need to be
correct to estimate the treatment effect without bias.

This page details the allowable arguments with `method = "quick"`. See
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
for an explanation of what each argument means in a general context and
how it can be specified.

Below is how
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
is used for generalized full matching:

    matchit(formula,
            data = NULL,
            method = "quick",
            distance = "glm",
            link = "logit",
            distance.options = list(),
            estimand = "ATT",
            exact = NULL,
            mahvars = NULL,
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

  set here to `"quick"`.

- distance:

  the distance measure to be used. See
  [`distance`](https://kosukeimai.github.io/MatchIt/reference/distance.md)
  for allowable options. Cannot be supplied as a matrix.

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
  when `distance` corresponds to a propensity score (e.g., to discard
  units for common support). If specified, the distance measure will not
  be used in matching.

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

  the width of the caliper used for caliper matching. A caliper can only
  be placed on the propensity score and cannot be negative.

- std.caliper:

  `logical`; when a caliper is specified, whether it is in standard
  deviation units (`TRUE`) or raw units (`FALSE`).

- verbose:

  `logical`; whether information about the matching process should be
  printed to the console.

- ...:

  additional arguments passed to
  [`quickmatch::quickmatch()`](https://rdrr.io/pkg/quickmatch/man/quickmatch.html).
  Allowed arguments include `treatment_constraints`, `size_constraint`,
  `target`, and other arguments passed to
  [`scclust::sc_clustering()`](https://rdrr.io/pkg/scclust/man/sc_clustering.html)
  (see
  [`quickmatch::quickmatch()`](https://rdrr.io/pkg/quickmatch/man/quickmatch.html)
  for details). In particular, changing `seed_method` from its default
  can improve performance. No arguments will be passed to
  [`distances::distances()`](https://rdrr.io/pkg/distances/man/distances.html).

  The arguments `replace`, `ratio`, `min.controls`, `max.controls`,
  `m.order`, and `antiexact` are ignored with a warning.

## Details

Generalized full matching is similar to optimal full matching, but has
some additional flexibility that can be controlled by some of the extra
arguments available. By default, `method = "quick"` performs a standard
full match in which all units are matched (unless restricted by the
caliper) and assigned to a subclass. Each subclass could contain
multiple units from each treatment group. The subclasses are chosen to
minimize the largest within-subclass distance between units (including
between units of the same treatment group). Notably, generalized full
matching requires less memory and can run much faster than optimal full
matching and optimal pair matching and, in some cases, even than nearest
neighbor matching, and it can be used with huge datasets (e.g., in the
millions) while running in under a minute.

## Outputs

All outputs described in
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
are returned with `method = "quick"` except for `match.matrix`. This is
because matching strata are not indexed by treated units as they are in
some other forms of matching. When `include.obj = TRUE` in the call to
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md),
the output of the call to
[`quickmatch::quickmatch()`](https://rdrr.io/pkg/quickmatch/man/quickmatch.html)
will be included in the output. When `exact` is specified, this will be
a list of such objects, one for each stratum of the `exact` variables.

## References

In a manuscript, be sure to cite the *quickmatch* package if using
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
with `method = "quick"`. A citation can be generated using
`citation("quickmatch")`.

For example, a sentence might read:

*Generalized full matching was performed using the MatchIt package (Ho,
Imai, King, & Stuart, 2011) in R, which calls functions from the
quickmatch package (Sävje, Sekhon, & Higgins, 2024).*

You should also cite the following paper, which develops and describes
the method:

Sävje, F., Higgins, M. J., & Sekhon, J. S. (2021). Generalized Full
Matching. *Political Analysis*, 29(4), 423–447.
[doi:10.1017/pan.2020.32](https://doi.org/10.1017/pan.2020.32)

## See also

[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
for a detailed explanation of the inputs and outputs of a call to
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md).

[`quickmatch::quickmatch()`](https://rdrr.io/pkg/quickmatch/man/quickmatch.html),
which is the workhorse.

[`method_full`](https://kosukeimai.github.io/MatchIt/reference/method_full.md)
for optimal full matching, which is nearly the same but offers more
customizability and more optimal solutions at the cost of speed.

## Examples

``` r
data("lalonde")

# Generalized full PS matching
m.out1 <- matchit(treat ~ age + educ + race + nodegree +
                    married + re74 + re75,
                  data = lalonde,
                  method = "quick")
m.out1
#> A `matchit` object
#>  - method: Generalized full matching
#>  - distance: Propensity score
#>              - estimated with logistic regression
#>  - number of obs.: 614 (original), 614 (matched)
#>  - target estimand: ATT
#>  - covariates: age, educ, race, nodegree, married, re74, re75
summary(m.out1)
#> 
#> Call:
#> matchit(formula = treat ~ age + educ + race + nodegree + married + 
#>     re74 + re75, data = lalonde, method = "quick")
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
#> distance          0.5774        0.5764          0.0047     0.9895    0.0044
#> age              25.8162       25.4476          0.0515     0.5099    0.0728
#> educ             10.3459       10.6781         -0.1652     0.5742    0.0251
#> raceblack         0.8432        0.8376          0.0155          .    0.0056
#> racehispan        0.0595        0.0607         -0.0052          .    0.0012
#> racewhite         0.0973        0.1017         -0.0148          .    0.0044
#> nodegree          0.7081        0.6602          0.1054          .    0.0479
#> married           0.1892        0.1355          0.1371          .    0.0537
#> re74           2095.5737     2873.4651         -0.1592     0.8886    0.0709
#> re75           1532.0553     1634.2731         -0.0318     1.8751    0.0736
#>            eCDF Max Std. Pair Dist.
#> distance     0.0541          0.0202
#> age          0.2659          1.2519
#> educ         0.0855          1.2506
#> raceblack    0.0056          0.0265
#> racehispan   0.0012          0.4970
#> racewhite    0.0044          0.3641
#> nodegree     0.0479          0.9324
#> married      0.0537          0.4870
#> re74         0.3004          0.8699
#> re75         0.2760          0.8428
#> 
#> Sample Sizes:
#>               Control Treated
#> All            429.       185
#> Matched (ESS)   50.74     185
#> Matched        429.       185
#> Unmatched        0.         0
#> Discarded        0.         0
#> 
```
