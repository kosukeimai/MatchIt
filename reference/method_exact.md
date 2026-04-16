# Exact Matching

In
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md),
setting `method = "exact"` performs exact matching. With exact matching,
a complete cross of the covariates is used to form subclasses defined by
each combination of the covariate levels. Any subclass that doesn't
contain both treated and control units is discarded, leaving only
subclasses containing treatment and control units that are exactly equal
on the included covariates. The benefits of exact matching are that
confounding due to the covariates included is completely eliminated,
regardless of the functional form of the treatment or outcome models.
The problem is that typically many units will be discarded, sometimes
dramatically reducing precision and changing the target population of
inference. To use exact matching in combination with another matching
method (i.e., to exact match on some covariates and some other form of
matching on others), use the `exact` argument with that method.

This page details the allowable arguments with `method = "exact"`. See
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
for an explanation of what each argument means in a general context and
how it can be specified.

Below is how
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
is used for exact matching:

    matchit(formula,
            data = NULL,
            method = "exact",
            estimand = "ATT",
            s.weights = NULL,
            verbose = FALSE,
            ...)

## Arguments

- formula:

  a two-sided [formula](https://rdrr.io/r/stats/formula.html) object
  containing the treatment and covariates to be used in creating the
  subclasses defined by a full cross of the covariate levels.

- data:

  a data frame containing the variables named in `formula`. If not found
  in `data`, the variables will be sought in the environment.

- method:

  set here to `"exact"`.

- estimand:

  a string containing the desired estimand. Allowable options include
  `"ATT"`, `"ATC"`, and `"ATE"`. The estimand controls how the weights
  are computed; see the Computing Weights section at
  [`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
  for details.

- s.weights:

  the variable containing sampling weights to be incorporated into
  balance statistics. These weights do not affect the matching process.

- verbose:

  `logical`; whether information about the matching process should be
  printed to the console.

- ...:

  ignored.

  The arguments `distance` (and related arguments), `exact`, `mahvars`,
  `discard` (and related arguments), `replace`, `m.order`, `caliper`
  (and related arguments), and `ratio` are ignored with a warning.

## Outputs

All outputs described in
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
are returned with `method = "exact"` except for `match.matrix`. This is
because matching strata are not indexed by treated units as they are in
some other forms of matching. `include.obj` is ignored.

## References

In a manuscript, you don't need to cite another package when using
`method = "exact"` because the matching is performed completely within
*MatchIt*. For example, a sentence might read:

*Exact matching was performed using the MatchIt package (Ho, Imai, King,
& Stuart, 2011) in R.*

## See also

[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
for a detailed explanation of the inputs and outputs of a call to
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md).
The `exact` argument can be used with other methods to perform exact
matching in combination with other matching methods.

[method_cem](https://kosukeimai.github.io/MatchIt/reference/method_cem.md)
for coarsened exact matching, which performs exact matching on coarsened
versions of the covariates.

## Examples

``` r
data("lalonde")

# Exact matching on age, race, married, and educ
m.out1 <- matchit(treat ~ age + race +
                    married + educ,
                  data = lalonde,
                  method = "exact")
m.out1
#> A `matchit` object
#>  - method: Exact matching
#>  - number of obs.: 614 (original), 113 (matched)
#>  - target estimand: ATT
#>  - covariates: age, race, married, educ
summary(m.out1)
#> 
#> Call:
#> matchit(formula = treat ~ age + race + married + educ, data = lalonde, 
#>     method = "exact")
#> 
#> Summary of Balance for All Data:
#>            Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean
#> age              25.8162       28.0303         -0.3094     0.4400    0.0813
#> raceblack         0.8432        0.2028          1.7615          .    0.6404
#> racehispan        0.0595        0.1422         -0.3498          .    0.0827
#> racewhite         0.0973        0.6550         -1.8819          .    0.5577
#> married           0.1892        0.5128         -0.8263          .    0.3236
#> educ             10.3459       10.2354          0.0550     0.4959    0.0347
#>            eCDF Max
#> age          0.1577
#> raceblack    0.6404
#> racehispan   0.0827
#> racewhite    0.5577
#> married      0.3236
#> educ         0.1114
#> 
#> Summary of Balance for Matched Data:
#>            Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean
#> age              19.9815       19.9815               0     0.9936         0
#> raceblack         0.7778        0.7778               0          .         0
#> racehispan        0.0370        0.0370               0          .         0
#> racewhite         0.1852        0.1852               0          .         0
#> married           0.0370        0.0370               0          .         0
#> educ             10.3333       10.3333               0     0.9936         0
#>            eCDF Max Std. Pair Dist.
#> age               0               0
#> raceblack         0               0
#> racehispan        0               0
#> racewhite         0               0
#> married           0               0
#> educ              0               0
#> 
#> Sample Sizes:
#>               Control Treated
#> All            429.       185
#> Matched (ESS)   40.36      54
#> Matched         59.        54
#> Unmatched      370.       131
#> Discarded        0.         0
#> 
```
