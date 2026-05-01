# Add sampling weights to a `matchit` object

Adds sampling weights to a `matchit` object so that they are
incorporated into balance assessment and creation of the weights. This
would typically only be used when an argument to `s.weights` was not
supplied to
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
(i.e., because they were not to be included in the estimation of the
propensity score) but sampling weights are required for generalizing an
effect to the correct population. Without adding sampling weights to the
`matchit` object, balance assessment tools (i.e.,
[`summary.matchit()`](https://kosukeimai.github.io/MatchIt/reference/summary.matchit.md)
and
[`plot.matchit()`](https://kosukeimai.github.io/MatchIt/reference/plot.matchit.md))
will not calculate balance statistics correctly, and the weights
produced by
[`match_data()`](https://kosukeimai.github.io/MatchIt/reference/match_data.md)
and
[`get_matches()`](https://kosukeimai.github.io/MatchIt/reference/match_data.md)
will not incorporate the sampling weights.

## Usage

``` r
add_s.weights(m, s.weights = NULL, data = NULL)
```

## Arguments

- m:

  a `matchit` object; the output of a call to
  [`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md),
  typically with the `s.weights` argument unspecified.

- s.weights:

  an numeric vector of sampling weights to be added to the `matchit`
  object. Can also be specified as a string containing the name of
  variable in `data` to be used or a one-sided formula with the variable
  on the right-hand side (e.g., `~ SW`).

- data:

  a data frame containing the sampling weights if given as a string or
  formula. If unspecified, `add_s.weights()` will attempt to find the
  dataset using the environment of the `matchit` object.

## Value

a `matchit` object with an `s.weights` component containing the supplied
sampling weights. If `s.weights = NULL`, the original `matchit` object
is returned.

## See also

[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md);
[`match_data()`](https://kosukeimai.github.io/MatchIt/reference/match_data.md)

## Author

Noah Greifer

## Examples

``` r

data("lalonde")

# Generate random sampling weights, just
# for this example
sw <- rchisq(nrow(lalonde), 2)

# NN PS match using logistic regression PS that doesn't
# include sampling weights
m.out <- matchit(treat ~ age + educ + race + nodegree +
                   married  + re74 + re75,
                 data = lalonde)

m.out
#> A `matchit` object
#>  - method: 1:1 nearest neighbor matching without replacement
#>  - distance: Propensity score
#>              - estimated with logistic regression
#>  - number of obs.: 614 (original), 370 (matched)
#>  - target estimand: ATT
#>  - covariates: age, educ, race, nodegree, married, re74, re75

# Add s.weights to the matchit object
m.out <- add_s.weights(m.out, sw)

m.out #note additional output
#> A `matchit` object
#>  - method: 1:1 nearest neighbor matching without replacement
#>  - distance: Propensity score
#>              - estimated with logistic regression
#>              - sampling weights not included in estimation
#>  - number of obs.: 614 (original), 370 (matched)
#>  - sampling weights: present
#>  - target estimand: ATT
#>  - covariates: age, educ, race, nodegree, married, re74, re75

# Check balance; note that sample sizes incorporate
# s.weights
summary(m.out, improvement = FALSE)
#> 
#> Call:
#> matchit(formula = treat ~ age + educ + race + nodegree + married + 
#>     re74 + re75, data = lalonde)
#> 
#> Summary of Balance for All Data:
#>            Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean
#> distance          0.5796        0.1821          1.7726     0.9548    0.3817
#> age              25.7457       27.6471         -0.2776     0.3916    0.0917
#> educ             10.2366       10.1874          0.0242     0.5375    0.0324
#> raceblack         0.8523        0.1965          1.8482          .    0.6558
#> racehispan        0.0520        0.1497         -0.4397          .    0.0977
#> racewhite         0.0957        0.6538         -1.8974          .    0.5581
#> nodegree          0.7233        0.5991          0.2777          .    0.1242
#> married           0.2114        0.5043         -0.7175          .    0.2929
#> re74           1901.2003     5724.2451         -0.7745     0.4560    0.2349
#> re75           1423.0173     2586.7575         -0.3786     0.6725    0.1432
#>            eCDF Max
#> distance     0.6558
#> age          0.1616
#> educ         0.1242
#> raceblack    0.6558
#> racehispan   0.0977
#> racewhite    0.5581
#> nodegree     0.1242
#> married      0.2929
#> re74         0.4843
#> re75         0.3168
#> 
#> Summary of Balance for Matched Data:
#>            Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean
#> distance          0.5796        0.3554          0.9999     0.7895    0.1383
#> age              25.7457       24.6876          0.1545     0.4187    0.0942
#> educ             10.2366       10.5470         -0.1526     0.6304    0.0230
#> raceblack         0.8523        0.4415          1.1577          .    0.4108
#> racehispan        0.0520        0.2361         -0.8285          .    0.1840
#> racewhite         0.0957        0.3225         -0.7710          .    0.2268
#> nodegree          0.7233        0.6463          0.1722          .    0.0770
#> married           0.2114        0.1909          0.0501          .    0.0205
#> re74           1901.2003     2115.7108         -0.0435     1.6183    0.0576
#> re75           1423.0173     1831.8247         -0.1330     1.0301    0.0653
#>            eCDF Max Std. Pair Dist.
#> distance     0.4342          0.9567
#> age          0.3460          1.4563
#> educ         0.0770          1.2332
#> raceblack    0.4108          1.0511
#> racehispan   0.1840          1.1438
#> racewhite    0.2268          0.8453
#> nodegree     0.0770          1.0271
#> married      0.0205          0.7944
#> re74         0.3418          0.7885
#> re75         0.2391          0.7731
#> 
#> Sample Sizes:
#>               Control Treated
#> All (ESS)      214.55   93.95
#> All            429.    185.  
#> Matched (ESS)   97.95   93.95
#> Matched        185.    185.  
#> Unmatched      244.      0.  
#> Discarded        0.      0.  
#> 
```
