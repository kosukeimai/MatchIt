# Matching with Sampling Weights

## Introduction

Sampling weights (also known as survey weights) frequently appear when
using large, representative datasets. They are required to ensure any
estimated quantities generalize to a target population defined by the
weights. Evidence suggests that sampling weights need to be incorporated
into a propensity score matching analysis to obtain valid and unbiased
estimates of the treatment effect in the sampling weighted population
([DuGoff, Schuler, and Stuart 2014](#ref-dugoff2014); [Austin, Jembere,
and Chiu 2016](#ref-austin2016); [Lenis et al. 2019](#ref-lenis2019)).
In this guide, we demonstrate how to use sampling weights with `MatchIt`
for propensity score estimation, balance assessment, and effect
estimation. Fortunately, doing so is not complicated, but some care must
be taken to ensure sampling weights are incorporated correctly. It is
assumed one has read the other vignettes explaining matching
([`vignette("matching-methods")`](https://kosukeimai.github.io/MatchIt/articles/matching-methods.md)),
balance assessment
([`vignette("assessing-balance")`](https://kosukeimai.github.io/MatchIt/articles/assessing-balance.md)),
and effect estimation
([`vignette("estimating-effects")`](https://kosukeimai.github.io/MatchIt/articles/estimating-effects.md).

We will use the same simulated toy dataset used in
[`vignette("estimating-effects")`](https://kosukeimai.github.io/MatchIt/articles/estimating-effects.md)
except with the addition of a sampling weights variable, `SW`, which is
used to generalize the sample to a specific target population with a
distribution of covariates different from that of the sample. Code to
generate the covariates, treatment, and outcome is at the bottom of
[`vignette("estimating-effects")`](https://kosukeimai.github.io/MatchIt/articles/estimating-effects.md)
and code to generate the sampling weights is at the end of this
document. We will consider the effect of binary treatment `A` on
continuous outcome `Y_C`, adjusting for confounders `X1`-`X9`.

``` r
head(d)
```

    ##   A      X1      X2      X3       X4 X5      X6      X7      X8       X9     Y_C     SW
    ## 1 0  0.1725 -1.4283 -0.4103 -2.36059  1 -1.1199  0.6398 -0.4840 -0.59385 -3.5907  1.675
    ## 2 0 -1.0959  0.8463  0.2456 -0.12333  1 -2.2687 -1.4491 -0.5514 -0.31439 -1.5481  1.411
    ## 3 0  0.1768  0.7905 -0.8436  0.82366  1 -0.2221  0.2971 -0.6966 -0.69516  6.0714  2.332
    ## 4 0 -0.4595  0.1726  1.9542 -0.62661  1 -0.4019 -0.8294 -0.5384  0.20729  2.4906  1.644
    ## 5 1  0.3563 -1.8121  0.8135 -0.67189  1 -0.8297  1.7297 -0.6439 -0.02648 -0.6687  2.722
    ## 6 0 -2.4313 -1.7984 -1.2940  0.04609  1 -1.2419 -1.1252 -1.8659 -0.56513 -9.8504 14.773

``` r
library("MatchIt")
```

## Matching

When using sampling weights with propensity score matching, one has the
option of including the sampling weights in the model used to estimate
the propensity scores. Although evidence is mixed on whether this is
required ([Austin, Jembere, and Chiu 2016](#ref-austin2016); [Lenis et
al. 2019](#ref-lenis2019)), it can be a good idea. The choice should
depend on whether including the sampling weights improves the quality of
the matches. Specifications including and excluding sampling weights
should be tried to determine which is preferred.

To supply sampling weights to the propensity score-estimating function
in
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md),
the sampling weights variable should be supplied to the `s.weights`
argument. It can be supplied either as a numerical vector containing the
sampling weights, or a string or one-sided formula with the name of the
sampling weights variable in the supplied dataset. Below we demonstrate
including sampling weights into propensity scores estimated using
logistic regression for optimal full matching for the average treatment
effect in the population (ATE) (note that all methods and steps apply
the same way to all forms of matching and all estimands).

``` r
mF_s <- matchit(A ~ X1 + X2 + X3 + X4 + X5 + 
                  X6 + X7 + X8 + X9, data = d,
                method = "full", distance = "glm",
                estimand = "ATE", s.weights = ~SW)
mF_s
```

Notice that the description of the matching specification when the
`matchit` object is printed includes lines indicating that the sampling
weights were included in the estimation of the propensity score and that
they are present in the `matchit` object. It is stored in the
`s.weights` component of the `matchit` object. Note that at this stage,
the matching weights (stored in the `weights` component of the `matchit`
object) do not incorporate the sampling weights; they are calculated
simply as a result of the matching.

Now let’s perform full matching on a propensity score that does not
include the sampling weights in its estimation. Here we use the same
specification as was used in
[`vignette("estimating-effects")`](https://kosukeimai.github.io/MatchIt/articles/estimating-effects.md).

``` r
mF <- matchit(A ~ X1 + X2 + X3 + X4 + X5 + 
                X6 + X7 + X8 + X9, data = d,
              method = "full", distance = "glm",
              estimand = "ATE")
mF
```

Notice that there is no mention of sampling weights in the description
of the matching specification. However, to properly assess balance and
estimate effects, we need the sampling weights to be included in the
`matchit` object, even if they were not used at all in the matching. To
do so, we use the function
[`add_s.weights()`](https://kosukeimai.github.io/MatchIt/reference/add_s.weights.md),
which adds sampling weights to the supplied `matchit` objects.

``` r
mF <- add_s.weights(mF, ~SW)

mF
```

Now when we print the `matchit` object, we can see lines have been added
identifying that sampling weights are present but they were not used in
the estimation of the propensity score used in the matching.

Note that not all methods can involve sampling weights in the
estimation. Only methods that use the propensity score will be affected
by sampling weights; coarsened exact matching or Mahalanobis distance
optimal pair matching, for example, ignore the sampling weights, and
some propensity score estimation methods, like `randomForest` and `bart`
(as presently implemented), cannot incorporate sampling weights.
Sampling weights should still be supplied to
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
even when using these methods to avoid having to use
[`add_s.weights()`](https://kosukeimai.github.io/MatchIt/reference/add_s.weights.md)
and remembering which methods do or do not involve sampling weights.

## Assessing Balance

Now we need to decide which matching specification is the best to use
for effect estimation. We do this by selecting the one that yields the
best balance without sacrificing remaining effective sample size.
Because the sampling weights are incorporated into the `matchit` object,
the balance assessment tools in
[`plot.matchit()`](https://kosukeimai.github.io/MatchIt/reference/plot.matchit.md)
and
[`summary.matchit()`](https://kosukeimai.github.io/MatchIt/reference/summary.matchit.md)
incorporate them into their output.

We’ll use [`summary()`](https://rdrr.io/r/base/summary.html) to examine
balance on the two matching specifications. With sampling weights
included, the balance statistics for the unmatched data are weighted by
the sampling weights. The balance statistics for the matched data are
weighted by the product of the sampling weights and the matching
weights. It is the product of these weights that will be used in
estimating the treatment effect. Below we use
[`summary()`](https://rdrr.io/r/base/summary.html) to display balance
for the two matching specifications. No additional arguments to
[`summary()`](https://rdrr.io/r/base/summary.html) are required for it
to use the sampling weights; as long as they are in the `matchit` object
(either due to being supplied with the `s.weights` argument in the call
to
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
or to being added afterward by
[`add_s.weights()`](https://kosukeimai.github.io/MatchIt/reference/add_s.weights.md)),
they will be correctly incorporated into the balance statistics.

``` r
#Balance before matching and for the SW propensity score full matching
summary(mF_s)

#Balance for the non-SW propensity score full matching
summary(mF, un = FALSE)
```

The results of the two matching specifications are similar. Balance
appears to be slightly better when using the sampling weight-estimated
propensity scores than when using the unweighted propensity scores.
However, the effective sample size for the control group is larger when
using the unweighted propensity scores. Neither propensity score
specification achieves excellent balance, and more fiddling with the
matching specification (e.g., by changing the method of estimating
propensity scores, the type of matching, or the options used with the
matching) might yield a better matched set. For the purposes of this
analysis, we will move forward with the matching that used the sampling
weight-estimated propensity scores (`mF_s`) because of its superior
balance. Some of the remaining imbalance may be eliminated by adjusting
for the covariates in the outcome model.

Note that had we not added sampling weights to `mF`, the matching
specification that did not include the sampling weights, our balance
assessment would be inaccurate because the balance statistics would not
include the sampling weights. In this case, in fact, assessing balance
on `mF` without incorporated the sampling weights would have yielded
radically different results and a different conclusion. It is critical
to incorporate sampling weights into the `matchit` object using
[`add_s.weights()`](https://kosukeimai.github.io/MatchIt/reference/add_s.weights.md)
even if they are not included in the propensity score estimation.

## Estimating the Effect

Estimating the treatment effect after matching is straightforward when
using sampling weights. Effects are estimated in the same way as when
sampling weights are excluded, except that the matching weights must be
multiplied by the sampling weights for use in the outcome model to yield
accurate, generalizable estimates.
[`match_data()`](https://kosukeimai.github.io/MatchIt/reference/match_data.md)
and
[`get_matches()`](https://kosukeimai.github.io/MatchIt/reference/match_data.md)
do this automatically, so the weights produced by these functions
already are a product of the matching weights and the sampling weights.
Note this will only be true if sampling weights are incorporated into
the `matchit` object. With
[`avg_comparisons()`](https://rdrr.io/pkg/marginaleffects/man/comparisons.html),
only the sampling weights should be included when estimating the
treatment effect.

Below we estimate the effect of `A` on `Y_C` in the matched and sampling
weighted sample, adjusting for the covariates to improve precision and
decrease bias.

``` r
md_F_s <- match_data(mF_s)

fit <- lm(Y_C ~ A * (X1 + X2 + X3 + X4 + X5 + 
                       X6 + X7 + X8 + X9), data = md_F_s,
          weights = weights)

library("marginaleffects")
avg_comparisons(fit,
                variables = "A",
                vcov = ~subclass,
                newdata = subset(A == 1),
                wts = "SW")
```

Note that
[`match_data()`](https://kosukeimai.github.io/MatchIt/reference/match_data.md)
and `get_weights()` have the option `include.s.weights`, which, when set
to `FALSE`, makes it so the returned weights do not incorporate the
sampling weights and are simply the matching weights. Because one might
to forget to multiply the two sets of weights together, it is easier to
just use the default of `include.s.weights = TRUE` and ignore the
sampling weights in the rest of the analysis (because they are already
included in the returned weights).

## Code to Generate Data used in Examples

``` r
#Generatng data similar to Austin (2009) for demonstrating 
#treatment effect estimation with sampling weights
gen_X <- function(n) {
  X <- matrix(rnorm(9 * n), nrow = n, ncol = 9)
  X[,5] <- as.numeric(X[,5] < .5)
  X
}

#~20% treated
gen_A <- function(X) {
  LP_A <- - 1.2 + log(2)*X[,1] - log(1.5)*X[,2] + log(2)*X[,4] - log(2.4)*X[,5] + 
    log(2)*X[,7] - log(1.5)*X[,8]
  P_A <- plogis(LP_A)
  rbinom(nrow(X), 1, P_A)
}

# Continuous outcome
gen_Y_C <- function(A, X) {
  2*A + 2*X[,1] + 2*X[,2] + 2*X[,3] + 1*X[,4] + 2*X[,5] + 1*X[,6] + rnorm(length(A), 0, 5)
}
#Conditional:
#  MD: 2
#Marginal:
#  MD: 2

gen_SW <- function(X) {
  e <- rbinom(nrow(X), 1, .3)
  1/plogis(log(1.4)*X[,2] + log(.7)*X[,4] + log(.9)*X[,6] + log(1.5)*X[,8] + log(.9)*e +
             -log(.5)*e*X[,2] + log(.6)*e*X[,4])
}

set.seed(19599)

n <- 2000
X <- gen_X(n)
A <- gen_A(X)
SW <- gen_SW(X)

Y_C <- gen_Y_C(A, X)

d <- data.frame(A, X, Y_C, SW)
```

## References

Austin, Peter C., Nathaniel Jembere, and Maria Chiu. 2016. “Propensity
Score Matching and Complex Surveys.” *Statistical Methods in Medical
Research* 27 (4): 1240–57. <https://doi.org/10.1177/0962280216658920>.

DuGoff, Eva H., Megan Schuler, and Elizabeth A. Stuart. 2014.
“Generalizing Observational Study Results: Applying Propensity Score
Methods to Complex Surveys.” *Health Services Research* 49 (1): 284–303.
<https://doi.org/10.1111/1475-6773.12090>.

Lenis, David, Trang Quynh Nguyen, Nianbo Dong, and Elizabeth A. Stuart.
2019. “It’s All about Balance: Propensity Score Matching in the Context
of Complex Survey Data.” *Biostatistics* 20 (1): 147–63.
<https://doi.org/10.1093/biostatistics/kxx063>.
