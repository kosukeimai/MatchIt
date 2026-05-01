# Estimating Effects After Matching

## Introduction

After assessing balance and deciding on a matching specification, it
comes time to estimate the effect of the treatment in the matched
sample. How the effect is estimated and interpreted depends on the
desired estimand and the type of model used (if any). In addition to
estimating effects, estimating the uncertainty of the effects is
critical in communicating them and assessing whether the observed effect
is compatible with there being no effect in the population. This guide
explains how to estimate effects after various forms of matching and
with various outcome types. There may be situations that are not covered
here for which additional methodological research may be required, but
some of the recommended methods here can be used to guide such
applications.

This guide is structured as follows: first, information on the concepts
related to effect and standard error (SE) estimation is presented below.
Then, instructions for how to estimate effects and SEs are described for
the standard case (matching for the ATT with a continuous outcome) and
some other common circumstances. Finally, recommendations for reporting
results and tips to avoid making common mistakes are presented.

### Identifying the estimand

Before an effect is estimated, the estimand must be specified and
clarified. Although some aspects of the estimand depend not only on how
the effect is estimated after matching but also on the matching method
itself, other aspects must be considered at the time of effect
estimation and interpretation. Here, we consider three aspects of the
estimand: the population the effect is meant to generalize to (the
target population), the effect measure, and whether the effect is
marginal or conditional.

**The target population.** Different matching methods allow you to
estimate effects that can generalize to different target populations.
The most common estimand in matching is the average treatment effect in
the treated (ATT), which is the average effect of treatment for those
who receive treatment. This estimand is estimable for matching methods
that do not change the treated units (i.e., by weighting or discarding
units) and is requested in
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
by setting `estimand = "ATT"` (which is the default). The average
treatment effect in the population (ATE) is the average effect of
treatment for the population from which the sample is a random sample.
This estimand is estimable only for methods that allow the ATE and
either do not discard units from the sample or explicit target full
sample balance, which in `MatchIt` is limited to full matching,
subclassification, and profile matching when setting `estimand = "ATE"`.
When treated units are discarded (e.g., through the use of common
support restrictions, calipers, cardinality matching, or \[coarsened\]
exact matching), the estimand corresponds to neither the population ATT
nor the population ATE, but rather to an average treatment effect in the
remaining matched sample (ATM), which may not correspond to any specific
target population. See Greifer and Stuart
([2021](#ref-greiferChoosingEstimandWhen2021)) for a discussion on the
substantive considerations involved when choosing the target population
of the estimand.

**Marginal and conditional effects.** A marginal effect is a comparison
between the expected potential outcome under treatment and the expected
potential outcome under control. This is the same quantity estimated in
randomized trials without blocking or covariate adjustment and is
particularly useful for quantifying the overall effect of a policy or
population-wide intervention. A conditional effect is the comparison
between the expected potential outcomes in the treatment groups within
strata. This is useful for identifying the effect of a treatment for an
individual patient or a subset of the population.

**Effect measures.** The outcome types we consider here are continuous,
with the effect measured by the mean difference; binary, with the effect
measured by the risk difference (RD), risk ratio (RR), or odds ratio
(OR); and time-to-event (i.e., survival), with the effect measured by
the hazard ratio (HR). The RR, OR, and HR are *noncollapsible* effect
measures, which means the marginal effect on that scale is not a
(possibly) weighted average of the conditional effects within strata,
even if the stratum-specific effects are of the same magnitude. For
these effect measures, it is critical to distinguish between marginal
and conditional effects because different statistical methods target
different types of effects. The mean difference and RD are *collapsible*
effect measures, so the same methods can be used to estimate marginal
and conditional effects.

Our primary focus will be on marginal effects, which are appropriate for
all effect measures, easily interpretable, and require few modeling
assumptions. The “Common Mistakes” section includes examples of commonly
used methods that estimate conditional rather than marginal effects and
should not be used when marginal effects are desired.

### G-computation

To estimate marginal effects, we use a method known as g-computation
([Snowden et al.
2011](#ref-snowdenImplementationGComputationSimulated2011)) or
regression estimation ([Schafer and Kang
2008](#ref-schaferAverageCausalEffects2008)). This involves first
specifying a model for the outcome as a function of the treatment and
covariates. Then, for each unit, we compute their predicted values of
the outcome setting their treatment status to treated, and then again
for control, leaving us with two predicted outcome values for each unit,
which are estimates of the potential outcomes under each treatment
level. We compute the mean of each of the estimated potential outcomes
across the entire sample, which leaves us with two average estimated
potential outcomes. Finally, the contrast of these average estimated
potential outcomes (e.g., their difference or ratio, depending on the
effect measure desired) is the estimate of the treatment effect.

When doing g-computation after matching, a few additional considerations
are required. First, when we take the average of the estimated potential
outcomes under each treatment level, this must be a weighted average
that incorporates the matching weights. Second, if we want to target the
ATT or ATC, we only estimate potential outcomes for the treated or
control group, respectively (though we still generate predicted values
under both treatment and control).

G-computation as a framework for estimating effects after matching has a
number of advantages over other approaches. It works the same regardless
of the form of the outcome model or type of outcome (e.g., whether a
linear model is used for a continuous outcome or a logistic model is
used for a binary outcome); the only difference might be how the average
expected potential outcomes are contrasted in the final step. In simple
cases, the estimated effect is numerically identical to effects
estimated using other methods; for example, if no covariates are
included in the outcome model, the g-computation estimate is equal to
the difference in means from a t-test or coefficient of the treatment in
a linear model for the outcome. There are analytic approximations to the
SEs of the g-computation estimate, and these SEs can incorporate
pair/subclass membership (described in more detail below).

For all these reasons, we use g-computation when possible for all effect
estimates, even if there are simpler methods that would yield the same
estimates. Using a single workflow (with some slight modifications
depending on the context; see below) facilitates implementing best
practices regardless of what choices a user makes.

### Modeling the Outcome

The goal of the outcome model is to generate good predictions for use in
the g-computation procedure described above. The type and form of the
outcome model should depend on the outcome type. For continuous
outcomes, one can use a linear model regressing the outcome on the
treatment; for binary outcomes, one can use a generalized linear model
with, e.g., a logistic link; for time-to-event outcomes, one can use a
Cox proportional hazards model.

An additional decision to make is whether (and how) to include
covariates in the outcome model. One may ask, why use matching at all if
you are going to model the outcome with covariates anyway? Matching
reduces the dependence of the effect estimate on correct specification
of the outcome model; this is the central thesis of Ho et al.
([2007](#ref-ho2007)). Including covariates in the outcome model after
matching has several functions: it can increase precision in the effect
estimate, reduce the bias due to residual imbalance, and make the effect
estimate “doubly robust”, which means it is consistent if either the
matching reduces sufficient imbalance in the covariates or if the
outcome model is correct. For these reasons, we recommend covariate
adjustment after matching when possible. There is some evidence that
covariate adjustment is most helpful for covariates with standardized
mean differences greater than .1 ([Nguyen et al.
2017](#ref-nguyen2017)), so these covariates and covariates thought to
be highly predictive of the outcome should be prioritized in treatment
effect models if not all can be included due to sample size constraints.

Although there are many possible ways to include covariates (e.g., not
just main effects but interactions, smoothing terms like splines, or
other nonlinear transformations), it is important not to engage in
specification search (i.e., trying many outcomes models in search of the
“best” one). Doing so can invalidate results and yield a conclusion that
fails to replicate. For this reason, we recommend only including the
same terms included in the propensity score model unless there is a
strong *a priori* and justifiable reason to model the outcome
differently.

It is important not to interpret the coefficients and tests of
covariates in the outcome model. These are not causal effects and their
estimates may be severely confounded. Only the treatment effect estimate
can be interpreted as causal assuming the relevant assumptions about
unconfoundedness are met. Inappropriately interpreting the coefficients
of covariates in the outcome model is known as the Table 2 fallacy
([Westreich and Greenland 2013](#ref-westreich2013)). To avoid this, we
only display the results of the g-computation procedure and do not
examine or interpret the outcome models themselves.

### Estimating Standard Errors and Confidence Intervals

Uncertainty estimation (i.e., of SEs, confidence intervals, and
p-values) may consider the variety of sources of uncertainty present in
the analysis, including (but not limited to!) estimation of the
propensity score (if used), matching (i.e., because treated units might
be matched to different control units if others had been sampled), and
estimation of the treatment effect (i.e., because of sampling error). In
general, there are no analytic solutions to all these issues, so much of
the research done on uncertainty estimation after matching has relied on
simulation studies. The two primary methods that have been shown to
perform well in matched samples are using cluster-robust SEs and the
bootstrap, described below.

To compute SEs after g-computation, a method known as the delta method
is used; this is a way to compute the SEs of the derived quantities (the
expected potential outcomes and their contrast) from the variance of the
coefficients of the outcome models. For nonlinear models (e.g., logistic
regression), the delta method is only an approximation subject to error
(though in many cases this error is small and shrinks in large samples).
Because the delta method relies on the variance of the coefficients from
the outcome model, it is important to correctly estimate these
variances, using either robust or cluster-robust methods as described
below.

#### Robust and Cluster-Robust Standard Errors

**Robust standard errors.** Also known as sandwich SEs (due to the form
of the formula for computing them), heteroscedasticity-consistent SEs,
or Huber-White SEs, robust SEs are an adjustment to the usual maximum
likelihood or ordinary least squares SEs that are robust to violations
of some of the assumptions required for usual SEs to be valid
([MacKinnon and White 1985](#ref-mackinnon1985)). Although there has
been some debate about their utility ([King and Roberts
2015](#ref-king2015)), robust SEs rarely degrade inferences and often
improve them. Generally, robust SEs **must** be used when any
non-uniform weights are included in the estimation (e.g., with matching
with replacement or inverse probability weighting).

**Cluster-robust standard errors.** A version of robust SEs known as
cluster-robust SEs ([Liang and Zeger 1986](#ref-liang1986)) can be used
to account for dependence between observations within clusters (e.g.,
matched pairs). Abadie and Spiess ([2019](#ref-abadie2019)) demonstrate
analytically that cluster-robust SEs are generally valid after matching,
whereas regular robust SEs can over- or under-estimate the true sampling
variability of the effect estimator depending on the specification of
the outcome model (if any) and degree of effect modification. A plethora
of simulation studies have further confirmed the validity of
cluster-robust SEs after matching (e.g., [Austin
2009](#ref-austin2009a), [2013a](#ref-austin2013); [Austin and Small
2014](#ref-austin2014); [Gayat et al. 2012](#ref-gayat2012); [Wan
2019](#ref-wan2019)). Given this evidence favoring the use of
cluster-robust SEs, we recommend them in most cases and use them
judiciously in this guide[^1].

#### Bootstrapping

One problem when using robust and cluster-robust SEs along with the
delta method is that the delta method is an approximation, as previously
mentioned. One solution to this problem is bootstrapping, which is a
technique used to simulate the sampling distribution of an estimator by
repeatedly drawing samples with replacement and estimating the effect in
each bootstrap sample ([Efron and Tibshirani 1993](#ref-efron1993)).
From the bootstrap distribution, SEs and confidence intervals can be
computed in several ways, including using the standard deviation of the
bootstrap estimates as the SE estimate or using the 2.5 and 97.5
percentiles as 95% confidence interval bounds. Bootstrapping tends to be
most useful when no analytic estimator of a SE is possible or has been
derived yet. Although Abadie and Imbens ([2008](#ref-abadie2008)) found
analytically that the bootstrap is inappropriate for matched samples,
simulation evidence has found it to be adequate in many cases ([Hill and
Reiter 2006](#ref-hill2006); [Austin and Small 2014](#ref-austin2014);
[Austin and Stuart 2017](#ref-austin2017)).

Typically, bootstrapping involves performing the entire estimation
process in each bootstrap sample, including propensity score estimation,
matching, and effect estimation. This tends to be the most
straightforward route, though intervals from this method may be
conservative in some cases (i.e., they are wider than necessary to
achieve nominal coverage) ([Austin and Small 2014](#ref-austin2014)).
Less conservative and more accurate intervals have been found when using
different forms of the bootstrap, including the wild bootstrap develop
by Bodory et al. ([2020](#ref-bodory2020)) and the matched/cluster
bootstrap described by Austin and Small ([2014](#ref-austin2014)) and
Abadie and Spiess ([2019](#ref-abadie2019)). The cluster bootstrap
involves sampling matched pairs/strata of units from the matched sample
and performing the analysis within each sample composed of the sampled
pairs. Abadie and Spiess ([2019](#ref-abadie2019)) derived analytically
that the cluster bootstrap is valid for estimating SEs and confidence
intervals in the same circumstances cluster robust SEs are; indeed, the
cluster bootstrap SE is known to approximate the cluster-robust SE
([Cameron and Miller 2015](#ref-cameron2015)).

With bootstrapping, more bootstrap replications are always better but
can take time and increase the chances that at least one error will
occur within the bootstrap analysis (e.g., a bootstrap sample with zero
treated units or zero units with an event). In general, numbers of
replications upwards of 999 are recommended, with values one less than a
multiple of 100 preferred to avoid interpolation when using the
percentiles as confidence interval limits ([MacKinnon
2006](#ref-mackinnon2006)). There are several methods of computing
bootstrap confidence intervals, but the bias-corrected accelerated (BCa)
bootstrap confidence interval often performs best ([Austin and Small
2014](#ref-austin2014); [Carpenter and Bithell
2000](#ref-carpenter2000)) and is easy to implement, simply by setting
`type = "bca"` in the call to
[`boot::boot.ci()`](https://rdrr.io/pkg/boot/man/boot.ci.html) after
running [`boot::boot()`](https://rdrr.io/pkg/boot/man/boot.html)[^2].

Most of this guide will consider analytic (i.e., non-bootstrapping)
approaches to estimating uncertainty; the section “Using Bootstrapping
to Estimate Confidence Intervals” describes broadly how to use
bootstrapping. Although analytic estimates are faster to compute, in
many cases bootstrap confidence intervals are more accurate.

## Estimating Treatment Effects and Standard Errors After Matching

Below, we describe effect estimation after matching. We’ll be using a
simulated toy dataset `d` with several outcome types. Code to generate
the dataset is at the end of this document. The focus here is not on
evaluating the methods but simply on demonstrating them. In all cases,
the correct propensity score model is used. Below we display the first
six rows of `d`:

``` r

head(d)
```

    ##   A      X1      X2      X3       X4 X5      X6      X7      X8       X9      Y_C Y_B     Y_S
    ## 1 0  0.1725 -1.4283 -0.4103 -2.36059  1 -1.1199  0.6398 -0.4840 -0.59385  0.07104   0  278.46
    ## 2 0 -1.0959  0.8463  0.2456 -0.12333  1 -2.2687 -1.4491 -0.5514 -0.31439  0.15619   0  330.63
    ## 3 0  0.1768  0.7905 -0.8436  0.82366  1 -0.2221  0.2971 -0.6966 -0.69516 -0.85180   1  369.94
    ## 4 0 -0.4595  0.1726  1.9542 -0.62661  1 -0.4019 -0.8294 -0.5384  0.20729 -2.35184   0   91.06
    ## 5 1  0.3563 -1.8121  0.8135 -0.67189  1 -0.8297  1.7297 -0.6439 -0.02648  0.68058   0  182.73
    ## 6 0 -2.4313 -1.7984 -1.2940  0.04609  1 -1.2419 -1.1252 -1.8659 -0.56513 -5.62260   0 2563.73

`A` is the treatment variable, `X1` through `X9` are covariates, `Y_C`
is a continuous outcome, `Y_B` is a binary outcome, and `Y_S` is a
survival outcome.

We will need to the following packages to perform the desired analyses:

- `marginaleffects` provides the `avg_comparisons()` function for
  performing g-computation and estimating the SEs and confidence
  intervals of the average estimate potential outcomes and treatment
  effects
- `sandwich` is used internally by `marginaleffects` to compute robust
  and cluster-robust SEs
- `survival` provides
  [`coxph()`](https://rdrr.io/pkg/survival/man/coxph.html) to estimate
  the coefficients in a Cox-proportional hazards model for the marginal
  hazard ratio, which we will use for survival outcomes.

Of course, we also need `MatchIt` to perform the matching.

``` r

library("MatchIt")
## library("marginaleffects")
```

All effect estimates will be computed using
`marginaleffects::avg_comparions()`, even when its use may be
superfluous (e.g., for performing a t-test in the matched set). As
previously mentioned, this is because it is useful to have a single
workflow that works no matter the situation, perhaps with very slight
modifications to accommodate different contexts. Using
`avg_comparions()` has several advantages, even when the alternatives
are simple: it only provides the effect estimate, and not other
coefficients; it automatically incorporates robust and cluster-robust
SEs if requested; and it always produces average marginal effects for
the correct population if requested.

Other packages may be of use but are not used here. There are
alternatives to the `marginaleffects` package for computing average
marginal effects, including `margins` and `stdReg`. The `survey` package
can be used to estimate robust SEs incorporating weights and provides
functions for survey-weighted generalized linear models and
Cox-proportional hazards models.

### The Standard Case

For almost all matching methods, whether a caliper, common support
restriction, exact matching specification, or $`k`$:1 matching
specification is used, estimating the effect in the matched dataset is
straightforward and involves fitting a model for the outcome that
incorporates the matching weights[^3], then estimating the treatment
effect using g-computation (i.e., using
[`marginaleffects::avg_comparisons()`](https://rdrr.io/pkg/marginaleffects/man/comparisons.html))
with a cluster-robust SE to account for pair membership. This procedure
is the same for continuous and binary outcomes with and without
covariates.

There are a few adjustments that need to be made for certain scenarios,
which we describe in the section “Adjustments to the Standard Case”.
These adjustments include for the following cases: when matching for the
ATE rather than the ATT, for matching with replacement, for matching
with a method that doesn’t involve creating pairs (e.g., cardinality and
profile matching and coarsened exact matching), for subclassification,
for estimating effects with binary outcomes, and for estimating effects
with survival outcomes. You must read the Standard Case to understand
the basic procedure before reading about these special scenarios.

Here, we demonstrate the faster analytic approach to estimating
confidence intervals; for the bootstrap approach, see the section “Using
Bootstrapping to Estimate Confidence Intervals” below.

First, we will perform variable-ratio nearest neighbor matching without
replacement on the propensity score for the ATT. Remember, all matching
methods use this exact procedure or a slight variation, so this section
is critical even if you are using a different matching method.

``` r

#Variable-ratio NN matching on the PS for the ATT
mV <- matchit(A ~ X1 + X2 + X3 + X4 + X5 + 
                X6 + X7 + X8 + X9,
              data = d,
              ratio = 2,
              max.controls = 4)
mV
```

    ## A `matchit` object
    ##  - method: Variable ratio 2:1 nearest neighbor matching without replacement
    ##  - distance: Propensity score
    ##              - estimated with logistic regression
    ##  - number of obs.: 2000 (original), 1323 (matched)
    ##  - target estimand: ATT
    ##  - covariates: X1, X2, X3, X4, X5, X6, X7, X8, X9

``` r

#Extract matched data
md <- match_data(mV)

head(md)
```

    ##    A      X1      X2      X3      X4 X5      X6      X7      X8       X9      Y_C Y_B    Y_S distance weights subclass
    ## 1  0  0.1725 -1.4283 -0.4103 -2.3606  1 -1.1199  0.6398 -0.4840 -0.59385  0.07104   0 278.46  0.08461     0.5      365
    ## 3  0  0.1768  0.7905 -0.8436  0.8237  1 -0.2221  0.2971 -0.6966 -0.69516 -0.85180   1 369.94  0.22210     0.5       42
    ## 5  1  0.3563 -1.8121  0.8135 -0.6719  1 -0.8297  1.7297 -0.6439 -0.02648  0.68058   0 182.73  0.43291     1.0        1
    ## 7  0  1.8402  1.7601 -1.0746 -1.6428  1  1.4482  0.7131  0.6972 -0.94673  4.28651   1  97.49  0.09274     0.5        6
    ## 9  0  0.7808  1.3137  0.6580  0.8540  1  0.9495 -0.5731 -0.2362 -0.14580 15.89771   1  67.53  0.15751     0.5      218
    ## 10 1 -0.5651 -0.1053 -0.1369  1.6233  1 -0.5304 -0.3342  0.4184  0.46308  1.07888   1 113.70  0.16697     1.0        2

Typically one would assess balance and ensure that this matching
specification works, but we will skip that step here to focus on effect
estimation. See
[`vignette("MatchIt")`](https://kosukeimai.github.io/MatchIt/articles/MatchIt.md)
and
[`vignette("assessing-balance")`](https://kosukeimai.github.io/MatchIt/articles/assessing-balance.md)
for more information on this necessary step. Because we did not use a
caliper, the target estimand is the ATT.

We perform all analyses using the matched dataset, `md`, which, for
matching methods that involve dropping units, contains only the units
retained in the sample.

First, we fit a model for the outcome given the treatment and
(optionally) the covariates. It’s usually a good idea to include
treatment-covariate interactions, which we do below, but this is not
always necessary, especially when excellent balance has been achieved.
You can also include the propensity score (usually labeled `distance` in
the
[`match_data()`](https://kosukeimai.github.io/MatchIt/reference/match_data.md)
output), which can add some robustness, especially when modeled flexibly
(e.g., with polynomial terms or splines) ([Austin
2017](#ref-austinDoublePropensityscoreAdjustment2017)); see
[here](https://stats.stackexchange.com/a/580174/116195) for an example.

``` r

#Linear model with covariates
fit1 <- lm(Y_C ~ A * (X1 + X2 + X3 + X4 + X5 + 
                        X6 + X7 + X8 + X9),
           data = md,
           weights = weights)
```

Next, we use
[`marginaleffects::avg_comparisons()`](https://rdrr.io/pkg/marginaleffects/man/comparisons.html)
to estimate the ATT.

``` r

avg_comparisons(fit1,
                variables = "A",
                vcov = ~subclass,
                newdata = subset(A == 1))
```

Let’s break down the call to `avg_comparisons()`: to the first argument,
we supply the model fit, `fit1`; to the `variables` argument, the name
of the treatment (`"A"`); to the `vcov` argument, a formula with
subclass membership (`~subclass`) to request cluster-robust SEs; and to
the `newdata` argument, a version of the matched dataset containing only
the treated units (`subset(A == 1)`) to request the ATT. Some of these
arguments differ depending on the specifics of the matching method and
outcome type; see the sections below for information.

If, in addition to the effect estimate, we want the average estimated
potential outcomes, we can use
[`marginaleffects::avg_predictions()`](https://rdrr.io/pkg/marginaleffects/man/predictions.html),
which we demonstrate below. Note the interpretation of the resulting
estimates as the expected potential outcomes is only valid if all
covariates present in the outcome model (if any) are interacted with the
treatment.

``` r

avg_predictions(fit1,
                variables = "A",
                vcov = ~subclass,
                newdata = subset(A == 1))
```

We can see that the difference in potential outcome means is equal to
the average treatment effect computed previously[^4]. All of the
arguments to `avg_predictions()` are the same as those to
`avg_comparisons()`.

### Adjustments to the Standard Case

This section explains how the procedure might differ if any of the
following special circumstances occur.

#### Matching for the ATE

When matching for the ATE (including \[coarsened\] exact matching, full
matching, subclassification, and cardinality matching), everything is
identical to the Standard Case except that in the calls to
`avg_comparisons()` and `avg_predictions()`, the `newdata` argument is
omitted. This is because the estimated potential outcomes are computed
for the full sample rather than just the treated units.

#### Matching with replacement

When matching with replacement (i.e., nearest neighbor or genetic
matching with `replace = TRUE`), effect and SE estimation need to
account for control unit multiplicity (i.e., repeated use) and
within-pair correlations ([Hill and Reiter 2006](#ref-hill2006); [Austin
and Cafri 2020](#ref-austin2020a)). Although Abadie and Imbens
([2008](#ref-abadie2008)) demonstrated analytically that bootstrap SEs
may be invalid for matching with replacement, simulation work by Hill
and Reiter ([2006](#ref-hill2006)) and Bodory et al.
([2020](#ref-bodory2020)) has found that bootstrap SEs are adequate and
generally slightly conservative. See the section “Using Bootstrapping to
Estimate Confidence Intervals” for instructions on using the bootstrap
and an example that use matching with replacement.

Because control units do not belong to unique pairs, there is no pair
membership in the
[`match_data()`](https://kosukeimai.github.io/MatchIt/reference/match_data.md)
output. One can simply change `vcov = ~subclass` to `vcov = "HC3"` in
the calls to `avg_comparisons()` and `avg_predictions()` to use robust
SEs instead of cluster-robust SEs, as recommended by Hill and Reiter
([2006](#ref-hill2006)). There is some evidence for an alternative
approach that incorporates pair membership and adjusts for reuse of
control units, though this has only been studied for survival outcomes
([Austin and Cafri 2020](#ref-austin2020a)). This adjustment involves
using two-way cluster-robust SEs with pair membership and unit ID as the
clustering variables. For continuous and binary outcomes, this involves
the following two changes: 1) replace
[`match_data()`](https://kosukeimai.github.io/MatchIt/reference/match_data.md)
with
[`get_matches()`](https://kosukeimai.github.io/MatchIt/reference/match_data.md),
which produces a dataset with one row per unit per pair, meaning control
units matched to multiple treated units will appear multiple times in
the dataset; 2) set `vcov = ~subclass + id` in the calls to
`avg_comparisons()` and `avg_predictions()`. For survival outcomes, a
special procedure must be used; see the section on survival outcomes
below.

#### Matching without pairing

Some matching methods do not involve creating pairs; these include
cardinality and profile matching with `mahvars = NULL` (the default),
exact matching, and coarsened exact matching with `k2k = FALSE` (the
default). The only change that needs to be made to the Standard Case is
that one should change `vcov = ~subclass` to `vcov = "HC3"` in the calls
to `avg_comparisons()` and `avg_predictions()` to use robust SEs instead
of cluster-robust SEs. Remember that if matching is done for the ATE
(even if units are dropped), the `newdata` argument should be dropped.

#### Propensity score subclassification

There are two natural ways to estimate marginal effects after
subclassification: the first is to estimate subclass-specific treatment
effects and pool them using an average marginal effects procedure, and
the second is to use the stratum weights to estimate a single average
marginal effect. This latter approach is also known as marginal mean
weighting through stratification (MMWS), and is described in detail by
Hong ([2010](#ref-hong2010))[^5]. When done properly, both methods
should yield similar or identical estimates of the treatment effect.

All of the methods described above for the Standard Case also work with
MMWS because the formation of the weights is the same; the only
difference is that it is not appropriate to use cluster-robust SEs with
MMWS because of how few clusters are present, so one should change
`vcov = ~subclass` to `vcov = "HC3"` in the calls to `avg_comparisons()`
and `avg_predictions()` to use robust SEs instead of cluster-robust SEs.
The subclasses can optionally be included in the outcome model
(optionally interacting with treatment) as an alternative to including
the propensity score.

The subclass-specific approach omits the weights and uses the subclasses
directly. It is only appropriate when there are a small number of
subclasses relative to the sample size. In the outcome model, `subclass`
should interact with all other predictors in the model (including the
treatment, covariates, and interactions, if any), and the `weights`
argument should be omitted. As with MMWS, one should change
`vcov = ~subclass` to `vcov = "HC3"` in the calls to `avg_comparisons()`
and `avg_predictions()`. See an example below:

``` r

#Subclassification on the PS for the ATT
mS <- matchit(A ~ X1 + X2 + X3 + X4 + X5 + 
                X6 + X7 + X8 + X9,
              data = d,
              method = "subclass",
              estimand = "ATT")

#Extract matched data
md <- match_data(mS)

fitS <- lm(Y_C ~ subclass * (A * (X1 + X2 + X3 + X4 + X5 + 
                                    X6 + X7 + X8 + X9)),
           data = md)

avg_comparisons(fitS,
                variables = "A",
                vcov = "HC3",
                newdata = subset(A == 1))
```

A model with fewer terms may be required when subclasses are small;
removing covariates or their interactions with treatment may be required
and can increase precision in smaller datasets. Remember that if
subclassification is done for the ATE (even if units are dropped), the
`newdata` argument should be dropped.

#### Binary outcomes

Estimating effects on binary outcomes is essentially the same as for
continuous outcomes. The main difference is that there are several
measures of the effect one can consider, which include the odds ratio
(OR), risk ratio/relative risk (RR), and risk difference (RD), and the
syntax to `avg_comparisons()` depends on which one is desired. The
outcome model should be one appropriate for binary outcomes (e.g.,
logistic regression) but is unrelated to the desired effect measure
because we can compute any of the above effect measures using
`avg_comparisons()` after the logistic regression.

To fit a logistic regression model, change
[`lm()`](https://rdrr.io/r/stats/lm.html) to
[`glm()`](https://rdrr.io/r/stats/glm.html) and set
`family = quasibinomial()`[^6]. To compute the marginal RD, we can use
exactly the same syntax as in the Standard Case; nothing needs to
change[^7].

To compute the marginal RR, we need to add `comparison = "lnratioavg"`
to `avg_comparisons()`; this computes the marginal log RR. To get the
marginal RR, we need to add `transform = "exp"` to `avg_comparisons()`,
which exponentiates the marginal log RR and its confidence interval. The
code below computes the effects and displays the statistics of interest:

``` r

#Logistic regression model with covariates
fit2 <- glm(Y_B ~ A * (X1 + X2 + X3 + X4 + X5 + 
                         X6 + X7 + X8 + X9),
            data = md,
            weights = weights,
            family = quasibinomial())

#Compute effects; RR and confidence interval
avg_comparisons(fit2,
                variables = "A",
                vcov = ~subclass,
                newdata = subset(A == 1),
                comparison = "lnratioavg",
                transform = "exp")
```

The output displays the marginal RR, its Z-value, the p-value for the
Z-test of the log RR against 0, and its confidence interval. (Note that
even though the `Contrast` label still suggests the log RR, the RR is
actually displayed.) To view the log RR and its standard error, omit the
`transform` argument.

For the marginal OR, the only thing that needs to change is that
`comparison` should be set to `"lnoravg"`. For the marginal RD, both the
`comparison` and `transform` arguments can be removed (yielding the same
call as in the standard case).

#### Survival outcomes

There are several measures of effect size for survival outcomes. When
using the Cox proportional hazards model, the quantity of interest is
the hazard ratio (HR) between the treated and control groups. As with
the OR, the HR is non-collapsible, which means the estimated HR will
only be a valid estimate of the marginal HR when no other covariates are
included in the model. Other effect measures, such as the difference in
mean survival times or probability of survival after a given time, can
be treated just like continuous and binary outcomes as previously
described.

For the HR, we cannot compute average marginal effects and must use the
coefficient on treatment in a Cox model fit without covariates[^8]. This
means that we cannot use the procedures from the Standard Case. Here we
describe estimating the marginal HR using
[`coxph()`](https://rdrr.io/pkg/survival/man/coxph.html) from the
`survival` package. (See
[`help("coxph", package = "survival")`](https://rdrr.io/pkg/survival/man/coxph.html)
for more information on this model.) To request cluster-robust SEs as
recommended by Austin ([2013b](#ref-austin2013a)), we need to supply
pair membership (stored in the `subclass` column of `md`) to the
`cluster` argument and set `robust = TRUE`. For matching methods that
don’t involve pairing (e.g., cardinality and profile matching and
\[coarsened\] exact matching), we can omit the `cluster` argument (but
keep `robust = TRUE`)[^9].

``` r

library("survival")

#Cox Regression for marginal HR
coxph(Surv(Y_S) ~ A,
      data = md,
      robust = TRUE, 
      weights = weights,
      cluster = subclass)
```

    ## Call:
    ## coxph(formula = Surv(Y_S) ~ A, data = md, weights = weights, 
    ##     robust = TRUE, cluster = subclass)
    ## 
    ##   coef exp(coef) se(coef) robust se z     p
    ## A 0.47      1.60     0.06      0.07 7 2e-12
    ## 
    ## Likelihood ratio test=61  on 1 df, p=7e-15
    ## n= 1323, number of events= 1323

The `coef` column contains the log HR, and `exp(coef)` contains the HR.
Remember to always use the `robust se` for the SE of the log HR. The
displayed z-test p-value results from using the robust SE.

For matching with replacement, a special procedure described by Austin
and Cafri ([2020](#ref-austin2020a)) can be necessary for valid
inference. According to the results of their simulation studies, when
the treatment prevalence is low (\<30%), a SE that does not involve pair
membership (i.e., the
[`match_data()`](https://kosukeimai.github.io/MatchIt/reference/match_data.md)
approach, as demonstrated above) is sufficient. When treatment
prevalence is higher, the SE that ignores pair membership may be too
low, and the authors recommend using a custom SE estimator that uses
information about both multiplicity and pairing.

Doing so must be done manually for survival models using
[`get_matches()`](https://kosukeimai.github.io/MatchIt/reference/match_data.md)
and several calls to
[`coxph()`](https://rdrr.io/pkg/survival/man/coxph.html) as demonstrated
in the appendix of Austin and Cafri ([2020](#ref-austin2020a)). We
demonstrate this below:

``` r

#get_matches() after matching with replacement
gm <- get_matches(mR)

#Austin & Cafri's (2020) SE estimator
fs <- coxph(Surv(Y_S) ~ A, data = gm, robust = TRUE, 
            weights = weights, cluster = subclass)
Vs <- fs$var
ks <- nlevels(gm$subclass)

fi <- coxph(Surv(Y_S) ~ A, data = gm, robust = TRUE, 
            weights = weights, cluster = id)
Vi <- fi$var
ki <- length(unique(gm$id))

fc <- coxph(Surv(Y_S) ~ A, data = gm, robust = TRUE, 
            weights = weights)
Vc <- fc$var
kc <- nrow(gm)

#Compute the variance and sneak it back into the fit object
fc$var <- (ks/(ks-1))*Vs + (ki/(ki-1))*Vi - (kc/(kc-1))*Vc

fc
```

The `robust se` column contains the computed SE, and the reported Z-test
uses this SE. The `se(coef)` column should be ignored.

### Using Bootstrapping to Estimate Confidence Intervals

The bootstrap is an alternative to the delta method for estimating
confidence intervals for estimated effects. See the section
Bootstrapping above for details. Here, we’ll demonstrate two forms of
the bootstrap: 1) the standard bootstrap, which involve resampling units
and performing matching and effect estimation within each bootstrap
sample, and 2) the cluster bootstrap, which involves resampling pairs
after matching and estimating the effect in each bootstrap sample. For
both, we will use functionality in the `boot` package. It is critical to
set a seed using [`set.seed()`](https://rdrr.io/r/base/Random.html)
prior to performing the bootstrap in order for results to be replicable.

#### The standard bootstrap

For the standard bootstrap, we need a function that takes in the
original dataset and a vector of sampled unit indices and returns the
estimated quantity of interest. This function should perform the
matching on the bootstrap sample, fit the outcome model, and estimate
the treatment effect using g-computation. In this example, we’ll use
matching with replacement, since the standard bootstrap has been found
to work well with it ([Bodory et al. 2020](#ref-bodory2020); [Hill and
Reiter 2006](#ref-hill2006)), despite some analytic results recommending
otherwise ([Abadie and Imbens 2008](#ref-abadie2008)). We’ll implement
g-computation manually rather than using `avg_comparisons()`, as this
dramatically improves the speed of the estimation since we don’t require
standard errors to be estimated in each sample (or other processing
`avg_comparisons()` does). We’ll consider the marginal RR ATT of `A` on
the binary outcome `Y_B`.

The first step is to write the estimation function, we call `boot_fun`.
This function returns the marginal RR. In it, we perform the matching,
estimate the effect, and return the estimate of interest.

``` r

boot_fun <- function(data, i) {
  boot_data <- data[i,]
  
  #Do 1:1 PS matching with replacement
  m <- matchit(A ~ X1 + X2 + X3 + X4 + X5 + 
                 X6 + X7 + X8 + X9,
               data = boot_data,
               replace = TRUE)
  
  #Extract matched dataset
  md <- match_data(m, data = boot_data)
  
  #Fit outcome model
  fit <- glm(Y_B ~ A * (X1 + X2 + X3 + X4 + X5 + 
                          X6 + X7 + X8 + X9),
             data = md, weights = weights,
             family = quasibinomial())
  
  ## G-computation ##
  #Subset to treated units for ATT; skip for ATE
  md1 <- subset(md, A == 1)
  
  #Estimated potential outcomes under treatment
  p1 <- predict(fit, type = "response",
                newdata = transform(md1, A = 1))
  Ep1 <- mean(p1)
  
  #Estimated potential outcomes under control
  p0 <- predict(fit, type = "response",
                newdata = transform(md1, A = 0))
  Ep0 <- mean(p0)
  
  #Risk ratio
  Ep1 / Ep0
}
```

Next, we call [`boot::boot()`](https://rdrr.io/pkg/boot/man/boot.html)
with this function and the original dataset supplied to perform the
bootstrapping. We’ll request 199 bootstrap replications here, but in
practice you should use many more, upwards of 999. More is always
better. Using more also allows you to use the bias-corrected and
accelerated (BCa) bootstrap confidence intervals (which you can request
by setting `type = "bca"` in the call to
[`boot.ci()`](https://rdrr.io/pkg/boot/man/boot.ci.html)), which are
known to be the most accurate. See
[`?boot.ci`](https://rdrr.io/pkg/boot/man/boot.ci.html) for details.
Here, we’ll just use a percentile confidence interval.

``` r

library("boot")
set.seed(54321)
boot_out <- boot(d, boot_fun, R = 199)

boot_out
```

    ## 
    ## ORDINARY NONPARAMETRIC BOOTSTRAP
    ## 
    ## 
    ## Call:
    ## boot(data = d, statistic = boot_fun, R = 199)
    ## 
    ## 
    ## Bootstrap Statistics :
    ##     original  bias    std. error
    ## t1*    1.347  0.1417      0.1937

``` r

boot.ci(boot_out, type = "perc")
```

    ## BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS
    ## Based on 199 bootstrap replicates
    ## 
    ## CALL : 
    ## boot.ci(boot.out = boot_out, type = "perc")
    ## 
    ## Intervals : 
    ## Level     Percentile     
    ## 95%   ( 1.144,  1.891 )  
    ## Calculations and Intervals on Original Scale
    ## Some percentile intervals may be unstable

We find a RR of 1.347 with a confidence interval of (1.144, 1.891). If
we had wanted a risk difference, we could have changed the final line in
`boot_fun()` to be `Ep1 - Ep0`.

#### The cluster bootstrap

For the cluster bootstrap, we need a function that takes in a vector of
subclass (e.g., pairs) and a vector of sampled pair indices and returns
the estimated quantity of interest. This function should fit the outcome
model and estimate the treatment effect using g-computation, but the
matching step occurs prior to the bootstrap. Here, we’ll use matching
without replacement, since the cluster bootstrap has been found to work
well with it ([Austin and Small 2014](#ref-austin2014); [Abadie and
Spiess 2019](#ref-abadie2019)). This could be used for any method that
returns pair membership, including other pair matching methods without
replacement and full matching.

As before, we’ll use g-computation to estimate the marginal RR ATT, and
we’ll do so manually rather than using `avg_comparisons()` for speed.
Note that the cluster bootstrap is already much faster than the standard
bootstrap because matching does not need to occur within each bootstrap
sample. First, we’ll do a round of matching.

``` r

mNN <- matchit(A ~ X1 + X2 + X3 + X4 + X5 + 
                 X6 + X7 + X8 + X9, data = d)
mNN
```

    ## A `matchit` object
    ##  - method: 1:1 nearest neighbor matching without replacement
    ##  - distance: Propensity score
    ##              - estimated with logistic regression
    ##  - number of obs.: 2000 (original), 882 (matched)
    ##  - target estimand: ATT
    ##  - covariates: X1, X2, X3, X4, X5, X6, X7, X8, X9

``` r

md <- match_data(mNN)
```

Next, we’ll write the function that takes in cluster membership and the
sampled indices and returns an estimate.

``` r

#Unique pair IDs
pair_ids <- levels(md$subclass)

#Unit IDs, split by pair membership
split_inds <- split(seq_len(nrow(md)), md$subclass)

cluster_boot_fun <- function(pairs, i) {
  
  #Extract units corresponding to selected pairs
  ids <- unlist(split_inds[pairs[i]])
  
  #Subset md with block bootstrapped indices
  boot_md <- md[ids,]
  
  #Fit outcome model
  fit <- glm(Y_B ~ A * (X1 + X2 + X3 + X4 + X5 + 
                          X6 + X7 + X8 + X9),
             data = boot_md, weights = weights,
             family = quasibinomial())
  
  ## G-computation ##
  #Subset to treated units for ATT; skip for ATE
  md1 <- subset(boot_md, A == 1)
  
  #Estimated potential outcomes under treatment
  p1 <- predict(fit, type = "response",
                newdata = transform(md1, A = 1))
  Ep1 <- mean(p1)
  
  #Estimated potential outcomes under control
  p0 <- predict(fit, type = "response",
                newdata = transform(md1, A = 0))
  Ep0 <- mean(p0)
  
  #Risk ratio
  Ep1 / Ep0
}
```

Next, we call [`boot::boot()`](https://rdrr.io/pkg/boot/man/boot.html)
with this function and the vector of pair membership supplied to perform
the bootstrapping. We’ll request 199 bootstrap replications, but in
practice you should use many more, upwards of 999. More is always
better. Using more also allows you to use the bias-corrected and
accelerated (BCa) boot strap confidence intervals, which are known to be
the most accurate. See
[`?boot.ci`](https://rdrr.io/pkg/boot/man/boot.ci.html) for details.
Here, we’ll just use a percentile confidence interval.

``` r

library("boot")
set.seed(54321)
cluster_boot_out <- boot(pair_ids, cluster_boot_fun,
                         R = 199)

cluster_boot_out
```

    ## 
    ## ORDINARY NONPARAMETRIC BOOTSTRAP
    ## 
    ## 
    ## Call:
    ## boot(data = pair_ids, statistic = cluster_boot_fun, R = 199)
    ## 
    ## 
    ## Bootstrap Statistics :
    ##     original   bias    std. error
    ## t1*    1.588 0.001319      0.1265

``` r

boot.ci(cluster_boot_out, type = "perc")
```

    ## BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS
    ## Based on 199 bootstrap replicates
    ## 
    ## CALL : 
    ## boot.ci(boot.out = cluster_boot_out, type = "perc")
    ## 
    ## Intervals : 
    ## Level     Percentile     
    ## 95%   ( 1.356,  1.857 )  
    ## Calculations and Intervals on Original Scale
    ## Some percentile intervals may be unstable

We find a RR of 1.588 with a confidence interval of (1.356, 1.857). If
we had wanted a risk difference, we could have changed the final line in
`cluster_boot_fun()` to be `Ep1 - Ep0`.

### Moderation Analysis

Moderation analysis involves determining whether a treatment effect
differs across levels of another variable. The use of matching with
moderation analysis is described in Green and Stuart
([2014](#ref-greenExaminingModerationAnalyses2014)). The goal is to
achieve balance within each subgroup of the potential moderating
variable, and there are several ways of doing so. Broadly, one can
either perform matching in the full dataset, requiring exact matching on
the moderator, or one can perform completely separate analyses in each
subgroup. We’ll demonstrate the first approach below; see the blog post
[“Subgroup Analysis After Propensity Score Matching Using
R”](https://ngreifer.github.io/blog/subgroup-analysis-psm/) by Noah
Greifer for an example of the other approach.

There are benefits to using either approach, and Green and Stuart
([2014](#ref-greenExaminingModerationAnalyses2014)) find that either can
be successful at balancing the subgroups. The first approach may be most
effective with small samples, where separate propensity score models
would be fit with greater uncertainty and an increased possibility of
perfect prediction or failure to converge ([Wang et al.
2018](#ref-wangRelativePerformancePropensity2018)). The second approach
may be more effective with larger samples or with matching methods that
target balance in the matched sample, such as genetic matching ([Kreif
et al. 2012](#ref-kreifMethodsEstimatingSubgroup2012)). With genetic
matching, separate subgroup analyses ensure balance is optimized within
each subgroup rather than just overall. The chosen approach should be
that which achieves the best balance, though we don’t demonstrate
assessing balance here to maintain focus on effect estimation.

The full dataset approach involves pooling information across subgroups.
This could involve estimating propensity scores using a single model for
both groups but exact matching on the potential moderator. The
propensity score model could include moderator-by-covariate interactions
to allow the propensity score model to vary across subgroups on some
covariates. It is critical that exact matching is done on the moderator
so that matched pairs are not split across subgroups.

We’ll consider the binary variable `X5` to be the potential moderator of
the effect of `A` on `Y_C`. Below, we’ll estimate a propensity score
using a single propensity score model with a few moderator-by-covariate
interactions. We’ll perform nearest neighbor matching on the propensity
score and exact matching on the moderator, `X5`.

``` r

mP <- matchit(A ~ X1 + X2 + X5*X3 + X4 + 
                X5*X6 + X7 + X5*X8 + X9,
              data = d,
              exact = ~X5)
mP
```

    ## A `matchit` object
    ##  - method: 1:1 nearest neighbor matching without replacement
    ##  - distance: Propensity score
    ##              - estimated with logistic regression
    ##  - number of obs.: 2000 (original), 882 (matched)
    ##  - target estimand: ATT
    ##  - covariates: X1, X2, X5, X3, X4, X6, X7, X8, X9

Although it is straightforward to assess balance overall using
[`summary()`](https://rdrr.io/r/base/summary.html), it is more
challenging to assess balance within subgroups. The easiest way to check
subgroup balance would be to use
[`cobalt::bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.html),
which has a `cluster` argument that can be used to assess balance within
subgroups, e.g., by `cobalt::bal.tab(mP, cluster = "X5")`. See the
vignette “Appendix 2: Using cobalt with Clustered, Multiply Imputed, and
Other Segmented Data” on the `cobalt`
[website](https://ngreifer.github.io/cobalt/index.html) for details.

If we are satisfied with balance, we can then model the outcome with an
interaction between the treatment and the moderator.

``` r

mdP <- match_data(mP)

fitP <- lm(Y_C ~ A * X5, data = mdP, weights = weights)
```

To estimate the subgroup ATTs, we can use `avg_comparisons()`, this time
specifying the `by` argument to signify that we want treatment effects
stratified by the moderator.

``` r

avg_comparisons(fitP,
                variables = "A",
                vcov = ~subclass,
                newdata = subset(A == 1),
                by = "X5")
```

We can see that the subgroup mean differences are quite similar to each
other. Finally, we can test for moderation using another call to
`avg_comparisons()`, this time using the `hypothesis` argument to
signify that we want to compare effects between subgroups:

``` r

avg_comparisons(fitP,
                variables = "A",
                vcov = ~subclass,
                newdata = subset(A == 1),
                by = "X5",
                hypothesis = ~pairwise)
```

As expected, the difference between the subgroup treatment effects is
small and nonsignificant, so there is no evidence of moderation by `X5`.

When the moderator has more than two levels, it is possible to run an
omnibus test for moderation by changing `hypothesis` to `~reference` and
supplying the output to `hypotheses()` with `joint = TRUE`, e.g.,

``` r

avg_comparisons(fitP,
                variables = "A",
                vcov = ~subclass,
                newdata = subset(A == 1),
                by = "X5",
                hypothesis = ~reference) |>
  hypotheses(joint = TRUE)
```

This produces a single p-value for the test that all pairwise
differences between subgroups are equal to zero.

### Reporting Results

It is important to be as thorough and complete as possible when
describing the methods of estimating the treatment effect and the
results of the analysis. This improves transparency and replicability of
the analysis. Results should at least include the following:

- a description of the outcome model used (e.g., logistic regression, a
  linear model with treatment-covariate interactions and covariates, a
  Cox proportional hazards model with the matching weights applied)
- the way the effect was estimated (e.g., using g-computation or as the
  coefficient in the outcome model)
- the way SEs and confidence intervals were estimated (e.g., using
  robust SEs, using cluster-robust SEs with pair membership as the
  cluster, using the BCa bootstrap with 4999 bootstrap replications and
  the entire process of matching and effect estimation included in each
  replication)
- R packages and functions used in estimating the effect and its SE
  (e.g., [`glm()`](https://rdrr.io/r/stats/glm.html) in base R,
  `avg_comparisons()` in `marginaleffects`,
  [`boot()`](https://rdrr.io/pkg/boot/man/boot.html) and
  [`boot.ci()`](https://rdrr.io/pkg/boot/man/boot.ci.html) in `boot`)
- The effect and its SE and confidence interval

All this is in addition to information about the matching method,
propensity score estimation procedure (if used), balance assessment,
etc. mentioned in the other vignettes.

## Common Mistakes

There are a few common mistakes that should be avoided. It is important
not only to avoid these mistakes in one’s own research but also to be
able to spot these mistakes in others’ analyses.

### 1. Failing to include weights

Several methods involve weights that are to be used in estimating the
treatment effect. With full matching and stratification matching (when
analyzed using MMWS), the weights do the entire work of balancing the
covariates across the treatment groups. Omitting weights essentially
ignores the entire purpose of matching. Some cases are less obvious.
When performing matching with replacement and estimating the treatment
effect using the
[`match_data()`](https://kosukeimai.github.io/MatchIt/reference/match_data.md)
output, weights must be included to ensure control units matched to
multiple treated units are weighted accordingly. Similarly, when
performing k:1 matching where not all treated units receive k matches,
weights are required to account for the differential weight of the
matched control units. The only time weights can be omitted after pair
matching is when performing 1:1 matching without replacement. Including
weights even in this scenario will not affect the analysis and it can be
good practice to always include weights to prevent this error from
occurring. There are some scenarios where weights are not useful because
the conditioning occurs through some other means, such as when using the
direct subclass strategy rather than MMWS for estimating marginal
effects after stratification.

### 2. Failing to use robust or cluster-robust standard errors

Robust SEs are required when using weights to estimate the treatment
effect. The model-based SEs resulting from weighted least squares or
maximum likelihood are inaccurate when using matching weights because
they assume weights are frequency weights rather than probability
weights. Cluster-robust SEs account for both the matching weights and
pair membership and should be used when appropriate. Sometimes,
researchers use functions in the `survey` package to estimate robust
SEs, especially with inverse probability weighting; this is a valid way
to compute robust SEs and will give similar results to
[`sandwich::vcovHC()`](https://rdrr.io/pkg/sandwich/man/vcovHC.html).[^10]

### 3. Interpreting conditional effects as marginal effects

The distinction between marginal and conditional effects is not always
clear both in methodological and applied papers. Some statistical
methods are valid only for estimating conditional effects and they
should not be used to estimate marginal effects (without further
modification). Sometimes conditional effects are desirable, and such
methods may be useful for them, but when marginal effects are the target
of inference, it is critical not to inappropriately interpret estimates
resulting from statistical methods aimed at estimating conditional
effects as marginal effects. Although this issue is particularly salient
with binary and survival outcomes due to the general noncollapsibility
of the OR, RR, and HR, this can also occur with linear models for
continuous outcomes or the RD.

The following methods estimate **conditional effects** for binary or
survival outcomes (with noncollapsible effect measures) and should
**not** be used to estimate marginal effects:

- Logistic regression or Cox proportional hazards model with covariates
  and/or the propensity score included, using the coefficient on
  treatment as the effect estimate
- Conditional logistic regression after matching (e.g., using
  [`survival::clogit()`](https://rdrr.io/pkg/survival/man/clogit.html))
- Stratified Cox regression after matching (e.g., using
  [`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html)
  with [`strata()`](https://rdrr.io/pkg/survival/man/strata.html) in the
  model formula)
- Averaging stratum-specific effect estimates after stratification,
  including using Mantel-Haenszel OR pooling
- Including pair or stratum fixed or random effects in a logistic
  regression model, using the coefficient on treatment as the effect
  estimate

In addition, with continuous outcomes, conditional effects can be
mistakenly interpreted as marginal effect estimates when
treatment-covariate interactions are present in the outcome model. If
the covariates are not centered at their mean in the target population
(e.g., the treated group for the ATT, the full sample for the ATE, or
the remaining matched sample for an ATM), the coefficient on treatment
will not correspond to the marginal effect in the target population; it
will correspond to the effect of treatment when the covariate values are
equal to zero, which may not be meaningful or plausible. G-computation
is always the safest way to estimate effects when including covariates
in the outcome model, especially in the presence of treatment-covariate
interactions.

## References

Abadie, Alberto, and Guido W. Imbens. 2008. “On the Failure of the
Bootstrap for Matching Estimators.” *Econometrica* 76 (6): 1537–57.
<https://doi.org/10.3982/ECTA6474>.

Abadie, Alberto, and Jann Spiess. 2019. *Robust Post-Matching
Inference*. January, 34.
<https://doi.org/10.1080/01621459.2020.1840383>.

Austin, Peter C. 2009. “Type i Error Rates, Coverage of Confidence
Intervals, and Variance Estimation in Propensity-Score Matched
Analyses.” *The International Journal of Biostatistics* 5 (1).
<https://doi.org/10.2202/1557-4679.1146>.

Austin, Peter C. 2013a. “The Performance of Different Propensity Score
Methods for Estimating Marginal Hazard Ratios.” *Statistics in Medicine*
32 (16): 2837–49. <https://doi.org/10.1002/sim.5705>.

Austin, Peter C. 2013b. “The Use of Propensity Score Methods with
Survival or Time-to-Event Outcomes: Reporting Measures of Effect Similar
to Those Used in Randomized Experiments.” *Statistics in Medicine* 33
(7): 1242–58. <https://doi.org/10.1002/sim.5984>.

Austin, Peter C. 2017. “Double Propensity-Score Adjustment: A Solution
to Design Bias or Bias Due to Incomplete Matching.” *Statistical Methods
in Medical Research* 26 (1): 201–22.
<https://doi.org/10.1177/0962280214543508>.

Austin, Peter C., and Guy Cafri. 2020. “Variance Estimation When Using
Propensity-Score Matching with Replacement with Survival or
Time-to-Event Outcomes.” *Statistics in Medicine* 39 (11): 1623–40.
<https://doi.org/10.1002/sim.8502>.

Austin, Peter C., and Dylan S. Small. 2014. “The Use of Bootstrapping
When Using Propensity-Score Matching Without Replacement: A Simulation
Study.” *Statistics in Medicine* 33 (24): 4306–19.
<https://doi.org/10.1002/sim.6276>.

Austin, Peter C., and Elizabeth A. Stuart. 2017. “Estimating the Effect
of Treatment on Binary Outcomes Using Full Matching on the Propensity
Score.” *Statistical Methods in Medical Research* 26 (6): 2505–25.
<https://doi.org/10.1177/0962280215601134>.

Austin, Peter C., Neal Thomas, and Donald B. Rubin. 2020.
“Covariate-Adjusted Survival Analyses in Propensity-Score Matched
Samples: Imputing Potential Time-to-Event Outcomes.” *Statistical
Methods in Medical Research* 29 (3): 728–51.
<https://doi.org/10.1177/0962280218817926>.

Bodory, Hugo, Lorenzo Camponovo, Martin Huber, and Michael Lechner.
2020. “The Finite Sample Performance of Inference Methods for Propensity
Score Matching and Weighting Estimators.” *Journal of Business &
Economic Statistics* 38 (1): 183–200.
<https://doi.org/10.1080/07350015.2018.1476247>.

Cameron, A. Colin, and Douglas L. Miller. 2015. “A Practitioner’s Guide
to Cluster-Robust Inference.” *Journal of Human Resources* 50 (2):
317–72. <https://doi.org/10.3368/jhr.50.2.317>.

Carpenter, James, and John Bithell. 2000. “Bootstrap Confidence
Intervals: When, Which, What? A Practical Guide for Medical
Statisticians.” *Statistics in Medicine* 19 (9): 1141–64.
<https://doi.org/10.1002/(SICI)1097-0258(20000515)19:9%3C1141::AID-SIM479%3E3.0.CO;2-F>.

Desai, Rishi J., Kenneth J. Rothman, Brian T. Bateman, Sonia
Hernandez-Diaz, and Krista F. Huybrechts. 2017. “A
Propensity-Score-Based Fine Stratification Approach for Confounding
Adjustment When Exposure Is Infrequent:” *Epidemiology* 28 (2): 249–57.
<https://doi.org/10.1097/EDE.0000000000000595>.

Efron, Bradley, and Robert J. Tibshirani. 1993. *An Introduction to the
Bootstrap*. Springer US.

Gayat, Etienne, Matthieu Resche-Rigon, Jean-Yves Mary, and Raphaël
Porcher. 2012. “Propensity Score Applied to Survival Data Analysis
Through Proportional Hazards Models: A Monte Carlo Study.”
*Pharmaceutical Statistics* 11 (3): 222–29.
<https://doi.org/10.1002/pst.537>.

Green, Kerry M., and Elizabeth A. Stuart. 2014. “Examining Moderation
Analyses in Propensity Score Methods: Application to Depression and
Substance Use.” *Journal of Consulting and Clinical Psychology*,
Advances in Data Analytic Methods, vol. 82 (5): 773–83.
<https://doi.org/10.1037/a0036515>.

Greifer, Noah, and Elizabeth A. Stuart. 2021. “Choosing the Estimand
When Matching or Weighting in Observational Studies.” *arXiv:2106.10577
\[Stat\]*, June. <https://arxiv.org/abs/2106.10577>.

Hill, Jennifer, and Jerome P. Reiter. 2006. “Interval Estimation for
Treatment Effects Using Propensity Score Matching.” *Statistics in
Medicine* 25 (13): 2230–56. <https://doi.org/10.1002/sim.2277>.

Ho, Daniel E., Kosuke Imai, Gary King, and Elizabeth A. Stuart. 2007.
“Matching as Nonparametric Preprocessing for Reducing Model Dependence
in Parametric Causal Inference.” *Political Analysis* 15 (3): 199–236.
<https://doi.org/10.1093/pan/mpl013>.

Hong, Guanglei. 2010. “Marginal Mean Weighting Through Stratification:
Adjustment for Selection Bias in Multilevel Data.” *Journal of
Educational and Behavioral Statistics* 35 (5): 499–531.
<https://doi.org/10.3102/1076998609359785>.

King, Gary, and Margaret E. Roberts. 2015. “How Robust Standard Errors
Expose Methodological Problems They Do Not Fix, and What to Do about
It.” *Political Analysis* 23 (2): 159–79.
<https://doi.org/10.1093/pan/mpu015>.

Kreif, Noemi, Richard Grieve, Rosalba Radice, Zia Sadique, Roland
Ramsahai, and Jasjeet S. Sekhon. 2012. “Methods for Estimating Subgroup
Effects in Cost-Effectiveness Analyses That Use Observational Data.”
*Medical Decision Making* 32 (6): 750–63.
<https://doi.org/10.1177/0272989X12448929>.

Liang, Kung-Yee, and Scott L. Zeger. 1986. “Longitudinal Data Analysis
Using Generalized Linear Models.” *Biometrika* 73 (1): 13–22.
<https://doi.org/10.1093/biomet/73.1.13>.

MacKinnon, James G. 2006. “Bootstrap Methods in Econometrics\*.”
*Economic Record* 82 (s1): S2–18.
<https://doi.org/10.1111/j.1475-4932.2006.00328.x>.

MacKinnon, James G., and Halbert White. 1985. “Some
Heteroskedasticity-Consistent Covariance Matrix Estimators with Improved
Finite Sample Properties.” *Journal of Econometrics* 29 (3): 305–25.
<https://doi.org/10.1016/0304-4076(85)90158-7>.

Nguyen, Tri-Long, Gary S. Collins, Jessica Spence, et al. 2017.
“Double-Adjustment in Propensity Score Matching Analysis: Choosing a
Threshold for Considering Residual Imbalance.” *BMC Medical Research
Methodology* 17: 78. <https://doi.org/10.1186/s12874-017-0338-0>.

Schafer, Joseph L., and Joseph Kang. 2008. “Average Causal Effects from
Nonrandomized Studies: A Practical Guide and Simulated Example.”
*Psychological Methods* 13 (4): 279–313.
<https://doi.org/10.1037/a0014268>.

Snowden, Jonathan M., Sherri Rose, and Kathleen M. Mortimer. 2011.
“Implementation of g-Computation on a Simulated Data Set: Demonstration
of a Causal Inference Technique.” *American Journal of Epidemiology* 173
(7): 731–38. <https://doi.org/10.1093/aje/kwq472>.

Wan, Fei. 2019. “Matched or Unmatched Analyses with
Propensity-Scorematched Data?” *Statistics in Medicine* 38 (2): 289–300.
<https://doi.org/10.1002/sim.7976>.

Wang, Shirley V., Yinzhu Jin, Bruce Fireman, et al. 2018. “Relative
Performance of Propensity Score Matching Strategies for Subgroup
Analyses.” *American Journal of Epidemiology* 187 (8): 1799–807.
<https://doi.org/10.1093/aje/kwy049>.

Westreich, D., and S. Greenland. 2013. “The Table 2 Fallacy: Presenting
and Interpreting Confounder and Modifier Coefficients.” *American
Journal of Epidemiology* 177 (4): 292–98.
<https://doi.org/10.1093/aje/kws412>.

## Code to Generate Data used in Examples

``` r

#Generating data similar to Austin (2009) for demonstrating treatment effect estimation
gen_X <- function(n) {
  X <- matrix(rnorm(9 * n), nrow = n, ncol = 9)
  X[,5] <- as.numeric(X[,5] < .5)
  X
}

#~20% treated
gen_A <- function(X) {
  LP_A <- - 1.2 + log(2)*X[,1] - log(1.5)*X[,2] + log(2)*X[,4] - log(2.4)*X[,5] + log(2)*X[,7] - log(1.5)*X[,8]
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

# Binary outcome
gen_Y_B <- function(A, X) {
  LP_B <- -2 + log(2.4)*A + log(2)*X[,1] + log(2)*X[,2] + log(2)*X[,3] + log(1.5)*X[,4] + log(2.4)*X[,5] + log(1.5)*X[,6]
  P_B <- plogis(LP_B)
  rbinom(length(A), 1, P_B)
}
#Conditional:
#  OR:   2.4
#  logOR: .875
#Marginal:
#  RD:    .144
#  RR:   1.54
#  logRR: .433
#  OR:   1.92
#  logOR  .655

# Survival outcome
gen_Y_S <- function(A, X) {
  LP_S <- -2 + log(2.4)*A + log(2)*X[,1] + log(2)*X[,2] + log(2)*X[,3] + log(1.5)*X[,4] + log(2.4)*X[,5] + log(1.5)*X[,6]
  sqrt(-log(runif(length(A)))*2e4*exp(-LP_S))
}
#Conditional:
#  HR:   2.4
#  logHR: .875
#Marginal:
#  HR:   1.57
#  logHR: .452

set.seed(19599)

n <- 2000
X <- gen_X(n)
A <- gen_A(X)

Y_C <- gen_Y_C(A, X)
Y_B <- gen_Y_B(A, X)
Y_S <- gen_Y_S(A, X)

d <- data.frame(A, X, Y_C, Y_B, Y_S)
```

[^1]: Because they are only appropriate with a large number of clusters,
    cluster-robust SEs are generally not used with subclassification
    methods. Regular robust SEs are valid with these methods when using
    the subclassification weights to estimate marginal effects.

[^2]: Sometimes, an error will occur with this method, which usually
    means more bootstrap replications are required. The number of
    replicates must be greater than the original sample size when using
    the full bootstrap and greater than the number of pairs/strata when
    using the block bootstrap.

[^3]: The matching weights are not necessary when performing 1:1
    matching, but we include them here for generality. When weights are
    not necessary, including them does not affect the estimates. Because
    it may not always be clear when weights are required, we recommend
    always including them.

[^4]: To verify that they are equal, supply the output of
    `avg_predictions()` to `hypotheses()`, e.g.,
    `avg_predictions(...) |> hypotheses(~pairwise)`; this explicitly
    compares the average potential outcomes and should yield identical
    estimates to the `avg_comparisons()` call.

[^5]: It is also known as fine stratification weighting, described by
    Desai et al. ([2017](#ref-desai2017)).

[^6]: We use [`quasibinomial()`](https://rdrr.io/r/stats/family.html)
    instead of [`binomial()`](https://rdrr.io/r/stats/family.html)
    simply to avoid a spurious warning that can occur with certain kinds
    of matching; the results will be identical regardless.

[^7]: Note that for low or high average expected risks computed with
    `avg_predictions()`, the confidence intervals may go below 0 or
    above 1; this is because an approximation is used. To avoid this
    problem, bootstrapping or simulation-based inference can be used
    instead.

[^8]: It is not immediately clear how to estimate a marginal HR when
    covariates are included in the outcome model; though Austin et al.
    ([2020](#ref-austin2020)) describe several ways of including
    covariates in a model to estimate the marginal HR, they do not
    develop SEs and little research has been done on this method, so we
    will not present it here. Instead, we fit a simple Cox model with
    the treatment as the sole predictor.

[^9]: For subclassification, only MMWS can be used; this is done simply
    by including the stratification weights in the Cox model and
    omitting the `cluster` argument.

[^10]: To use `survey` to adjust for pair membership, one can use the
    following code to specify the survey design to be used with
    `svyglm()`:
    `svydesign(ids = ~subclass, weights = ~weights, data = md)` where
    `md` is the output of
    [`match_data()`](https://kosukeimai.github.io/MatchIt/reference/match_data.md).
    After `svyglm()`, `avg_comparisons()` can be used, and the `vcov`
    argument does not need to be specified.
