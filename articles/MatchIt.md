# MatchIt: Getting Started

## Introduction

`MatchIt` implements the suggestions of Ho et al. ([2007](#ref-ho2007))
for improving parametric statistical models for estimating treatment
effects in observational studies and reducing model dependence by
preprocessing data with semi-parametric and non-parametric matching
methods. After appropriately preprocessing with `MatchIt`, researchers
can use whatever parametric model they would have used without `MatchIt`
and produce inferences that are more robust and less sensitive to
modeling assumptions. `MatchIt` reduces the dependence of causal
inferences on commonly made, but hard-to-justify, statistical modeling
assumptions using a large range of sophisticated matching methods. The
package includes several popular approaches to matching and provides
access to methods implemented in other packages through its single,
unified, and easy-to-use interface.

Matching is used in the context of estimating the causal effect of a
binary treatment or exposure on an outcome while controlling for
measured pre-treatment variables, typically confounding variables or
variables prognostic of the outcome. Here and throughout the `MatchIt`
documentation we use the word “treatment” to refer to the focal causal
variable of interest, with “treated” and “control” reflecting the names
of the treatment groups. The goal of matching is to produce *covariate
balance*, that is, for the distributions of covariates in the two groups
to be approximately equal to each other, as they would be in a
successful randomized experiment. The importance of covariate balance is
that it allows for increased robustness to the choice of model used to
estimate the treatment effect; in perfectly balanced samples, a simple
difference in means can be a valid treatment effect estimate. Here we do
not aim to provide a full introduction to matching or causal inference
theory, but simply to explain how to use `MatchIt` to perform
nonparametric preprocessing. For excellent and accessible introductions
to matching, see Stuart ([2010](#ref-stuart2010)) and Austin
([2011](#ref-austin2011b)).

A matching analysis involves four primary steps: 1) planning, 2)
matching, 3) assessing the quality of matches, and 4) estimating the
treatment effect and its uncertainty. Here we briefly discuss these
steps and how they can be implemented with `MatchIt`; in the other
included vignettes, these steps are discussed in more detail.

We will use Lalonde’s data on the evaluation of the National Supported
Work program to demonstrate `MatchIt`’s capabilities. First, we load
`MatchIt` and bring in the `lalonde` dataset.

``` r
library("MatchIt")
data("lalonde")

head(lalonde)
```

    ##      treat age educ   race married nodegree re74 re75       re78
    ## NSW1     1  37   11  black       1        1    0    0  9930.0460
    ## NSW2     1  22    9 hispan       0        1    0    0  3595.8940
    ## NSW3     1  30   12  black       0        0    0    0 24909.4500
    ## NSW4     1  27   11  black       0        1    0    0  7506.1460
    ## NSW5     1  33    8  black       0        1    0    0   289.7899
    ## NSW6     1  22    9  black       0        1    0    0  4056.4940

The statistical quantity of interest is the causal effect of the
treatment (`treat`) on 1978 earnings (`re78`). The other variables are
pre-treatment covariates. See
[`?lalonde`](https://kosukeimai.github.io/MatchIt/reference/lalonde.md)
for more information on this dataset. In particular, the analysis is
concerned with the marginal, total effect of the treatment for those who
actually received the treatment.

In what follows, we briefly describe the four steps of a matching
analysis and how to implement them in `MatchIt`. For more details, we
recommend reading the other vignettes,
[`vignette("matching-methods")`](https://kosukeimai.github.io/MatchIt/articles/matching-methods.md),
[`vignette("assessing-balance")`](https://kosukeimai.github.io/MatchIt/articles/assessing-balance.md),
and
[`vignette("estimating-effects")`](https://kosukeimai.github.io/MatchIt/articles/estimating-effects.md),
especially for users less familiar with matching methods. For the use of
`MatchIt` with sampling weights, also see
[`vignette("sampling-weights")`](https://kosukeimai.github.io/MatchIt/articles/sampling-weights.md).
It is important to recognize that the ease of using `MatchIt` does not
imply the simplicity of matching methods; advanced statistical methods
like matching that require many decisions to be made and caution in
their use should only be performed by those with statistical training.

## Planning

The planning phase of a matching analysis involves selecting the type of
effect to be estimated, selecting the target population to which the
treatment effect is to generalize, and selecting the covariates for
which balance is required for an unbiased estimate of the treatment
effect. Each of these are theoretical steps that do not involve
performing analyses on the data. Ideally, they should be considered
prior to data collection in the planning stage of a study. Thinking
about them early can aid in performing a complete and cost-effective
analysis.

**Selecting the type of effect to be estimated.** There are a few
different types of effects to be estimated. In the presence of mediating
variables, one might be interested in the direct effect of the treatment
that does not pass through the mediating variables or the total effect
of the treatment across all causal pathways. Matching is well suited for
estimating total effects, and specific mediation methods may be better
suited for other mediation-related quantities. One may be interested in
a conditional effect or a marginal effect. A conditional effect is the
effect of a treatment within some strata of other prognostic variables
(e.g., at the patient level), and a marginal effect is the average
effect of a treatment in a population (e.g., for implementing a broad
policy change). Different types of matching are well suited for each of
these, but the most common forms are best used for estimating marginal
treatment effects; for conditional treatment effects, typically modeling
assumptions are required or matching must be done within strata of the
conditioning variables. Matching can reduce the reliance on correct
model specification for conditional effects.

**Selecting a target population.** The target population is the
population to which the effect estimate is to generalize. Typically, an
effect estimated in a sample generalizes to the population from which
the sample is a probability sample. If the sample is not a probability
sample from any population (e.g., it is a convenience sample or involves
patients from an arbitrary hospital), the target population can be
unclear. Often, the target population is a group of units who are
eligible for the treatment (or a subset thereof). Causal estimands are
defined by the target population to which they generalize.

The average treatment effect in the population (ATE) is the average
effect of the treatment for all units in the target population. The
average treatment effect in the treated (ATT) is the average effect of
the treatment for units like those who actually were treated. The most
common forms of matching are best suited for estimating the ATT, though
some are also available for estimating the ATE. Some matching methods
distort the sample in such a way that the estimated treatment effect
corresponds neither to the ATE nor to the ATT, but rather to the effect
in an unspecified population (sometimes called the ATM, or average
treatment effect in the remaining matched sample). When the target
population is not so important (e.g., in the case of treatment effect
discovery), such methods may be attractive; otherwise, care should be
taken in ensuring the effect generalizes to the target population of
interest. Different matching methods allow for different target
populations, so it is important to choose a matching method that allows
one to estimate the desired effect. See Greifer and Stuart
([2021](#ref-greiferChoosingEstimandWhen2021)) for guidance on making
this choice.

**Selecting covariates to balance.** Selecting covariates carefully is
critical for ensuring the resulting treatment effect estimate is free of
confounding and can be validly interpreted as a causal effect. To
estimate total causal effects, all covariates must be measured prior to
treatment (or otherwise not be affected by the treatment). Covariates
should be those that cause variation in the outcome and selection into
treatment group; these are known as confounding variables. See
VanderWeele ([2019](#ref-vanderweele2019)) for a guide on covariate
selection. Ideally these covariates are measured without error and are
free of missingness.

## Check Initial Imbalance

After planning and prior to matching, it can be a good idea to view the
initial imbalance in one’s data that matching is attempting to
eliminate. We can do this using the code below:

``` r
# No matching; constructing a pre-match matchit object
m.out0 <- matchit(treat ~ age + educ + race + married + 
                    nodegree + re74 + re75,
                  data = lalonde,
                  method = NULL,
                  distance = "glm")
```

The first argument is a `formula` relating the treatment to the
covariates used in estimating the propensity score and for which balance
is to be assessed. The `data` argument specifies the dataset where these
variables exist. Typically, the `method` argument specifies the method
of matching to be performed; here, we set it to `NULL` so we can assess
balance prior to matching[¹](#fn1). The `distance` argument specifies
the method for estimating the propensity score, a one-dimensional
summary of all the included covariates, computed as the predicted
probability of being the treated group given the covariates; here, we
set it to `"glm"` for generalized linear model, which implements
logistic regression by default[²](#fn2) (see
[`?distance`](https://kosukeimai.github.io/MatchIt/reference/distance.md)
for other options).

Below we assess balance on the unmatched data using
[`summary()`](https://rdrr.io/r/base/summary.html):

``` r
# Checking balance prior to matching
summary(m.out0)
```

    ## 
    ## Call:
    ## matchit(formula = treat ~ age + educ + race + married + nodegree + 
    ##     re74 + re75, data = lalonde, method = NULL, distance = "glm")
    ## 
    ## Summary of Balance for All Data:
    ##            Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean eCDF Max
    ## distance          0.5774        0.1822          1.7941     0.9211    0.3774   0.6444
    ## age              25.8162       28.0303         -0.3094     0.4400    0.0813   0.1577
    ## educ             10.3459       10.2354          0.0550     0.4959    0.0347   0.1114
    ## raceblack         0.8432        0.2028          1.7615          .    0.6404   0.6404
    ## racehispan        0.0595        0.1422         -0.3498          .    0.0827   0.0827
    ## racewhite         0.0973        0.6550         -1.8819          .    0.5577   0.5577
    ## married           0.1892        0.5128         -0.8263          .    0.3236   0.3236
    ## nodegree          0.7081        0.5967          0.2450          .    0.1114   0.1114
    ## re74           2095.5737     5619.2365         -0.7211     0.5181    0.2248   0.4470
    ## re75           1532.0553     2466.4844         -0.2903     0.9563    0.1342   0.2876
    ## 
    ## Sample Sizes:
    ##           Control Treated
    ## All           429     185
    ## Matched       429     185
    ## Unmatched       0       0
    ## Discarded       0       0

We can see severe imbalances as measured by the standardized mean
differences (`Std. Mean Diff.`), variance ratios (`Var. Ratio`), and
empirical cumulative distribution function (eCDF) statistics. Values of
standardized mean differences and eCDF statistics close to zero and
values of variance ratios close to one indicate good balance, and here
many of them are far from their ideal values.

## Matching

Now, matching can be performed. There are several different classes and
methods of matching, described in
[`vignette("matching-methods")`](https://kosukeimai.github.io/MatchIt/articles/matching-methods.md).
Here, we begin by briefly demonstrating 1:1 nearest neighbor (NN)
matching on the propensity score, which is appropriate for estimating
the ATT. One by one, each treated unit is paired with an available
control unit that has the closest propensity score to it. Any remaining
control units are left unmatched and excluded from further analysis. Due
to the theoretical balancing properties of the propensity score
described by Rosenbaum and Rubin ([1983](#ref-rosenbaum1983)),
propensity score matching can be an effective way to achieve covariate
balance in the treatment groups. Below we demonstrate the use of
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
to perform nearest neighbor propensity score matching.

``` r
# 1:1 NN PS matching w/o replacement
m.out1 <- matchit(treat ~ age + educ + race + married + 
                    nodegree + re74 + re75,
                  data = lalonde,
                  method = "nearest",
                  distance = "glm")
```

We use the same syntax as before, but this time specify
`method = "nearest"` to implement nearest neighbor matching, again using
a logistic regression propensity score. Many other arguments are
available for tuning the matching method and method of propensity score
estimation.

The matching outputs are contained in the `m.out1` object. Printing this
object gives a description of the type of matching performed:

``` r
m.out1
```

    ## A `matchit` object
    ##  - method: 1:1 nearest neighbor matching without replacement
    ##  - distance: Propensity score
    ##              - estimated with logistic regression
    ##  - number of obs.: 614 (original), 370 (matched)
    ##  - target estimand: ATT
    ##  - covariates: age, educ, race, married, nodegree, re74, re75

The key components of the `m.out1` object are `weights` (the computed
matching weights), `subclass` (matching pair membership), `distance`
(the estimated propensity score), and `match.matrix` (which control
units are matched to each treated unit). How these can be used for
estimating the effect of the treatment after matching is detailed in
[`vignette("estimating-effects")`](https://kosukeimai.github.io/MatchIt/articles/estimating-effects.md).

## Assessing the Quality of Matches

Although matching on the propensity score is often effective at
eliminating differences between the treatment groups to achieve
covariate balance, its performance in this regard must be assessed. If
covariates remain imbalanced after matching, the matching is considered
unsuccessful, and a different matching specification should be tried.
`MatchIt` offers a few tools for the assessment of covariate balance
after matching. These include graphical and statistical methods. More
detail on the interpretation of the included plots and statistics can be
found in
[`vignette("assessing-balance")`](https://kosukeimai.github.io/MatchIt/articles/assessing-balance.md).

In addition to covariate balance, the quality of the match is determined
by how many units remain after matching. Matching often involves
discarding units that are not paired with other units, and some matching
options, such as setting restrictions for common support or calipers,
can further decrease the number of remaining units. If, after matching,
the remaining sample size is small, the resulting effect estimate may be
imprecise. In many cases, there will be a trade-off between balance and
remaining sample size. How to optimally choose among them is an instance
of the fundamental bias-variance trade-off problem that cannot be
resolved without substantive knowledge of the phenomena under study.
Prospective power analyses can be used to determine how small a sample
can be before necessary precision is sacrificed.

To assess the quality of the resulting matches numerically, we can use
the [`summary()`](https://rdrr.io/r/base/summary.html) function on
`m.out1` as before. Here we set `un = FALSE` to suppress display of the
balance before matching for brevity and because we already saw it.
(Leaving it as `TRUE`, its default, would display balance both before
and after matching.)

``` r
# Checking balance after NN matching
summary(m.out1, un = FALSE)
```

    ## 
    ## Call:
    ## matchit(formula = treat ~ age + educ + race + married + nodegree + 
    ##     re74 + re75, data = lalonde, method = "nearest", distance = "glm")
    ## 
    ## Summary of Balance for Matched Data:
    ##            Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean eCDF Max Std. Pair Dist.
    ## distance          0.5774        0.3629          0.9739     0.7566    0.1321   0.4216          0.9740
    ## age              25.8162       25.3027          0.0718     0.4568    0.0847   0.2541          1.3938
    ## educ             10.3459       10.6054         -0.1290     0.5721    0.0239   0.0757          1.2474
    ## raceblack         0.8432        0.4703          1.0259          .    0.3730   0.3730          1.0259
    ## racehispan        0.0595        0.2162         -0.6629          .    0.1568   0.1568          1.0743
    ## racewhite         0.0973        0.3135         -0.7296          .    0.2162   0.2162          0.8390
    ## married           0.1892        0.2108         -0.0552          .    0.0216   0.0216          0.8281
    ## nodegree          0.7081        0.6378          0.1546          .    0.0703   0.0703          1.0106
    ## re74           2095.5737     2342.1076         -0.0505     1.3289    0.0469   0.2757          0.7965
    ## re75           1532.0553     1614.7451         -0.0257     1.4956    0.0452   0.2054          0.7381
    ## 
    ## Sample Sizes:
    ##           Control Treated
    ## All           429     185
    ## Matched       185     185
    ## Unmatched     244       0
    ## Discarded       0       0

At the top is a summary of covariate balance after matching. Although
balance has improved for some covariates, in general balance is still
quite poor, indicating that nearest neighbor propensity score matching
is not sufficient for removing confounding in this dataset. The final
column, `Std. Pair Diff`, displays the average absolute within-pair
difference of each covariate. When these values are small, better
balance is typically achieved and estimated effects are more robust to
misspecification of the outcome model ([King and Nielsen
2019](#ref-king2019); [Rubin 1973](#ref-rubin1973a)).

Next is a table of the sample sizes before and after matching. The
matching procedure left 244 control units unmatched. Ideally, unmatched
units would be those far from the treated units and would require
greater extrapolation were they to have been retained. We can visualize
the distribution of propensity scores of those who were matched using
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) with
`type = "jitter"`:

``` r
plot(m.out1, type = "jitter", interactive = FALSE)
```

![Jitter plot of the propensity scores, which shows that no treated unit
were dropped, and a large number of control units with low propensity
scores were dropped.](MatchIt_files/figure-html/unnamed-chunk-7-1.png)

We can visually examine balance on the covariates using
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) with
`type = "density"`:

``` r
plot(m.out1, type = "density", interactive = FALSE,
     which.xs = ~age + married + re75)
```

![Density plots of age, married and re75 in the unmatched and matched
samples.](MatchIt_files/figure-html/unnamed-chunk-8-1.png)

Imbalances are represented by the differences between the black
(treated) and gray (control) distributions. Although `married` and
`re75` appear to have improved balance after matching, the case is mixed
for `age`.

### Trying a Different Matching Specification

Given the poor performance of nearest neighbor matching in this example,
we can try a different matching method or make other changes to the
matching algorithm or distance specification. Below, we’ll try full
matching, which matches every treated unit to at least one control and
every control to at least one treated unit ([Hansen
2004](#ref-hansen2004); [Stuart and Green 2008](#ref-stuart2008a)).
We’ll also try a different link (probit) for the propensity score model.

``` r
# Full matching on a probit PS
m.out2 <- matchit(treat ~ age + educ + race + married + 
                    nodegree + re74 + re75,
                  data = lalonde,
                  method = "full",
                  distance = "glm",
                  link = "probit")
m.out2
```

    ## A `matchit` object
    ##  - method: Optimal full matching
    ##  - distance: Propensity score
    ##              - estimated with probit regression
    ##  - number of obs.: 614 (original), 614 (matched)
    ##  - target estimand: ATT
    ##  - covariates: age, educ, race, married, nodegree, re74, re75

We can examine balance on this new matching specification.

``` r
# Checking balance after full matching
summary(m.out2, un = FALSE)
```

    ## 
    ## Call:
    ## matchit(formula = treat ~ age + educ + race + married + nodegree + 
    ##     re74 + re75, data = lalonde, method = "full", distance = "glm", 
    ##     link = "probit")
    ## 
    ## Summary of Balance for Matched Data:
    ##            Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean eCDF Max Std. Pair Dist.
    ## distance          0.5773        0.5764          0.0045     0.9949    0.0043   0.0486          0.0198
    ## age              25.8162       25.5347          0.0393     0.4790    0.0787   0.2742          1.2843
    ## educ             10.3459       10.5381         -0.0956     0.6192    0.0253   0.0730          1.2179
    ## raceblack         0.8432        0.8389          0.0119          .    0.0043   0.0043          0.0162
    ## racehispan        0.0595        0.0492          0.0435          .    0.0103   0.0103          0.4412
    ## racewhite         0.0973        0.1119         -0.0493          .    0.0146   0.0146          0.3454
    ## married           0.1892        0.1633          0.0660          .    0.0259   0.0259          0.4473
    ## nodegree          0.7081        0.6577          0.1110          .    0.0504   0.0504          0.9872
    ## re74           2095.5737     2100.2150         -0.0009     1.3467    0.0314   0.1881          0.8387
    ## re75           1532.0553     1561.4420         -0.0091     1.5906    0.0536   0.1984          0.8240
    ## 
    ## Sample Sizes:
    ##               Control Treated
    ## All            429.       185
    ## Matched (ESS)   50.76     185
    ## Matched        429.       185
    ## Unmatched        0.         0
    ## Discarded        0.         0

Balance is far better, as determined by the lower standardized mean
differences and eCDF statistics. The balance should be reported when
publishing the results of a matching analysis. This can be done either
in a table, using the values resulting from
[`summary()`](https://rdrr.io/r/base/summary.html), or in a plot, such
as a Love plot, which we can make by calling
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) on the
[`summary()`](https://rdrr.io/r/base/summary.html) output:

``` r
plot(summary(m.out2))
```

![A love plot with matched dots below the threshold lines, indicaitng
good balance after matching, in contrast to the unmatched dots far from
the treshold lines, indicating poor balance before
matching.](MatchIt_files/figure-html/unnamed-chunk-12-1.png)

Love plots are a simple and straightforward way to summarize balance
visually. See
[`vignette("assessing-balance")`](https://kosukeimai.github.io/MatchIt/articles/assessing-balance.md)
for more information on how to customize `MatchIt`’s Love plot and how
to use `cobalt`, a package designed specifically for balance assessment
and reporting that is compatible with `MatchIt`.

## Estimating the Treatment Effect

How treatment effects are estimated depends on what form of matching was
performed. See
[`vignette("estimating-effects")`](https://kosukeimai.github.io/MatchIt/articles/estimating-effects.md)
for information on how to estimate treatment effects in a variety of
scenarios (i.e., different matching methods and outcome types). After
full matching and most other matching methods, we can run a regression
of the outcome on the treatment and covariates in the matched sample
(i.e., including the matching weights) and estimate the treatment effect
using g-computation as implemented in
[`marginaleffects::avg_comparisons()`](https://rdrr.io/pkg/marginaleffects/man/comparisons.html)[³](#fn3).
Including the covariates used in the matching in the effect estimation
can provide additional robustness to slight imbalances remaining after
matching and can improve precision.

Because full matching was successful at balancing the covariates, we’ll
demonstrate here how to estimate a treatment effect after performing
such an analysis. First, we’ll extract the matched dataset from the
`matchit` object using
[`match_data()`](https://kosukeimai.github.io/MatchIt/reference/match_data.md).
This dataset only contains the matched units and adds columns for
`distance`, `weights`, and `subclass` (described previously).

``` r
m.data <- match_data(m.out2)

head(m.data)
```

    ##      treat age educ   race married nodegree re74 re75       re78  distance weights subclass
    ## NSW1     1  37   11  black       1        1    0    0  9930.0460 0.6356769       1        1
    ## NSW2     1  22    9 hispan       0        1    0    0  3595.8940 0.2298151       1       55
    ## NSW3     1  30   12  black       0        0    0    0 24909.4500 0.6813558       1       63
    ## NSW4     1  27   11  black       0        1    0    0  7506.1460 0.7690590       1       70
    ## NSW5     1  33    8  black       0        1    0    0   289.7899 0.6954138       1       79
    ## NSW6     1  22    9  black       0        1    0    0  4056.4940 0.6943658       1       86

We can then model the outcome in this dataset using the standard
regression functions in R, like
[`lm()`](https://rdrr.io/r/stats/lm.html) or
[`glm()`](https://rdrr.io/r/stats/glm.html), being sure to include the
matching weights (stored in the `weights` variable of the
[`match_data()`](https://kosukeimai.github.io/MatchIt/reference/match_data.md)
output) in the estimation[⁴](#fn4). Finally, we use
[`marginaleffects::avg_comparisons()`](https://rdrr.io/pkg/marginaleffects/man/comparisons.html)
to perform g-computation to estimate the ATT. We recommend using
cluster-robust standard errors for most analyses, with pair membership
as the clustering variable;
[`avg_comparisons()`](https://rdrr.io/pkg/marginaleffects/man/comparisons.html)
makes this straightforward.

``` r
library("marginaleffects")

fit <- lm(re78 ~ treat * (age + educ + race + married +
                            nodegree + re74 + re75),
          data = m.data,
          weights = weights)

avg_comparisons(fit,
                variables = "treat",
                vcov = ~subclass,
                newdata = subset(treat == 1))
```

The outcome model coefficients and tests should not be interpreted or
reported. See
[`vignette("estimating-effects")`](https://kosukeimai.github.io/MatchIt/articles/estimating-effects.md)
for more information on how to estimate effects and standard errors with
different forms of matching and with different outcome types.

A benefit of matching is that the outcome model used to estimate the
treatment effect is robust to misspecification when balance has been
achieved. With full matching, we were able to achieve balance, so the
effect estimate should depend less on the form of the outcome model used
than had we used 1:1 matching without replacement or no matching at all.

## Reporting Results

To report matching results in a manuscript or research report, a few key
pieces of information are required. One should be as detailed as
possible about the matching procedure and the decisions made to ensure
the analysis is replicable and can be adequately assessed for soundness
by the audience. Key pieces of information to include are 1) the
matching specification used (including the method and any additional
options, like calipers or common support restrictions), 2) the distance
measure used (including how it was estimated e.g., using logistic
regression for propensity scores), 3) which other matching methods were
tried prior to settling on a final specification and how the choices
were made, 4) the balance of the final matching specification (including
standardized mean differences and other balance statistics for the
variables, their powers, and their interactions; some of these can be
reported as summaries rather than in full detail), 5) the number of
matched, unmatched, and discarded units included in the effect
estimation, and 6) the method of estimating the treatment effect and
standard error or confidence interval (including the specific model used
and the specific type of standard error). See Thoemmes and Kim
([2011](#ref-thoemmes2011)) for a complete list of specific details to
report. Below is an example of how we might write up the prior analysis:

> We used propensity score matching to estimate the average marginal
> effect of the treatment on 1978 earnings on those who received it
> accounting for confounding by the included covariates. We first
> attempted 1:1 nearest neighbor propensity score matching without
> replacement with a propensity score estimated using logistic
> regression of the treatment on the covariates. This matching
> specification yielded poor balance, so we instead tried full matching
> on the propensity score, which yielded adequate balance, as indicated
> in Table 1 and Figure 1. The propensity score was estimated using a
> probit regression of the treatment on the covariates, which yielded
> better balance than did a logistic regression. After matching, all
> standardized mean differences for the covariates were below 0.1 and
> all standardized mean differences for squares and two-way interactions
> between covariates were below .15, indicating adequate balance. Full
> matching uses all treated and all control units, so no units were
> discarded by the matching.
>
> To estimate the treatment effect and its standard error, we fit a
> linear regression model with 1978 earnings as the outcome and the
> treatment, covariates, and their interaction as predictors and
> included the full matching weights in the estimation. The
> [`lm()`](https://rdrr.io/r/stats/lm.html) function was used to fit the
> outcome, and the
> [`avg_comparisons()`](https://rdrr.io/pkg/marginaleffects/man/comparisons.html)
> function in the `marginaleffects` package was used to perform
> g-computation in the matched sample to estimate the ATT. A
> cluster-robust variance was used to estimate its standard error with
> matching stratum membership as the clustering variable.
>
> The estimated effect was \$2114 (SE = 646, p = 0.001), indicating that
> the average effect of the treatment for those who received it is to
> increase earnings.

## Conclusion

Although we have covered the basics of performing a matching analysis
here, to use matching to its full potential, the more advanced methods
available in `MatchIt` should be considered. We recommend reading the
other vignettes included here to gain a better understand of all the
`MatchIt` has to offer and how to use it responsibly and effectively. As
previously stated, the ease of using `MatchIt` does not imply that
matching or causal inference in general are simple matters; matching is
an advanced statistical technique that should be used with care and
caution. We hope the capabilities of `MatchIt` ease and encourage the
use of nonparametric preprocessing for estimating causal effects in a
robust and well-justified way.

## References

Austin, Peter C. 2011. “An Introduction to Propensity Score Methods for
Reducing the Effects of Confounding in Observational Studies.”
*Multivariate Behavioral Research* 46 (3): 399–424.
<https://doi.org/10.1080/00273171.2011.568786>.

Greifer, Noah, and Elizabeth A. Stuart. 2021. “Choosing the Estimand
When Matching or Weighting in Observational Studies.” *arXiv:2106.10577
\[Stat\]*, June. <https://arxiv.org/abs/2106.10577>.

Hansen, Ben B. 2004. “Full Matching in an Observational Study of
Coaching for the SAT.” *Journal of the American Statistical Association*
99 (467): 609–18. <https://doi.org/10.1198/016214504000000647>.

Ho, Daniel E., Kosuke Imai, Gary King, and Elizabeth A. Stuart. 2007.
“Matching as Nonparametric Preprocessing for Reducing Model Dependence
in Parametric Causal Inference.” *Political Analysis* 15 (3): 199–236.
<https://doi.org/10.1093/pan/mpl013>.

King, Gary, and Richard Nielsen. 2019. “Why Propensity Scores Should Not
Be Used for Matching.” *Political Analysis*, May, 1–20.
<https://doi.org/10.1017/pan.2019.11>.

Rosenbaum, Paul R., and Donald B. Rubin. 1983. “The Central Role of the
Propensity Score in Observational Studies for Causal Effects.”
*Biometrika* 70 (1): 41–55. <https://doi.org/10.1093/biomet/70.1.41>.

Rubin, Donald B. 1973. “Matching to Remove Bias in Observational
Studies.” *Biometrics* 29 (1): 159. <https://doi.org/10.2307/2529684>.

Stuart, Elizabeth A. 2010. “Matching Methods for Causal Inference: A
Review and a Look Forward.” *Statistical Science* 25 (1): 1–21.
<https://doi.org/10.1214/09-STS313>.

Stuart, Elizabeth A., and Kerry M. Green. 2008. “Using Full Matching to
Estimate Causal Effects in Nonexperimental Studies: Examining the
Relationship Between Adolescent Marijuana Use and Adult Outcomes.”
*Developmental Psychology* 44 (2): 395–406.
<https://doi.org/10.1037/0012-1649.44.2.395>.

Thoemmes, Felix J., and Eun Sook Kim. 2011. “A Systematic Review of
Propensity Score Methods in the Social Sciences.” *Multivariate
Behavioral Research* 46 (1): 90–118.
<https://doi.org/10.1080/00273171.2011.540475>.

VanderWeele, Tyler J. 2019. “Principles of Confounder Selection.”
*European Journal of Epidemiology* 34 (3): 211–19.
<https://doi.org/10.1007/s10654-019-00494-6>.

------------------------------------------------------------------------

1.  Note that the default for `method` is `"nearest"` to perform nearest
    neighbor matching. To prevent any matching from taking place in
    order to assess pre-matching imbalance, `method` must be set to
    `NULL`.

2.  Note that setting `distance = "logit"`, which was the default in
    `MatchIt` version prior to 4.0.0, or `"ps"`, which was the default
    prior to version 4.5.0, will also estimate logistic regression
    propensity scores. Because it is the default, the `distance`
    argument can actually be omitted if logistic regression propensity
    scores are desired.

3.  In some cases, the coefficient on the treatment variable in the
    outcome model can be used as the effect estimate, but g-computation
    always yields a valid effect estimate regardless of the form of the
    outcome model and its use is the same regardless of the outcome
    model type or matching method (with some slight variations), so we
    always recommend performing g-computation after fitting the outcome
    model. G-computation is explained in detail in
    [`vignette("estimating-effects")`](https://kosukeimai.github.io/MatchIt/articles/estimating-effects.md).

4.  With 1:1 nearest neighbor matching without replacement, excluding
    the matching weights does not change the estimates. For all other
    forms of matching, they are required, so we recommend always
    including them for consistency.
