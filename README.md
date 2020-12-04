
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MatchIt: Nonparametric Preprocessing for Parametric Causal Inference <img src="man/figures/logo.png" align="right" width="150"/>

## [![CRAN\_Status\_Badge](https://img.shields.io/cran/v/MatchIt?color=952100)](https://cran.r-project.org/package=MatchIt) [![CRAN\_Downloads\_Badge](https://cranlogs.r-pkg.org/badges/MatchIt?color=952100)](https://cran.r-project.org/package=MatchIt)

### Overview

`MatchIt` provides a simple and straightforward interface to various
methods of matching for covariate balance in observational studies.
Matching is one way to reduce confounding and model dependence when
estimating treatment effects. Several matching methods are available,
including nearest neighbor matching, optimal pair matching, optimal full
matching, genetic matching, exact matching, coarsened exact matching,
and subclassification, some of which rely on functions from other R
packages. A variety of methods to estimate propensity scores for
propensity score matching are included. Below is an example of the use
of `MatchIt` to perform full matching and assessing overlap and balance:

``` r
library("MatchIt")
data("lalonde", package = "MatchIt")

#Full matching on the propensity score
m.out <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75, 
                 data = lalonde, method = "full")
```

Printing the `MatchIt` object provides details of the kind of matching
performed.

``` r
m.out
```

    #> A matchit object
    #>  - method: Optimal full matching
    #>  - distance: Propensity score
    #>              - estimated with logistic regression
    #>  - number of obs.: 614 (original), 614 (matched)
    #>  - target estimand: ATT
    #>  - covariates: age, educ, race, married, nodegree, re74, re75

We can check covariate balance for the original and matched samples
using `summary()`:

``` r
#Checking balance before and after matching:
summary(m.out)
```

    #> 
    #> Call:
    #> matchit(formula = treat ~ age + educ + race + married + nodegree + 
    #>     re74 + re75, data = lalonde, method = "full")
    #> 
    #> Summary of Balance for All Data:
    #>            Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean eCDF Max
    #> distance          0.5774        0.1822          1.7941     0.9211    0.3774   0.6444
    #> age              25.8162       28.0303         -0.3094     0.4400    0.0813   0.1577
    #> educ             10.3459       10.2354          0.0550     0.4959    0.0347   0.1114
    #> raceblack         0.8432        0.2028          1.7615          .    0.6404   0.6404
    #> racehispan        0.0595        0.1422         -0.3498          .    0.0827   0.0827
    #> racewhite         0.0973        0.6550         -1.8819          .    0.5577   0.5577
    #> married           0.1892        0.5128         -0.8263          .    0.3236   0.3236
    #> nodegree          0.7081        0.5967          0.2450          .    0.1114   0.1114
    #> re74           2095.5737     5619.2365         -0.7211     0.5181    0.2248   0.4470
    #> re75           1532.0553     2466.4844         -0.2903     0.9563    0.1342   0.2876
    #> 
    #> 
    #> Summary of Balance for Matched Data:
    #>            Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean eCDF Max Std. Pair Dist.
    #> distance          0.5774        0.5761          0.0060     0.9918    0.0039   0.0486          0.0192
    #> age              25.8162       24.6928          0.1570     0.4853    0.0838   0.3220          1.2606
    #> educ             10.3459       10.3227          0.0116     0.5577    0.0235   0.0620          1.2200
    #> raceblack         0.8432        0.8347          0.0236          .    0.0086   0.0086          0.0378
    #> racehispan        0.0595        0.0583          0.0049          .    0.0012   0.0012          0.5638
    #> racewhite         0.0973        0.1071         -0.0329          .    0.0098   0.0098          0.4168
    #> married           0.1892        0.1285          0.1549          .    0.0607   0.0607          0.4806
    #> nodegree          0.7081        0.7040          0.0090          .    0.0041   0.0041          0.9143
    #> re74           2095.5737     2199.7126         -0.0213     1.2008    0.0383   0.2350          0.8668
    #> re75           1532.0553     1524.8362          0.0022     2.0048    0.0651   0.2308          0.7932
    #> 
    #> Percent Balance Improvement:
    #>            Std. Mean Diff. Var. Ratio eCDF Mean eCDF Max
    #> distance              99.7       90.0      99.0     92.5
    #> age                   49.3       11.9      -3.0   -104.2
    #> educ                  78.9       16.8      32.4     44.3
    #> raceblack             98.7          .      98.7     98.7
    #> racehispan            98.6          .      98.6     98.6
    #> racewhite             98.2          .      98.2     98.2
    #> married               81.2          .      81.2     81.2
    #> nodegree              96.3          .      96.3     96.3
    #> re74                  97.0       72.2      83.0     47.4
    #> re75                  99.2    -1456.4      51.5     19.8
    #> 
    #> Sample Sizes:
    #>               Control Treated
    #> All            429.       185
    #> Matched (ESS)   53.33     185
    #> Matched        429.       185
    #> Unmatched        0.         0
    #> Discarded        0.         0

At the top is balance for the original sample. Below that is balance in
the matched sample, followed by the percent reduction in imbalance and
the sample sizes before and after matching. Smaller values for the
balance statistics indicate better balance. We can plot the standardized
mean differences in a Love plot for a clean, visual display of balance
across the sample:

``` r
#Plot balance
plot(summary(m.out))
```

<img src="man/figures/README-unnamed-chunk-5-1.png" style="display: block; margin: auto;" />

Although much has been written about matching theory, most of the theory
relied upon in `MatchIt` is described well in [Ho, Imai, King, and
Stuart (2007)](https//:doi.org/10.1093/pan/mpl013) and [Stuart
(2010)](https://doi.org/10.1214/09-STS313). The *Journal of Statistical
Software* article for `MatchIt` can be accessed
[here](https://doi.org/10.18637/jss.v042.i08), though note that some
options have changed, so the `MatchIt` reference pages and included
vignettes should be used for understanding the functions and methods
available. Further references for individual methods are present in
their respective help pages. The `MatchIt`
[website](https://kosukeimai.github.io/MatchIt/) provides access to
vignettes and documentation files.

### Citing `MatchIt`

Please cite `MatchIt` when using it for analysis presented in
publications, which you can do by citing the *Journal of Statistical
Software* article below:

Ho, D. E., Imai, K., King, G., & Stuart, E. A. (2011). MatchIt:
Nonparametric Preprocessing for Parametric Causal Inference. *Journal of
Statistical Software*, 42(8).
[doi:10.18637/jss.v042.i08](https://doi.org/10.18637/jss.v042.i08)

This citation can also be accessed using `citation("MatchIt")` in R. For
reproducibility purposes, it is also important to include the version
number for the version used.
