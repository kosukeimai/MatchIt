
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MatchIt: Nonparametric Preprocessing for Parametric Causal Inference

[![Build
Status](https://travis-ci.org/kosukeimai/MatchIt.svg?branch=master)](https://travis-ci.org/kosukeimai/MatchIt)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/MatchIt)](https://cran.r-project.org/package=MatchIt)
![CRAN downloads](http://cranlogs.r-pkg.org/badges/grand-total/MatchIt)

`MatchIt` provides a simple and straightforward interface to various
methods of matching for covariate balance in observational studies.
Matching is one way to reduce confounding and model dependence when
estimating treatment effects. Several matching methods are available,
including nearest neighbor matching, optimal pair matching, optimal full
matching, genetic matching, exact matching, coarsened exact matching,
and subclassification, some of which rely on functions from other R
packages. A variety of methods to estimate propensity scores for
propensity score matching are included. Below is an example of the use
of `MatchIt` to perform 2:1 nearest neighbor propensity score matching
with a propensity score caliper and assessing overlap and balance:

``` r
library("MatchIt")
data("lalonde", package = "MatchIt")

#Nearest neighbor PS matching with replacement and a caliper
m.out <- matchit(treat ~ age + educ + race + married + nodegree +
                     re74 + re75, data = lalonde, 
                 ratio = 2, caliper = .025)
```

Printing the `MatchIt` object provides details of the kind of matching
performed.

``` r
print(m.out)
```

    #> A matchit object
    #>  - method: 2:1 nearest neighbor matching without replacement
    #>  - distance: Propensity score [matching, caliper]
    #>               - estimated with logistic regression
    #>  - caliper: 0.025 SD
    #>  - number of obs.: 614 (original), 234 (matched)
    #>  - target estimand: ATT
    #>  - covariates: age, educ, race, married, nodegree, re74, re75

We can view propensity score overlap and see which observations were
matched and unmatched using a jitter plot:

``` r
#Checking for PS overlap
plot(m.out, type = "jitter", interactive = FALSE)
```

<img src="inst/figures/README-unnamed-chunk-4-1.png" style="display: block; margin: auto;" />

With this we can see that most of the unmatched control units had small
propensity scores, making them unlike the treated group. Other plots are
available to view the distributions of propensity scores and covariates.

We can check covariate balance for the original and matched samples
using `summary()`:

``` r
#Checking balance before and after matching:
summary(m.out)
```

    #> 
    #> Call:
    #> matchit(formula = treat ~ age + educ + race + married + nodegree + 
    #>     re74 + re75, data = lalonde, caliper = 0.025, ratio = 2)
    #> 
    #> Summary of Balance for All Data:
    #>            Means Treated Means Control Std. Mean Diff. eCDF Med eCDF Mean eCDF Max
    #> distance          0.5774        0.1822          1.7219   0.3964    0.3774   0.6444
    #> age              25.8162       28.0303         -0.2053   0.0827    0.0813   0.1577
    #> educ             10.3459       10.2354          0.0387   0.0228    0.0347   0.1114
    #> raceblack         0.8432        0.2028          1.5910   0.6404    0.6404   0.6404
    #> racehispan        0.0595        0.1422         -0.2366   0.0827    0.0827   0.0827
    #> racewhite         0.0973        0.6550         -1.1719   0.5577    0.5577   0.5577
    #> married           0.1892        0.5128         -0.6467   0.3236    0.3236   0.3236
    #> nodegree          0.7081        0.5967          0.2268   0.1114    0.1114   0.1114
    #> re74           2095.5737     5619.2365         -0.5190   0.2335    0.2248   0.4470
    #> re75           1532.0553     2466.4844         -0.2838   0.1355    0.1342   0.2876
    #> 
    #> 
    #> Summary of Balance for Matched Data:
    #>            Means Treated Means Control Std. Mean Diff. eCDF Med eCDF Mean eCDF Max
    #> distance          0.4962        0.4937          0.0110   0.0000    0.0043   0.0392
    #> age              25.3627       25.5637         -0.0186   0.0711    0.0800   0.2157
    #> educ             10.2549       10.4216         -0.0584   0.0196    0.0237   0.0784
    #> raceblack         0.7157        0.7157          0.0000   0.0000    0.0000   0.0000
    #> racehispan        0.1078        0.1078          0.0000   0.0000    0.0000   0.0000
    #> racewhite         0.1765        0.1765          0.0000   0.0000    0.0000   0.0000
    #> married           0.2059        0.2157         -0.0196   0.0098    0.0098   0.0098
    #> nodegree          0.6765        0.6422          0.0699   0.0343    0.0343   0.0343
    #> re74           2409.1384     2573.6987         -0.0242   0.0147    0.0528   0.2598
    #> re75           1704.2015     1721.3424         -0.0052   0.0294    0.0391   0.1275
    #> 
    #> Percent Balance Improvement:
    #>            Std. Mean Diff. eCDF Med eCDF Mean eCDF Max
    #> distance              99.4    100.0      98.9     93.9
    #> age                   90.9     14.1       1.6    -36.7
    #> educ                 -50.8     13.9      31.6     29.6
    #> raceblack            100.0    100.0     100.0    100.0
    #> racehispan           100.0    100.0     100.0    100.0
    #> racewhite            100.0    100.0     100.0    100.0
    #> married               97.0     97.0      97.0     97.0
    #> nodegree              69.2     69.2      69.2     69.2
    #> re74                  95.3     93.7      76.5     41.9
    #> re75                  98.2     78.3      70.8     55.7
    #> 
    #> Sample Sizes:
    #>           Control Treated
    #> All           429     185
    #> Matched       132     102
    #> Unmatched     297      83
    #> Discarded       0       0

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

<img src="inst/figures/README-unnamed-chunk-6-1.png" style="display: block; margin: auto;" />

Although much has been written about matching theory, most of the theory
relied upon in `MatchIt` is described well in [Ho, Imai, King, and
Stuart (2007)](https//:doi.org/10.1093/pan/mpl013) and [Stuart
(2010)](https://doi.org/10.1214/09-STS313). The *Journal of Statistical
Software* article for `MatchIt` can be accessed
[here](https://doi.org/10.18637/jss.v042.i08), though note that some
options have changed, so the `MatchIt` reference pages should be used
for understanding the functions and methods available. Further
references for individual methods are present in their respective help
pages.

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
