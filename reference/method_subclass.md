# Subclassification

In
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md),
setting `method = "subclass"` performs subclassification on the distance
measure (i.e., propensity score). Treatment and control units are placed
into subclasses based on quantiles of the propensity score in the
treated group, in the control group, or overall, depending on the
desired estimand. Weights are computed based on the proportion of
treated units in each subclass. Subclassification implemented here does
not rely on any other package.

This page details the allowable arguments with `method = "subclass"`.
See
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
for an explanation of what each argument means in a general context and
how it can be specified.

Below is how
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
is used for subclassification:


    matchit(formula,
            data = NULL,
            method = "subclass",
            distance = "glm",
            link = "logit",
            distance.options = list(),
            estimand = "ATT",
            discard = "none",
            reestimate = FALSE,
            s.weights = NULL,
            verbose = FALSE,
            ...) 

## Arguments

- formula:

  a two-sided [formula](https://rdrr.io/r/stats/formula.html) object
  containing the treatment and covariates to be used in creating the
  distance measure used in the subclassification.

- data:

  a data frame containing the variables named in `formula`. If not found
  in `data`, the variables will be sought in the environment.

- method:

  set here to `"subclass"`.

- distance:

  the distance measure to be used. See
  [`distance`](https://kosukeimai.github.io/MatchIt/reference/distance.md)
  for allowable options. Must be a vector of distance scores or the name
  of a method of estimating propensity scores.

- link:

  when `distance` is specified as a string, an additional argument
  controlling the link function used in estimating the distance measure.
  See
  [`distance`](https://kosukeimai.github.io/MatchIt/reference/distance.md)
  for allowable options with each option.

- distance.options:

  a named list containing additional arguments supplied to the function
  that estimates the distance measure as determined by the argument to
  `distance`.

- estimand:

  the target `estimand`. If `"ATT"`, the default, subclasses are formed
  based on quantiles of the distance measure in the treated group; if
  `"ATC"`, subclasses are formed based on quantiles of the distance
  measure in the control group; if `"ATE"`, subclasses are formed based
  on quantiles of the distance measure in the full sample. The estimand
  also controls how the subclassification weights are computed; see the
  Computing Weights section at
  [`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
  for details.

- discard:

  a string containing a method for discarding units outside a region of
  common support.

- reestimate:

  if `discard` is not `"none"`, whether to re-estimate the propensity
  score in the remaining sample prior to subclassification.

- s.weights:

  the variable containing sampling weights to be incorporated into
  propensity score models and balance statistics.

- verbose:

  `logical`; whether information about the matching process should be
  printed to the console.

- ...:

  additional arguments that control the subclassification:

  `subclass`

  :   either the number of subclasses desired or a vector of quantiles
      used to divide the distance measure into subclasses. Default is 6.

  `min.n`

  :   the minimum number of units of each treatment group that are to be
      assigned each subclass. If the distance measure is divided in such
      a way that fewer than `min.n` units of a treatment group are
      assigned a given subclass, units from other subclasses will be
      reassigned to fill the deficient subclass. Default is 1.

  The arguments `exact`, `mahvars`, `replace`, `m.order`, `caliper` (and
  related arguments), and `ratio` are ignored with a warning.

## Details

After subclassification, effect estimates can be computed separately in
the subclasses and combined, or a single marginal effect can be
estimated by using the weights in the full sample. When using the
weights, the method is sometimes referred to as marginal mean weighting
through stratification (MMWS; Hong, 2010) or fine stratification
weighting (Desai et al., 2017). The weights can be interpreted just like
inverse probability weights. See
[`vignette("estimating-effects")`](https://kosukeimai.github.io/MatchIt/articles/estimating-effects.md)
for details.

Changing `min.n` can change the quality of the weights. Generally, a low
`min.w` will yield better balance because subclasses only contain units
with relatively similar distance values, but may yield higher variance
because extreme weights can occur due to there being few members of a
treatment group in some subclasses. When `min.n = 0`, some subclasses
may fail to contain units from both treatment groups, in which case all
units in such subclasses will be dropped.

Note that subclassification weights can also be estimated using
*WeightIt*, which provides some additional methods for estimating
propensity scores. Where propensity score-estimation methods overlap,
both packages will yield the same weights.

## Outputs

All outputs described in
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
are returned with `method = "subclass"` except that `match.matrix` is
excluded and one additional component, `q.cut`, is included, containing
a vector of the distance measure cutpoints used to define the
subclasses. Note that when `min.n > 0`, the subclass assignments may not
strictly obey the quantiles listed in `q.cut`. `include.obj` is ignored.

## References

In a manuscript, you don't need to cite another package when using
`method = "subclass"` because the subclassification is performed
completely within *MatchIt*. For example, a sentence might read:

*Propensity score subclassification was performed using the MatchIt
package (Ho, Imai, King, & Stuart, 2011) in R.*

It may be a good idea to cite Hong (2010) or Desai et al. (2017) if the
treatment effect is estimated using the subclassification weights.

Desai, R. J., Rothman, K. J., Bateman, B. . T., Hernandez-Diaz, S., &
Huybrechts, K. F. (2017). A Propensity-score-based Fine Stratification
Approach for Confounding Adjustment When Exposure Is Infrequent:
Epidemiology, 28(2), 249–257.
[doi:10.1097/EDE.0000000000000595](https://doi.org/10.1097/EDE.0000000000000595)

Hong, G. (2010). Marginal mean weighting through stratification:
Adjustment for selection bias in multilevel data. Journal of Educational
and Behavioral Statistics, 35(5), 499–531.
[doi:10.3102/1076998609359785](https://doi.org/10.3102/1076998609359785)

## See also

[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
for a detailed explanation of the inputs and outputs of a call to
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md).

[`method_full`](https://kosukeimai.github.io/MatchIt/reference/method_full.md)
for optimal full matching and
[`method_quick`](https://kosukeimai.github.io/MatchIt/reference/method_quick.md)
for generalized full matching, which are similar to subclassification
except that the number of subclasses and subclass membership are chosen
to optimize the within-subclass distance.

## Examples

``` r

data("lalonde")

# PS subclassification for the ATT with 7 subclasses
s.out1 <- matchit(treat ~ age + educ + race + nodegree +
                    married + re74 + re75,
                  data = lalonde,
                  method = "subclass",
                  subclass = 7)
s.out1
#> A `matchit` object
#>  - method: Subclassification (7 subclasses)
#>  - distance: Propensity score
#>              - estimated with logistic regression
#>  - number of obs.: 614 (original), 614 (matched)
#>  - target estimand: ATT
#>  - covariates: age, educ, race, nodegree, married, re74, re75
summary(s.out1, subclass = TRUE)
#> 
#> Call:
#> matchit(formula = treat ~ age + educ + race + nodegree + married + 
#>     re74 + re75, data = lalonde, method = "subclass", subclass = 7)
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
#> Summary of Balance by Subclass:
#> 
#> - Subclass 1
#>            Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean
#> distance          0.1141        0.0674          0.9626     1.0379    0.2771
#> age              24.4815       28.6536         -0.6501     0.3615    0.1150
#> educ             10.7407       10.2590          0.2544     0.4213    0.0401
#> raceblack         0.0370        0.0030          0.1802          .    0.0340
#> racehispan        0.2963        0.1506          0.3191          .    0.1457
#> racewhite         0.6667        0.8464         -0.3812          .    0.1797
#> nodegree          0.5556        0.5723         -0.0337          .    0.0167
#> married           0.2593        0.5904         -0.7555          .    0.3311
#> re74           3205.1533     6429.1274         -0.4426     1.0648    0.2177
#> re75           1983.7216     2631.3105         -0.2054     0.9066    0.0980
#>            eCDF Max
#> distance     0.4934
#> age          0.2446
#> educ         0.1169
#> raceblack    0.0340
#> racehispan   0.1457
#> racewhite    0.1797
#> nodegree     0.0167
#> married      0.3311
#> re74         0.4266
#> re75         0.2263
#> 
#> - Subclass 2
#>            Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean
#> distance          0.4522        0.3958          0.6193     0.8134    0.1707
#> age              29.1154       29.5135         -0.0603     0.3486    0.1308
#> educ              9.1538        9.7568         -0.2064     0.9233    0.0599
#> raceblack         0.8846        0.7027          0.5694          .    0.1819
#> racehispan        0.1154        0.2973         -0.5694          .    0.1819
#> racewhite         0.0000        0.0000          0.0000          .    0.0000
#> nodegree          0.6923        0.6757          0.0360          .    0.0166
#> married           0.5000        0.5405         -0.0811          .    0.0405
#> re74           6754.9455     5617.3241          0.1467     1.5744    0.0817
#> re75           3312.9595     3257.5161          0.0116     1.5170    0.0581
#>            eCDF Max
#> distance     0.2994
#> age          0.3129
#> educ         0.1674
#> raceblack    0.1819
#> racehispan   0.1819
#> racewhite    0.0000
#> nodegree     0.0166
#> married      0.0405
#> re74         0.2453
#> re75         0.1746
#> 
#> - Subclass 3
#>            Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean
#> distance          0.5998        0.6097         -0.4590     1.2565    0.1421
#> age              26.6923       25.9000          0.1037     0.2696    0.2119
#> educ             10.3846        8.8000          0.7416     0.3682    0.1731
#> raceblack         1.0000        1.0000          0.0000          .    0.0000
#> racehispan        0.0000        0.0000          0.0000          .    0.0000
#> racewhite         0.0000        0.0000          0.0000          .    0.0000
#> nodegree          0.6923        0.5000          0.4167          .    0.1923
#> married           0.5000        0.0000          1.0000          .    0.5000
#> re74           2416.9731     1748.8415          0.1851     2.6567    0.1324
#> re75           1350.8111     1004.3709          0.1863     1.8285    0.0677
#>            eCDF Max
#> distance     0.3615
#> age          0.5077
#> educ         0.3846
#> raceblack    0.0000
#> racehispan   0.0000
#> racewhite    0.0000
#> nodegree     0.1923
#> married      0.5000
#> re74         0.2769
#> re75         0.2077
#> 
#> - Subclass 4
#>            Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean
#> distance          0.6582        0.6588         -0.0428     0.9551    0.0452
#> age              22.2963       24.2941         -0.4195     0.2181    0.1053
#> educ             10.4074       11.0588         -0.3613     0.5596    0.0787
#> raceblack         1.0000        1.0000          0.0000          .    0.0000
#> racehispan        0.0000        0.0000          0.0000          .    0.0000
#> racewhite         0.0000        0.0000          0.0000          .    0.0000
#> nodegree          0.4815        0.5294         -0.0959          .    0.0479
#> married           0.0370        0.2353         -1.0498          .    0.1983
#> re74            663.6473     1733.6054         -0.6576     0.1297    0.0617
#> re75            498.2558     1882.9978         -1.0580     0.1094    0.1153
#>            eCDF Max
#> distance     0.1307
#> age          0.1808
#> educ         0.1765
#> raceblack    0.0000
#> racehispan   0.0000
#> racewhite    0.0000
#> nodegree     0.0479
#> married      0.1983
#> re74         0.1678
#> re75         0.2048
#> 
#> - Subclass 5
#>            Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean
#> distance          0.7013        0.6979          0.2860     0.8887    0.0765
#> age              24.3846       22.2941          0.2995     0.7142    0.1296
#> educ             10.1154        9.8824          0.1297     1.3668    0.0637
#> raceblack         1.0000        1.0000          0.0000          .    0.0000
#> racehispan        0.0000        0.0000          0.0000          .    0.0000
#> racewhite         0.0000        0.0000          0.0000          .    0.0000
#> nodegree          0.7308        0.8235         -0.2091          .    0.0928
#> married           0.0385        0.0000          0.2000          .    0.0385
#> re74            932.1795      812.5195          0.0403     3.9223    0.0899
#> re75           1551.4350      586.9098          0.1945    23.9287    0.0601
#>            eCDF Max
#> distance     0.2805
#> age          0.3032
#> educ         0.0928
#> raceblack    0.0000
#> racehispan   0.0000
#> racewhite    0.0000
#> nodegree     0.0928
#> married      0.0385
#> re74         0.2783
#> re75         0.2014
#> 
#> - Subclass 6
#>            Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean
#> distance          0.7396        0.7425         -0.2126     1.4075    0.1214
#> age              23.8846       22.4615          0.1747     0.6952    0.1058
#> educ             10.3462       11.3846         -0.8436     0.4905    0.1202
#> raceblack         1.0000        1.0000          0.0000          .    0.0000
#> racehispan        0.0000        0.0000          0.0000          .    0.0000
#> racewhite         0.0000        0.0000          0.0000          .    0.0000
#> nodegree          0.9231        0.7692          0.5774          .    0.1538
#> married           0.0000        0.0000          0.0000          .    0.0000
#> re74            591.5328      517.2464          0.0596     2.0185    0.0545
#> re75            974.8401      479.1179          0.2851     4.6722    0.0769
#>            eCDF Max
#> distance     0.2692
#> age          0.3462
#> educ         0.3077
#> raceblack    0.0000
#> racehispan   0.0000
#> racewhite    0.0000
#> nodegree     0.1538
#> married      0.0000
#> re74         0.1538
#> re75         0.1923
#> 
#> - Subclass 7
#>            Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean
#> distance          0.7837        0.7774          0.2999     2.7433    0.1407
#> age              29.8889       25.6667          0.6648     0.3002    0.2333
#> educ             11.2222       10.6667          0.3763     6.5385    0.0679
#> raceblack         1.0000        1.0000          0.0000          .    0.0000
#> racehispan        0.0000        0.0000          0.0000          .    0.0000
#> racewhite         0.0000        0.0000          0.0000          .    0.0000
#> nodegree          0.8889        1.0000         -0.3536          .    0.1111
#> married           0.0000        0.0000          0.0000          .    0.0000
#> re74            190.2635      281.4813         -0.1585     1.3943    0.0889
#> re75           1091.6910     1912.6611         -0.3734     2.5440    0.2189
#>            eCDF Max
#> distance     0.2963
#> age          0.6667
#> educ         0.1481
#> raceblack    0.0000
#> racehispan   0.0000
#> racewhite    0.0000
#> nodegree     0.1111
#> married      0.0000
#> re74         0.2222
#> re75         0.7407
#> 
#> Sample Sizes by Subclass:
#>           1  2  3  4  5  6  7 All
#> Control 332 37 10 17 17 13  3 429
#> Treated  27 26 26 27 26 26 27 185
#> Total   359 63 36 44 43 39 30 614
#> 

# PS subclassification for the ATE with 10 subclasses
# and at least 2 units in each group per subclass
s.out2 <- matchit(treat ~ age + educ + race + nodegree +
                    married + re74 + re75,
                  data = lalonde,
                  method = "subclass",
                  subclass = 10,
                  estimand = "ATE",
                  min.n = 2)
s.out2
#> A `matchit` object
#>  - method: Subclassification (10 subclasses)
#>  - distance: Propensity score
#>              - estimated with logistic regression
#>  - number of obs.: 614 (original), 614 (matched)
#>  - target estimand: ATE
#>  - covariates: age, educ, race, nodegree, married, re74, re75
summary(s.out2)
#> 
#> Call:
#> matchit(formula = treat ~ age + educ + race + nodegree + married + 
#>     re74 + re75, data = lalonde, method = "subclass", estimand = "ATE", 
#>     subclass = 10, min.n = 2)
#> 
#> Summary of Balance for All Data:
#>            Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean
#> distance          0.5774        0.1822          1.7569     0.9211    0.3774
#> age              25.8162       28.0303         -0.2419     0.4400    0.0813
#> educ             10.3459       10.2354          0.0448     0.4959    0.0347
#> raceblack         0.8432        0.2028          1.6708          .    0.6404
#> racehispan        0.0595        0.1422         -0.2774          .    0.0827
#> racewhite         0.0973        0.6550         -1.4080          .    0.5577
#> nodegree          0.7081        0.5967          0.2355          .    0.1114
#> married           0.1892        0.5128         -0.7208          .    0.3236
#> re74           2095.5737     5619.2365         -0.5958     0.5181    0.2248
#> re75           1532.0553     2466.4844         -0.2870     0.9563    0.1342
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
#> Summary of Balance Across Subclasses
#>            Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean
#> distance          0.3104        0.2972          0.0585     1.0249    0.0454
#> age              25.1099       27.1896         -0.2272     0.3247    0.0996
#> educ             10.5836       10.3105          0.1106     0.5537    0.0332
#> raceblack         0.4040        0.3898          0.0371          .    0.0142
#> racehispan        0.0894        0.1233         -0.1139          .    0.0340
#> racewhite         0.5066        0.4869          0.0498          .    0.0197
#> nodegree          0.5194        0.6163         -0.2049          .    0.0969
#> married           0.2688        0.4150         -0.3257          .    0.1462
#> re74           2489.2647     4627.1084         -0.3614     0.5922    0.1342
#> re75           1532.5437     2199.1568         -0.2047     0.9287    0.0962
#>            eCDF Max
#> distance     0.1674
#> age          0.2203
#> educ         0.0969
#> raceblack    0.0142
#> racehispan   0.0340
#> racewhite    0.0197
#> nodegree     0.0969
#> married      0.1462
#> re74         0.3399
#> re75         0.1839
#> 
#> Sample Sizes:
#>               Control Treated
#> All            429.    185.  
#> Matched (ESS)  332.11   42.66
#> Matched        429.    185.  
#> Unmatched        0.      0.  
#> Discarded        0.      0.  
#> 
```
