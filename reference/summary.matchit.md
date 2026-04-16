# View a balance summary of a `matchit` object

Computes and prints balance statistics for `matchit` and
`matchit.subclass` objects. Balance should be assessed to ensure the
matching or subclassification was effective at eliminating treatment
group imbalance and should be reported in the write-up of the results of
the analysis.

## Usage

``` r
# S3 method for class 'matchit'
summary(
  object,
  interactions = FALSE,
  addlvariables = NULL,
  standardize = TRUE,
  data = NULL,
  pair.dist = TRUE,
  un = TRUE,
  improvement = FALSE,
  ...
)

# S3 method for class 'matchit.subclass'
summary(
  object,
  interactions = FALSE,
  addlvariables = NULL,
  standardize = TRUE,
  data = NULL,
  pair.dist = FALSE,
  subclass = FALSE,
  un = TRUE,
  improvement = FALSE,
  ...
)

# S3 method for class 'summary.matchit'
print(x, digits = max(3, getOption("digits") - 3), ...)
```

## Arguments

- object:

  a `matchit` object; the output of a call to
  [`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md).

- interactions:

  `logical`; whether to compute balance statistics for two-way
  interactions and squares of covariates. Default is `FALSE`.

- addlvariables:

  additional variable for which balance statistics are to be computed
  along with the covariates in the `matchit` object. Can be entered in
  one of three ways: as a data frame of covariates with as many rows as
  there were units in the original
  [`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
  call, as a string containing the names of variables in `data`, or as a
  right-sided `formula` with the additional variables (and possibly
  their transformations) found in `data`, the environment, or the
  `matchit` object. Balance on squares and interactions of the
  additional variables will be included if `interactions = TRUE`.

- standardize:

  `logical`; whether to compute standardized (`TRUE`) or unstandardized
  (`FALSE`) statistics. The standardized statistics are the standardized
  mean difference and the mean and maximum of the difference in the
  (weighted) empirical cumulative distribution functions (ECDFs). The
  unstandardized statistics are the raw mean difference and the mean and
  maximum of the quantile-quantile (QQ) difference. Variance ratios are
  produced either way. See Details below. Default is `TRUE`.

- data:

  a optional data frame containing variables named in `addlvariables` if
  specified as a string or formula.

- pair.dist:

  `logical`; whether to compute average absolute pair distances. For
  matching methods that don't include a `match.matrix` component in the
  output (i.e., exact matching, coarsened exact matching, full matching,
  and subclassification), computing pair differences can take a long
  time, especially for large datasets and with many covariates. For
  other methods (i.e., nearest neighbor, optimal, and genetic matching),
  computation is fairly quick. Default is `FALSE` for subclassification
  and `TRUE` otherwise.

- un:

  `logical`; whether to compute balance statistics for the unmatched
  sample. Default `TRUE`; set to `FALSE` for more concise output.

- improvement:

  `logical`; whether to compute the percent reduction in imbalance.
  Default `FALSE`. Ignored if `un = FALSE`.

- ...:

  ignored.

- subclass:

  after subclassification, whether to display balance for individual
  subclasses, and, if so, for which ones. Can be `TRUE` (display balance
  for all subclasses), `FALSE` (display balance only in aggregate), or
  the indices (e.g., `1:6`) of the specific subclasses for which to
  display balance. When anything other than `FALSE`, aggregate balance
  statistics will not be displayed. Default is `FALSE`.

- x:

  a `summay.matchit` or `summary.matchit.subclass` object; the output of
  a call to [`summary()`](https://rdrr.io/r/base/summary.html).

- digits:

  the number of digits to round balance statistics to.

## Value

For `matchit` objects, a `summary.matchit` object, which is a list with
the following components:

- call:

  the original call to
  [`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)

- nn:

  a matrix of the sample sizes in the original (unmatched) and matched
  samples

- sum.all:

  if `un = TRUE`, a matrix of balance statistics for each covariate in
  the original (unmatched) sample

- sum.matched:

  a matrix of balance statistics for each covariate in the matched
  sample

- reduction:

  if `improvement = TRUE`, a matrix of the percent reduction in
  imbalance for each covariate in the matched sample

For `match.subclass` objects, a `summary.matchit.subclass` object, which
is a list as above containing the following components:

- call:

  the original call to
  [`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)

- sum.all:

  if `un = TRUE`, a matrix of balance statistics for each covariate in
  the original sample

- sum.subclass:

  if `subclass` is not `FALSE`, a list of matrices of balance statistics
  for each subclass

- sum.across:

  a matrix of balance statistics for each covariate computed using the
  subclassification weights

- reduction:

  if `improvement = TRUE`, a matrix of the percent reduction in
  imbalance for each covariate in the matched sample

- qn:

  a matrix of sample sizes within each subclass

- nn:

  a matrix of the sample sizes in the original (unmatched) and matched
  samples

## Details

[`summary()`](https://rdrr.io/r/base/summary.html) computes a balance
summary of a `matchit` object. This include balance before and after
matching or subclassification, as well as the percent improvement in
balance. The variables for which balance statistics are computed are
those included in the `formula`, `exact`, and `mahvars` arguments to
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md),
as well as the distance measure if `distance` is was supplied as a
numeric vector or method of estimating propensity scores. The `X`
component of the `matchit` object is used to supply the covariates.

The standardized mean differences are computed both before and after
matching or subclassification as the difference in treatment group means
divided by a standardization factor computed in the unmatched (original)
sample. The standardization factor depends on the argument supplied to
`estimand` in
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md):
for `"ATT"`, it is the standard deviation in the treated group; for
`"ATC"`, it is the standard deviation in the control group; for `"ATE"`,
it is the square root of the average of the variances within each
treatment group. The post-matching mean difference is computed with
weighted means in the treatment groups using the matching or
subclassification weights.

The variance ratio is computed as the ratio of the treatment group
variances. Variance ratios are not computed for binary variables because
their variance is a function solely of their mean. After matching,
weighted variances are computed using the formula used in
[`cov.wt()`](https://rdrr.io/r/stats/cov.wt.html). The percent reduction
in bias is computed using the log of the variance ratios.

The eCDF difference statistics are computed by creating a (weighted)
eCDF for each group and taking the difference between them for each
covariate value. The eCDF is a function that outputs the (weighted)
proportion of units with covariate values at or lower than the input
value. The maximum eCDF difference is the same thing as the
Kolmogorov-Smirnov statistic. The values are bounded at zero and one,
with values closer to zero indicating good overlap between the covariate
distributions in the treated and control groups. For binary variables,
all eCDF differences are equal to the (weighted) difference in
proportion and are computed that way.

The QQ difference statistics are computed by creating two samples of the
same size by interpolating the values of the larger one. The values are
arranged in order for each sample. The QQ difference for each quantile
is the difference between the observed covariate values at that quantile
between the two groups. The difference is on the scale of the original
covariate. Values close to zero indicate good overlap between the
covariate distributions in the treated and control groups. A weighted
interpolation is used for post-matching QQ differences. For binary
variables, all QQ differences are equal to the (weighted) difference in
proportion and are computed that way.

The pair distance is the average of the absolute differences of a
variable between pairs. For example, if a treated unit was paired with
four control units, that set of units would contribute four absolute
differences to the average. Within a subclass, each combination of
treated and control unit forms a pair that contributes once to the
average. The pair distance is described in Stuart and Green (2008) and
is the value that is minimized when using optimal (full) matching. When
`standardize = TRUE`, the standardized versions of the variables are
used, where the standardization factor is as described above for the
standardized mean differences. Pair distances are not computed in the
unmatched sample (because there are no pairs). Because pair distance can
take a while to compute, especially with large datasets or for many
covariates, setting `pair.dist = FALSE` is one way to speed up
[`summary()`](https://rdrr.io/r/base/summary.html).

The effective sample size (ESS) is a measure of the size of a
hypothetical unweighted sample with roughly the same precision as a
weighted sample. When non-uniform matching weights are computed (e.g.,
as a result of full matching, matching with replacement, or
subclassification), the ESS can be used to quantify the potential
precision remaining in the matched sample. The ESS will always be less
than or equal to the matched sample size, reflecting the loss in
precision due to using the weights. With non-uniform weights, it is
printed in the sample size table; otherwise, it is removed because it
does not contain additional information above the matched sample size.

After subclassification, the aggregate balance statistics are computed
using the subclassification weights rather than averaging across
subclasses.

All balance statistics (except pair differences) are computed
incorporating the sampling weights supplied to
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md),
if any. The unadjusted balance statistics include the sampling weights
and the adjusted balance statistics use the matching weights multiplied
by the sampling weights.

When printing, `NA` values are replaced with periods (`.`), and the pair
distance column in the unmatched and percent balance improvement
components of the output are omitted.

## See also

[`summary()`](https://rdrr.io/r/base/summary.html) for the generic
method;
[`plot.summary.matchit()`](https://kosukeimai.github.io/MatchIt/reference/plot.summary.matchit.md)
for making a Love plot from
[`summary()`](https://rdrr.io/r/base/summary.html) output.

[`cobalt::bal.tab.matchit()`](https://ngreifer.github.io/cobalt/reference/bal.tab.matchit.html),
which also displays balance for `matchit` objects.

## Examples

``` r
data("lalonde")
m.out <- matchit(treat ~ age + educ + married +
                   race + re74,
                 data = lalonde,
                 method = "nearest",
                 exact = ~ married,
                 replace = TRUE)

summary(m.out, interactions = TRUE)
#> 
#> Call:
#> matchit(formula = treat ~ age + educ + married + race + re74, 
#>     data = lalonde, method = "nearest", exact = ~married, replace = TRUE)
#> 
#> Summary of Balance for All Data:
#>                      Means Treated Means Control Std. Mean Diff. Var. Ratio
#> distance                    0.5710        0.1850          1.7864     0.8760
#> age                        25.8162       28.0303         -0.3094     0.4400
#> educ                       10.3459       10.2354          0.0550     0.4959
#> married                     0.1892        0.5128         -0.8263          .
#> raceblack                   0.8432        0.2028          1.7615          .
#> racehispan                  0.0595        0.1422         -0.3498          .
#> racewhite                   0.0973        0.6550         -1.8819          .
#> re74                     2095.5737     5619.2365         -0.7211     0.5181
#> age²                      717.3946      901.7786         -0.4276     0.3627
#> age * educ                266.9784      282.3636         -0.1663     0.4912
#> age * married               5.5568       16.4872         -0.9147     0.4615
#> age * raceblack            21.9081        5.2867          1.4327     1.0055
#> age * racehispan            1.3568        3.7646         -0.4368     0.3127
#> age * racewhite             2.5514       18.9790         -2.0322     0.2424
#> age * re74              54074.0365   185650.1507         -0.9974     0.2539
#> educ²                     111.0595      112.8974         -0.0468     0.5173
#> educ * married              1.9622        5.0816         -0.7475     0.5927
#> educ * raceblack            8.6973        2.0466          1.5803     0.9801
#> educ * racehispan           0.5784        1.2634         -0.2940     0.4869
#> educ * racewhite            1.0703        6.9254         -1.7671     0.3652
#> educ * re74             22898.7264    60430.2774         -0.6539     0.5188
#> married * raceblack         0.1568        0.0583          0.2709          .
#> married * racehispan        0.0162        0.0676         -0.4068          .
#> married * racewhite         0.0162        0.3869         -2.9352          .
#> married * re74            760.6329     4324.5356         -0.9734     0.2918
#> raceblack * re74         1817.2003      632.1307          0.2490     3.1701
#> racehispan * re74         151.3968      678.4817         -0.4400     0.1844
#> racewhite * re74          126.9766     4308.6240         -4.5164     0.0198
#> re74²                28141411.5686 77555527.0664         -0.4331     0.6548
#>                      eCDF Mean eCDF Max
#> distance                0.3765   0.6419
#> age                     0.0813   0.1577
#> educ                    0.0347   0.1114
#> married                 0.3236   0.3236
#> raceblack               0.6404   0.6404
#> racehispan              0.0827   0.0827
#> racewhite               0.5577   0.5577
#> re74                    0.2248   0.4470
#> age²                    0.0813   0.1577
#> age * educ              0.0570   0.1187
#> age * married           0.1517   0.3236
#> age * raceblack         0.1831   0.6521
#> age * racehispan        0.0396   0.0827
#> age * racewhite         0.1966   0.5577
#> age * re74              0.2338   0.4470
#> educ²                   0.0347   0.1114
#> educ * married          0.1732   0.3166
#> educ * raceblack        0.3537   0.6451
#> educ * racehispan       0.0457   0.0781
#> educ * racewhite        0.2791   0.5554
#> educ * re74             0.2185   0.4400
#> married * raceblack     0.0985   0.0985
#> married * racehispan    0.0514   0.0514
#> married * racewhite     0.3707   0.3707
#> married * re74          0.1889   0.3626
#> raceblack * re74        0.0861   0.1523
#> racehispan * re74       0.0405   0.0794
#> racewhite * re74        0.2539   0.4959
#> re74²                   0.2248   0.4470
#> 
#> Summary of Balance for Matched Data:
#>                      Means Treated Means Control Std. Mean Diff. Var. Ratio
#> distance                    0.5710        0.5705          0.0023     0.9769
#> age                        25.8162       25.9568         -0.0196     0.4500
#> educ                       10.3459       10.6054         -0.1290     0.5297
#> married                     0.1892        0.1892          0.0000          .
#> raceblack                   0.8432        0.8324          0.0297          .
#> racehispan                  0.0595        0.0595          0.0000          .
#> racewhite                   0.0973        0.1081         -0.0365          .
#> re74                     2095.5737     1882.1670          0.0437     1.5002
#> age²                      717.3946      785.1351         -0.1571     0.3693
#> age * educ                266.9784      268.1514         -0.0127     0.6110
#> age * married               5.5568        5.3459          0.0176     1.0543
#> age * raceblack            21.9081       21.1243          0.0676     0.7602
#> age * racehispan            1.3568        1.7297         -0.0677     0.5429
#> age * racewhite             2.5514        3.1027         -0.0682     0.6518
#> age * re74              54074.0365    56189.4636         -0.0160     0.9095
#> educ²                     111.0595      119.9459         -0.2261     0.5337
#> educ * married              1.9622        2.0432         -0.0194     0.9014
#> educ * raceblack            8.6973        8.7730         -0.0180     0.7955
#> educ * racehispan           0.5784        0.6054         -0.0116     0.8764
#> educ * racewhite            1.0703        1.2270         -0.0473     0.8114
#> educ * re74             22898.7264    21438.6573          0.0254     1.4256
#> married * raceblack         0.1568        0.1514          0.0149          .
#> married * racehispan        0.0162        0.0216         -0.0428          .
#> married * racewhite         0.0162        0.0162          0.0000          .
#> married * re74            760.6329      626.0440          0.0368     2.2366
#> raceblack * re74         1817.2003     1629.7393          0.0394     1.4586
#> racehispan * re74         151.3968       90.2110          0.0511     3.1377
#> racewhite * re74          126.9766      162.2167         -0.0381     1.0624
#> re74²                28141411.5686 19125662.0900          0.0790     2.9783
#>                      eCDF Mean eCDF Max Std. Pair Dist.
#> distance                0.0032   0.0432          0.0165
#> age                     0.0765   0.2216          1.1241
#> educ                    0.0353   0.1892          0.9194
#> married                 0.0000   0.0000          0.0000
#> raceblack               0.0108   0.0108          0.0297
#> racehispan              0.0000   0.0000          0.0216
#> racewhite               0.0108   0.0108          0.0365
#> re74                    0.0544   0.2919          0.3888
#> age²                    0.0765   0.2216          1.1358
#> age * educ              0.0342   0.0973          0.6831
#> age * married           0.0082   0.0649          0.1172
#> age * raceblack         0.0611   0.1946          0.5549
#> age * racehispan        0.0112   0.0270          0.2265
#> age * racewhite         0.0132   0.0270          0.1698
#> age * re74              0.0551   0.2919          0.4300
#> educ²                   0.0353   0.1892          0.9382
#> educ * married          0.0081   0.0811          0.0816
#> educ * raceblack        0.0337   0.1568          0.3956
#> educ * racehispan       0.0047   0.0216          0.1415
#> educ * racewhite        0.0079   0.0216          0.0963
#> educ * re74             0.0507   0.2919          0.4119
#> married * raceblack     0.0054   0.0054          0.0149
#> married * racehispan    0.0054   0.0054          0.0428
#> married * racewhite     0.0000   0.0000          0.0000
#> married * re74          0.0129   0.0865          0.1146
#> raceblack * re74        0.0448   0.2486          0.3194
#> racehispan * re74       0.0055   0.0162          0.1888
#> racewhite * re74        0.0069   0.0378          0.1735
#> re74²                   0.0544   0.2919          0.2207
#> 
#> Sample Sizes:
#>               Control Treated
#> All             429.      185
#> Matched (ESS)    47.6     185
#> Matched          86.      185
#> Unmatched       343.        0
#> Discarded         0.        0
#> 

s.out <- matchit(treat ~ age + educ + married +
                   race + nodegree + re74 + re75,
                 data = lalonde,
                 method = "subclass")

summary(s.out, addlvariables = ~log(age) + I(re74==0))
#> 
#> Call:
#> matchit(formula = treat ~ age + educ + married + race + nodegree + 
#>     re74 + re75, data = lalonde, method = "subclass")
#> 
#> Summary of Balance for All Data:
#>                  Means Treated Means Control Std. Mean Diff. Var. Ratio
#> distance                0.5774        0.1822          1.7941     0.9211
#> age                    25.8162       28.0303         -0.3094     0.4400
#> educ                   10.3459       10.2354          0.0550     0.4959
#> married                 0.1892        0.5128         -0.8263          .
#> raceblack               0.8432        0.2028          1.7615          .
#> racehispan              0.0595        0.1422         -0.3498          .
#> racewhite               0.0973        0.6550         -1.8819          .
#> nodegree                0.7081        0.5967          0.2450          .
#> re74                 2095.5737     5619.2365         -0.7211     0.5181
#> re75                 1532.0553     2466.4844         -0.2903     0.9563
#> log(age)                3.2167        3.2659         -0.1913     0.5093
#> I(re74 == 0)TRUE        0.7081        0.2611          0.9833          .
#>                  eCDF Mean eCDF Max
#> distance            0.3774   0.6444
#> age                 0.0813   0.1577
#> educ                0.0347   0.1114
#> married             0.3236   0.3236
#> raceblack           0.6404   0.6404
#> racehispan          0.0827   0.0827
#> racewhite           0.5577   0.5577
#> nodegree            0.1114   0.1114
#> re74                0.2248   0.4470
#> re75                0.1342   0.2876
#> log(age)            0.0813   0.1577
#> I(re74 == 0)TRUE    0.4470   0.4470
#> 
#> Summary of Balance Across Subclasses
#>                  Means Treated Means Control Std. Mean Diff. Var. Ratio
#> distance                0.5774        0.5610          0.0744     0.8183
#> age                    25.8162       26.2522         -0.0609     0.4375
#> educ                   10.3459       10.2809          0.0323     0.6392
#> married                 0.1892        0.2425         -0.1361          .
#> raceblack               0.8432        0.8279          0.0423          .
#> racehispan              0.0595        0.0360          0.0990          .
#> racewhite               0.0973        0.1361         -0.1309          .
#> nodegree                0.7081        0.6984          0.0214          .
#> re74                 2095.5737     2811.9421         -0.1466     0.8806
#> re75                 1532.0553     1834.5600         -0.0940     1.2233
#> log(age)                3.2167        3.1959          0.0808     0.4883
#> I(re74 == 0)TRUE        0.7081        0.4443          0.5803          .
#>                  eCDF Mean eCDF Max
#> distance            0.0366   0.0815
#> age                 0.0886   0.2558
#> educ                0.0147   0.0499
#> married             0.0533   0.0533
#> raceblack           0.0154   0.0154
#> racehispan          0.0234   0.0234
#> racewhite           0.0388   0.0388
#> nodegree            0.0097   0.0097
#> re74                0.0611   0.2638
#> re75                0.0644   0.2178
#> log(age)            0.0886   0.2558
#> I(re74 == 0)TRUE    0.2638   0.2638
#> 
#> Sample Sizes:
#>               Control Treated
#> All            429.       185
#> Matched (ESS)   66.49     185
#> Matched        429.       185
#> Unmatched        0.         0
#> Discarded        0.         0
#> 

summary(s.out, subclass = TRUE)
#> 
#> Call:
#> matchit(formula = treat ~ age + educ + married + race + nodegree + 
#>     re74 + re75, data = lalonde, method = "subclass")
#> 
#> Summary of Balance for All Data:
#>            Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean
#> distance          0.5774        0.1822          1.7941     0.9211    0.3774
#> age              25.8162       28.0303         -0.3094     0.4400    0.0813
#> educ             10.3459       10.2354          0.0550     0.4959    0.0347
#> married           0.1892        0.5128         -0.8263          .    0.3236
#> raceblack         0.8432        0.2028          1.7615          .    0.6404
#> racehispan        0.0595        0.1422         -0.3498          .    0.0827
#> racewhite         0.0973        0.6550         -1.8819          .    0.5577
#> nodegree          0.7081        0.5967          0.2450          .    0.1114
#> re74           2095.5737     5619.2365         -0.7211     0.5181    0.2248
#> re75           1532.0553     2466.4844         -0.2903     0.9563    0.1342
#>            eCDF Max
#> distance     0.6444
#> age          0.1577
#> educ         0.1114
#> married      0.3236
#> raceblack    0.6404
#> racehispan   0.0827
#> racewhite    0.5577
#> nodegree     0.1114
#> re74         0.4470
#> re75         0.2876
#> 
#> Summary of Balance by Subclass:
#> 
#> - Subclass 1
#>            Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean
#> distance          0.1376        0.0763          0.7782     1.5087    0.2760
#> age              25.5806       28.4595         -0.4173     0.4163    0.0915
#> educ             10.4839       10.2543          0.1056     0.5625    0.0291
#> married           0.2903        0.5780         -0.6339          .    0.2877
#> raceblack         0.0968        0.0145          0.2784          .    0.0823
#> racehispan        0.3226        0.1734          0.3191          .    0.1492
#> racewhite         0.5806        0.8121         -0.4691          .    0.2315
#> nodegree          0.5806        0.5809         -0.0006          .    0.0003
#> re74           3445.7771     6311.1947         -0.3806     1.1631    0.2055
#> re75           2080.7043     2637.3832         -0.1629     1.0676    0.1057
#>            eCDF Max
#> distance     0.4866
#> age          0.1600
#> educ         0.0959
#> married      0.2877
#> raceblack    0.0823
#> racehispan   0.1492
#> racewhite    0.2315
#> nodegree     0.0003
#> re74         0.4293
#> re75         0.2623
#> 
#> - Subclass 2
#>            Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean
#> distance          0.5081        0.4669          0.6017     1.3761    0.1939
#> age              26.9032       32.2500         -0.8559     0.3404    0.1742
#> educ              9.6129        9.3750          0.0872     0.6596    0.0370
#> married           0.4194        0.6667         -0.5012          .    0.2473
#> raceblack         0.9677        0.9583          0.0533          .    0.0094
#> racehispan        0.0323        0.0417         -0.0533          .    0.0094
#> racewhite         0.0000        0.0000          0.0000          .    0.0000
#> nodegree          0.6129        0.6250         -0.0248          .    0.0121
#> re74           6437.9113     6609.9134         -0.0259     0.9023    0.0549
#> re75           3089.0804     3401.3890         -0.0749     0.9876    0.0819
#>            eCDF Max
#> distance     0.3858
#> age          0.4220
#> educ         0.1035
#> married      0.2473
#> raceblack    0.0094
#> racehispan   0.0094
#> racewhite    0.0000
#> nodegree     0.0121
#> re74         0.1882
#> re75         0.2325
#> 
#> - Subclass 3
#>            Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean
#> distance          0.6262        0.6292         -0.1585     1.0943    0.0541
#> age              25.5517       25.1765          0.0452     0.4654    0.1030
#> educ             10.2759       10.0000          0.1312     0.3845    0.0986
#> married           0.4138        0.0588          0.7207          .    0.3550
#> raceblack         1.0000        1.0000          0.0000          .    0.0000
#> racehispan        0.0000        0.0000          0.0000          .    0.0000
#> racewhite         0.0000        0.0000          0.0000          .    0.0000
#> nodegree          0.6897        0.4706          0.4735          .    0.2191
#> re74            806.7050     2110.7067         -0.5707     0.2635    0.1379
#> re75            628.5508      940.9724         -0.2348     0.6585    0.0723
#>            eCDF Max
#> distance     0.1197
#> age          0.2677
#> educ         0.2191
#> married      0.3550
#> raceblack    0.0000
#> racehispan   0.0000
#> racewhite    0.0000
#> nodegree     0.2191
#> re74         0.3225
#> re75         0.1460
#> 
#> - Subclass 4
#>            Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean
#> distance          0.6823        0.6824         -0.0044     1.4231    0.0509
#> age              23.5000       23.8571         -0.0793     0.2048    0.1491
#> educ             10.1875       10.4762         -0.1570     1.0063    0.0599
#> married           0.0312        0.1429         -0.6414          .    0.1116
#> raceblack         1.0000        1.0000          0.0000          .    0.0000
#> racehispan        0.0000        0.0000          0.0000          .    0.0000
#> racewhite         0.0000        0.0000          0.0000          .    0.0000
#> nodegree          0.6250        0.6667         -0.0861          .    0.0417
#> re74           1099.4532     1134.9752         -0.0121     2.2051    0.0857
#> re75           1446.2640     1704.0458         -0.0568     1.7074    0.0762
#>            eCDF Max
#> distance     0.1860
#> age          0.2902
#> educ         0.1711
#> married      0.1116
#> raceblack    0.0000
#> racehispan   0.0000
#> racewhite    0.0000
#> nodegree     0.0417
#> re74         0.3051
#> re75         0.2902
#> 
#> - Subclass 5
#>            Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean
#> distance          0.7296        0.7346         -0.3468     0.7370    0.1180
#> age              24.0000       22.1111          0.2361     0.7413    0.1059
#> educ             10.2903       10.8889         -0.4814     0.5077    0.0737
#> married           0.0000        0.0000          0.0000          .    0.0000
#> raceblack         1.0000        1.0000          0.0000          .    0.0000
#> racehispan        0.0000        0.0000          0.0000          .    0.0000
#> racewhite         0.0000        0.0000          0.0000          .    0.0000
#> nodegree          0.8710        0.8333          0.1123          .    0.0376
#> re74            525.2020      432.2360          0.0781     2.3025    0.0602
#> re75            754.4831      357.4678          0.2722     4.2576    0.0773
#>            eCDF Max
#> distance     0.2832
#> age          0.2885
#> educ         0.1774
#> married      0.0000
#> raceblack    0.0000
#> racehispan   0.0000
#> racewhite    0.0000
#> nodegree     0.0376
#> re74         0.1398
#> re75         0.1613
#> 
#> - Subclass 6
#>            Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean
#> distance          0.7805        0.7774          0.1490     2.8080    0.1250
#> age              29.4194       25.6667          0.5095     0.4039    0.1945
#> educ             11.2258       10.6667          0.3908     6.1419    0.0707
#> married           0.0000        0.0000          0.0000          .    0.0000
#> raceblack         1.0000        1.0000          0.0000          .    0.0000
#> racehispan        0.0000        0.0000          0.0000          .    0.0000
#> racewhite         0.0000        0.0000          0.0000          .    0.0000
#> nodegree          0.8710        1.0000         -0.3849          .    0.1290
#> re74            207.3736      281.4813         -0.1288     1.3924    0.0878
#> re75           1137.7261     1912.6611         -0.3450     2.6542    0.2079
#>            eCDF Max
#> distance     0.3441
#> age          0.6022
#> educ         0.1398
#> married      0.0000
#> raceblack    0.0000
#> racehispan   0.0000
#> racewhite    0.0000
#> nodegree     0.1290
#> re74         0.2043
#> re75         0.7419
#> 
#> Sample Sizes by Subclass:
#>           1  2  3  4  5  6 All
#> Control 346 24 17 21 18  3 429
#> Treated  31 31 29 32 31 31 185
#> Total   377 55 46 53 49 34 614
#> 
```
