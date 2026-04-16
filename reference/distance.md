# Propensity scores and other distance measures

Several matching methods require or can involve the distance between
treated and control units. Options include the Mahalanobis distance,
propensity score distance, or distance between user-supplied values.
Propensity scores are also used for common support via the `discard`
options and for defining calipers. This page documents the options that
can be supplied to the `distance` argument to
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md).

## Note

In versions of *MatchIt* prior to 4.0.0, `distance` was specified in a
slightly different way. When specifying arguments using the old syntax,
they will automatically be converted to the corresponding method in the
new syntax but a warning will be thrown. `distance = "logit"`, the old
default, will still work in the new syntax, though
`distance = "glm", link = "logit"` is preferred (note that these are the
default settings and don't need to be made explicit).

## Allowable options

There are four ways to specify the `distance` argument: 1) as a string
containing the name of a method for estimating propensity scores, 2) as
a string containing the name of a method for computing pairwise
distances from the covariates, 3) as a vector of values whose pairwise
differences define the distance between units, or 4) as a distance
matrix containing all pairwise distances. The options are detailed
below.

### Propensity score estimation methods

When `distance` is specified as the name of a method for estimating
propensity scores (described below), a propensity score is estimated
using the variables in `formula` and the method corresponding to the
given argument. This propensity score can be used to compute the
distance between units as the absolute difference between the propensity
scores of pairs of units. Propensity scores can also be used to create
calipers and common support restrictions, whether or not they are used
in the actual distance measure used in the matching, if any.

In addition to the `distance` argument, two other arguments can be
specified that relate to the estimation and manipulation of the
propensity scores. The `link` argument allows for different links to be
used in models that require them such as generalized linear models, for
which the logit and probit links are allowed, among others. In addition
to specifying the link, the `link` argument can be used to specify
whether the propensity score or the linearized version of the propensity
score should be used (i.e., the linear predictor of the propensity score
model); by specifying `link = "linear.{link}"`, the linearized version
will be used. When `link = "linear.logit"`, for example, this requests
the logit of a propensity score estimated with a logistic link.

The `distance.options` argument can also be specified, which should be a
list of values passed to the propensity score-estimating function, for
example, to choose specific options or tuning parameters for the
estimation method. If `formula`, `data`, or `verbose` are not supplied
to `distance.options`, the corresponding arguments from
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
will be automatically supplied. See the Examples for demonstrations of
the uses of `link` and `distance.options`. When `s.weights` is supplied
in the call to
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md),
it will automatically be passed to the propensity score-estimating
function as the `weights` argument unless otherwise described below.

The following methods for estimating propensity scores are allowed:

- `"glm"`:

  The propensity scores are estimated using a generalized linear model
  (e.g., logistic regression). The `formula` supplied to
  [`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
  is passed directly to [`glm()`](https://rdrr.io/r/stats/glm.html), and
  [`predict.glm()`](https://rdrr.io/r/stats/predict.glm.html) is used to
  compute the propensity scores. The `link` argument can be specified as
  a link function supplied to
  [`binomial()`](https://rdrr.io/r/stats/family.html), e.g., `"logit"`,
  which is the default. When `link` is prepended by `"linear."`, the
  linear predictor is used instead of the predicted probabilities.
  `distance = "glm"` with `link = "logit"` (logistic regression) is the
  default in
  [`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md).
  (This used to be able to be requested as `distance = "ps"`, which
  still works.)

- `"gam"`:

  The propensity scores are estimated using a generalized additive
  model. The `formula` supplied to
  [`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
  is passed directly to
  [`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html), and
  [`mgcv::predict.gam()`](https://rdrr.io/pkg/mgcv/man/predict.gam.html)
  is used to compute the propensity scores. The `link` argument can be
  specified as a link function supplied to
  [`binomial()`](https://rdrr.io/r/stats/family.html), e.g., `"logit"`,
  which is the default. When `link` is prepended by `"linear."`, the
  linear predictor is used instead of the predicted probabilities. Note
  that unless the smoothing functions
  [`mgcv::s()`](https://rdrr.io/pkg/mgcv/man/s.html),
  [`mgcv::te()`](https://rdrr.io/pkg/mgcv/man/te.html),
  [`mgcv::ti()`](https://rdrr.io/pkg/mgcv/man/te.html), or
  [`mgcv::t2()`](https://rdrr.io/pkg/mgcv/man/t2.html) are used in
  `formula`, a generalized additive model is identical to a generalized
  linear model and will estimate the same propensity scores as
  [`glm()`](https://rdrr.io/r/stats/glm.html). See the documentation for
  [`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html),
  [`mgcv::formula.gam()`](https://rdrr.io/pkg/mgcv/man/formula.gam.html),
  and
  [`mgcv::gam.models()`](https://rdrr.io/pkg/mgcv/man/gam.models.html)
  for more information on how to specify these models. Also note that
  the formula returned in the
  [`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
  output object will be a simplified version of the supplied formula
  with smoothing terms removed (but all named variables present).

- `"gbm"`:

  The propensity scores are estimated using a generalized boosted model.
  The `formula` supplied to
  [`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
  is passed directly to
  [`gbm::gbm()`](https://rdrr.io/pkg/gbm/man/gbm.html), and
  [`gbm::predict.gbm()`](https://rdrr.io/pkg/gbm/man/predict.gbm.html)
  is used to compute the propensity scores. The optimal tree is chosen
  using 5-fold cross-validation by default, and this can be changed by
  supplying an argument to `method` to `distance.options`; see
  [`gbm::gbm.perf()`](https://rdrr.io/pkg/gbm/man/gbm.perf.html) for
  details. The `link` argument can be specified as `"linear"` to use the
  linear predictor instead of the predicted probabilities. No other
  links are allowed. The tuning parameter defaults differ from
  [`gbm::gbm()`](https://rdrr.io/pkg/gbm/man/gbm.html); they are as
  follows: `n.trees = 1e4`, `interaction.depth = 3`, `shrinkage = .01`,
  `bag.fraction = 1`, `cv.folds = 5`, `keep.data = FALSE`. These are the
  same defaults as used in *WeightIt* and *twang*, except for `cv.folds`
  and `keep.data`. Note this is not the same use of generalized boosted
  modeling as in *twang*; here, the number of trees is chosen based on
  cross-validation or out-of-bag error, rather than based on optimizing
  balance. twang should not be cited when using this method to estimate
  propensity scores. Note that because there is a random component to
  choosing the tuning parameter, results will vary across runs unless a
  [seed](https://rdrr.io/r/base/Random.html) is set.

- `"lasso"`, `"ridge"`, `"elasticnet"`:

  The propensity scores are estimated using a lasso, ridge, or elastic
  net model, respectively. The `formula` supplied to
  [`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
  is processed with
  [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) and
  passed to
  [`glmnet::cv.glmnet()`](https://rdrr.io/pkg/glmnet/man/cv.glmnet.html),
  and
  [`glmnet::predict.cv.glmnet()`](https://rdrr.io/pkg/glmnet/man/predict.cv.glmnet.html)
  is used to compute the propensity scores. The `link` argument can be
  specified as a link function supplied to
  [`binomial()`](https://rdrr.io/r/stats/family.html), e.g., `"logit"`,
  which is the default. When `link` is prepended by `"linear."`, the
  linear predictor is used instead of the predicted probabilities. When
  `link = "log"`, a Poisson model is used. For
  `distance = "elasticnet"`, the `alpha` argument, which controls how to
  prioritize the lasso and ridge penalties in the elastic net, is set to
  .5 by default and can be changed by supplying an argument to `alpha`
  in `distance.options`. For `"lasso"` and `"ridge"`, `alpha` is set to
  1 and 0, respectively, and cannot be changed. The `cv.glmnet()`
  defaults are used to select the tuning parameters and generate
  predictions and can be modified using `distance.options`. If the `s`
  argument is passed to `distance.options`, it will be passed to
  `predict.cv.glmnet()`. Note that because there is a random component
  to choosing the tuning parameter, results will vary across runs unless
  a [seed](https://rdrr.io/r/base/Random.html) is set.

- `"rpart"`:

  The propensity scores are estimated using a classification tree. The
  `formula` supplied to
  [`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
  is passed directly to
  [`rpart::rpart()`](https://rdrr.io/pkg/rpart/man/rpart.html), and
  [`rpart::predict.rpart()`](https://rdrr.io/pkg/rpart/man/predict.rpart.html)
  is used to compute the propensity scores. The `link` argument is
  ignored, and predicted probabilities are always returned as the
  distance measure.

- `"randomforest"`:

  The propensity scores are estimated using a random forest. The
  `formula` supplied to
  [`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
  is passed directly to
  [`randomForest::randomForest()`](https://rdrr.io/pkg/randomForest/man/randomForest.html),
  and
  [`randomForest::predict.randomForest()`](https://rdrr.io/pkg/randomForest/man/predict.randomForest.html)
  is used to compute the propensity scores. The `link` argument is
  ignored, and predicted probabilities are always returned as the
  distance measure. Note that because there is a random component,
  results will vary across runs unless a
  [seed](https://rdrr.io/r/base/Random.html) is set.

- `"nnet"`:

  The propensity scores are estimated using a single-hidden-layer neural
  network. The `formula` supplied to
  [`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
  is passed directly to
  [`nnet::nnet()`](https://rdrr.io/pkg/nnet/man/nnet.html), and
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) is used to
  compute the propensity scores. The `link` argument is ignored, and
  predicted probabilities are always returned as the distance measure.
  An argument to `size` must be supplied to `distance.options` when
  using `method = "nnet"`.

- `"cbps"`:

  The propensity scores are estimated using the covariate balancing
  propensity score (CBPS) algorithm, which is a form of logistic
  regression where balance constraints are incorporated to a generalized
  method of moments estimation of of the model coefficients. The
  `formula` supplied to
  [`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
  is passed directly to
  [`CBPS::CBPS()`](https://rdrr.io/pkg/CBPS/man/CBPS.html), and
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) is used to
  compute the propensity scores. The `link` argument can be specified as
  `"linear"` to use the linear predictor instead of the predicted
  probabilities. No other links are allowed. The `estimand` argument
  supplied to
  [`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
  will be used to select the appropriate estimand for use in defining
  the balance constraints, so no argument needs to be supplied to `ATT`
  in `CBPS`.

- `"bart"`:

  The propensity scores are estimated using Bayesian additive regression
  trees (BART). The `formula` supplied to
  [`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
  is passed directly to
  [`dbarts::bart2()`](https://rdrr.io/pkg/dbarts/man/bart.html), and
  [`dbarts::fitted.bart()`](https://rdrr.io/pkg/dbarts/man/bart.html) is
  used to compute the propensity scores. The `link` argument can be
  specified as `"linear"` to use the linear predictor instead of the
  predicted probabilities. When `s.weights` is supplied to
  [`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md),
  it will not be passed to `bart2` because the `weights` argument in
  `bart2` does not correspond to sampling weights. Note that because
  there is a random component to choosing the tuning parameter, results
  will vary across runs unless the `seed` argument is supplied to
  `distance.options`. Note that setting a seed using
  [`set.seed()`](https://rdrr.io/r/base/Random.html) is not sufficient
  to guarantee reproducibility unless single-threading is used. See
  [`dbarts::bart2()`](https://rdrr.io/pkg/dbarts/man/bart.html) for
  details.

### Methods for computing distances from covariates

The following methods involve computing a distance matrix from the
covariates themselves without estimating a propensity score. Calipers on
the distance measure and common support restrictions cannot be used, and
the `distance` component of the output object will be empty because no
propensity scores are estimated. The `link` and `distance.options`
arguments are ignored with these methods. See the individual matching
methods pages for whether these distances are allowed and how they are
used. Each of these distance measures can also be calculated outside
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
using its [corresponding
function](https://kosukeimai.github.io/MatchIt/reference/mahalanobis_dist.md).

- `"euclidean"`:

  The Euclidean distance is the raw distance between units, computed as
  \$\$d\_{ij} = \sqrt{(x_i - x_j)(x_i - x_j)'}\$\$ It is sensitive to
  the scale of the covariates, so covariates with larger scales will
  take higher priority.

- `"scaled_euclidean"`:

  The scaled Euclidean distance is the Euclidean distance computed on
  the scaled (i.e., standardized) covariates. This ensures the
  covariates are on the same scale. The covariates are standardized
  using the pooled within-group standard deviations, computed by
  treatment group-mean centering each covariate before computing the
  standard deviation in the full sample.

- `"mahalanobis"`:

  The Mahalanobis distance is computed as \$\$d\_{ij} = \sqrt{(x_i -
  x_j)\Sigma^{-1}(x_i - x_j)'}\$\$ where \\\Sigma\\ is the pooled
  within-group covariance matrix of the covariates, computed by
  treatment group-mean centering each covariate before computing the
  covariance in the full sample. This ensures the variables are on the
  same scale and accounts for the correlation between covariates.

- `"robust_mahalanobis"`:

  The robust rank-based Mahalanobis distance is the Mahalanobis distance
  computed on the ranks of the covariates with an adjustment for ties.
  It is described in Rosenbaum (2010, ch. 8) as an alternative to the
  Mahalanobis distance that handles outliers and rare categories better
  than the standard Mahalanobis distance but is not affinely invariant.

To perform Mahalanobis distance matching *and* estimate propensity
scores to be used for a purpose other than matching, the `mahvars`
argument should be used along with a different specification to
`distance`. See the individual matching method pages for details on how
to use `mahvars`.

### Distances supplied as a numeric vector or matrix

`distance` can also be supplied as a numeric vector whose values will be
taken to function like propensity scores; their pairwise difference will
define the distance between units. This might be useful for supplying
propensity scores computed outside
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
or resupplying
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
with propensity scores estimated previously without having to recompute
them.

`distance` can also be supplied as a matrix whose values represent the
pairwise distances between units. The matrix should either be a square,
with a row and column for each unit (e.g., as the output of a call to
`as.matrix([`dist`](.))`), or have as many rows as there are treated
units and as many columns as there are control units (e.g., as the
output of a call to
[`mahalanobis_dist()`](https://kosukeimai.github.io/MatchIt/reference/mahalanobis_dist.md)
or
[`optmatch::match_on()`](https://rdrr.io/pkg/optmatch/man/match_on-methods.html)).
Distance values of `Inf` will disallow the corresponding units to be
matched. When `distance` is a supplied as a numeric vector or matrix,
`link` and `distance.options` are ignored.

## Examples

``` r
data("lalonde")

# Matching on logit of a PS estimated with logistic
# regression:
m.out1 <- matchit(treat ~ age + educ + race + married +
                    nodegree + re74 + re75,
                  data = lalonde,
                  distance = "glm",
                  link = "linear.logit")
# GAM logistic PS with smoothing splines (s()):
m.out2 <- matchit(treat ~ s(age) + s(educ) +
                    race + married +
                    nodegree + re74 + re75,
                  data = lalonde,
                  distance = "gam")
summary(m.out2$model)
#> 
#> Family: quasibinomial 
#> Link function: logit 
#> 
#> Formula:
#> treat ~ s(age) + s(educ) + race + married + nodegree + re74 + 
#>     re75
#> 
#> Parametric coefficients:
#>               Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)  5.436e-01  3.950e-01   1.376  0.16928    
#> racehispan  -2.447e+00  4.323e-01  -5.661 2.34e-08 ***
#> racewhite   -2.995e+00  3.136e-01  -9.552  < 2e-16 ***
#> married     -1.644e+00  3.438e-01  -4.781 2.20e-06 ***
#> nodegree     7.894e-01  4.800e-01   1.645  0.10058    
#> re74        -9.838e-05  3.245e-05  -3.031  0.00254 ** 
#> re75         5.113e-05  5.001e-05   1.022  0.30702    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Approximate significance of smooth terms:
#>           edf Ref.df     F p-value    
#> s(age)  7.489  8.144 6.781  <2e-16 ***
#> s(educ) 2.647  3.359 2.311  0.0628 .  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> R-sq.(adj) =    0.5   Deviance explained = 46.1%
#> GCV = 0.69813  Scale est. = 1.0287    n = 614
# CBPS for ATC matching w/replacement, using the just-
# identified version of CBPS (setting method = "exact"):
m.out3 <- matchit(treat ~ age + educ + race + married +
                    nodegree + re74 + re75,
                  data = lalonde,
                  distance = "cbps",
                  estimand = "ATC",
                  distance.options = list(method = "exact"),
                  replace = TRUE)
# Mahalanobis distance matching - no PS estimated
m.out4 <- matchit(treat ~ age + educ + race + married +
                    nodegree + re74 + re75,
                  data = lalonde,
                  distance = "mahalanobis")

m.out4$distance #NULL
#> NULL

# Mahalanobis distance matching with PS estimated
# for use in a caliper; matching done on mahvars
m.out5 <- matchit(treat ~ age + educ + race + married +
                    nodegree + re74 + re75,
                  data = lalonde,
                  distance = "glm",
                  caliper = .1,
                  mahvars = ~ age + educ + race + married +
                                nodegree + re74 + re75)

summary(m.out5)
#> 
#> Call:
#> matchit(formula = treat ~ age + educ + race + married + nodegree + 
#>     re74 + re75, data = lalonde, distance = "glm", mahvars = ~age + 
#>     educ + race + married + nodegree + re74 + re75, caliper = 0.1)
#> 
#> Summary of Balance for All Data:
#>            Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean
#> distance          0.5774        0.1822          1.7941     0.9211    0.3774
#> age              25.8162       28.0303         -0.3094     0.4400    0.0813
#> educ             10.3459       10.2354          0.0550     0.4959    0.0347
#> raceblack         0.8432        0.2028          1.7615          .    0.6404
#> racehispan        0.0595        0.1422         -0.3498          .    0.0827
#> racewhite         0.0973        0.6550         -1.8819          .    0.5577
#> married           0.1892        0.5128         -0.8263          .    0.3236
#> nodegree          0.7081        0.5967          0.2450          .    0.1114
#> re74           2095.5737     5619.2365         -0.7211     0.5181    0.2248
#> re75           1532.0553     2466.4844         -0.2903     0.9563    0.1342
#>            eCDF Max
#> distance     0.6444
#> age          0.1577
#> educ         0.1114
#> raceblack    0.6404
#> racehispan   0.0827
#> racewhite    0.5577
#> married      0.3236
#> nodegree     0.1114
#> re74         0.4470
#> re75         0.2876
#> 
#> Summary of Balance for Matched Data:
#>            Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean
#> distance          0.5096        0.4916          0.0817     1.0782    0.0225
#> age              25.9459       25.3874          0.0781     0.4201    0.0869
#> educ             10.4865       10.2883          0.0986     0.6448    0.0199
#> raceblack         0.7387        0.7207          0.0496          .    0.0180
#> racehispan        0.0991        0.1171         -0.0762          .    0.0180
#> racewhite         0.1622        0.1622          0.0000          .    0.0000
#> married           0.2072        0.2162         -0.0230          .    0.0090
#> nodegree          0.6486        0.6306          0.0396          .    0.0180
#> re74           2667.1135     2357.0686          0.0634     1.7471    0.0387
#> re75           1811.2988     1506.3709          0.0947     1.9333    0.0253
#>            eCDF Max Std. Pair Dist.
#> distance     0.1441          0.0922
#> age          0.3063          0.9922
#> educ         0.0811          0.7707
#> raceblack    0.0180          0.0496
#> racehispan   0.0180          0.0762
#> racewhite    0.0000          0.0000
#> married      0.0090          0.4370
#> nodegree     0.0180          0.4756
#> re74         0.2432          0.5712
#> re75         0.0901          0.5442
#> 
#> Sample Sizes:
#>           Control Treated
#> All           429     185
#> Matched       111     111
#> Unmatched     318      74
#> Discarded       0       0
#> 

# User-supplied propensity scores
p.score <- fitted(glm(treat ~ age + educ + race + married +
                        nodegree + re74 + re75,
                      data = lalonde,
                      family = binomial))

m.out6 <- matchit(treat ~ age + educ + race + married +
                    nodegree + re74 + re75,
                  data = lalonde,
                  distance = p.score)

# User-supplied distance matrix using rank_mahalanobis()
dist_mat <- robust_mahalanobis_dist(
              treat ~ age + educ + race + nodegree +
                married + re74 + re75,
              data = lalonde)

m.out7 <- matchit(treat ~ age + educ + race + nodegree +
                    married + re74 + re75,
                  data = lalonde,
                  distance = dist_mat)
```
