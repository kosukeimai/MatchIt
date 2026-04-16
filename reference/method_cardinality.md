# Cardinality Matching

In
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md),
setting `method = "cardinality"` performs cardinality matching and other
forms of matching that use mixed integer programming. Rather than
forming pairs, cardinality matching selects the largest subset of units
that satisfies user-supplied balance constraints on mean differences.
One of several available optimization programs can be used to solve the
mixed integer program. The default is the HiGHS library as implemented
in the *highs* package, both of which are free, but performance can be
improved using Gurobi and the *gurobi* package, for which there is a
free academic license.

This page details the allowable arguments with `method = "cardinality"`.
See
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
for an explanation of what each argument means in a general context and
how it can be specified.

Below is how
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
is used for cardinality matching:

    matchit(formula,
            data = NULL,
            method = "cardinality",
            estimand = "ATT",
            exact = NULL,
            mahvars = NULL,
            s.weights = NULL,
            ratio = 1,
            verbose = FALSE,
            tols = .05,
            std.tols = TRUE,
            solver = "highs",
            ...) 

## Arguments

- formula:

  a two-sided [formula](https://rdrr.io/r/stats/formula.html) object
  containing the treatment and covariates to be balanced.

- data:

  a data frame containing the variables named in `formula`. If not found
  in `data`, the variables will be sought in the environment.

- method:

  set here to `"cardinality"`.

- estimand:

  a string containing the desired estimand. Allowable options include
  `"ATT"`, `"ATC"`, and `"ATE"`. See Details.

- exact:

  for which variables exact matching should take place. Separate
  optimization will occur within each subgroup of the exact matching
  variables.

- mahvars:

  which variables should be used for pairing after subset selection. Can
  only be set when `ratio` is a whole number. See Details.

- s.weights:

  the variable containing sampling weights to be incorporated into the
  optimization. The balance constraints refer to the product of the
  sampling weights and the matching weights, and the sum of the product
  of the sampling and matching weights will be maximized.

- ratio:

  the desired ratio of control to treated units. Can be set to `NA` to
  maximize sample size without concern for this ratio. See Details.

- verbose:

  `logical`; whether information about the matching process should be
  printed to the console.

- ...:

  additional arguments that control the matching specification:

  `tols`

  :   `numeric`; a vector of imbalance tolerances for mean differences,
      one for each covariate in `formula`. If only one value is
      supplied, it is applied to all. See `std.tols` below. Default is
      `.05` for standardized mean differences of at most .05 for all
      covariates between the treatment groups in the matched sample.

  `std.tols`

  :   `logical`; whether each entry in `tols` corresponds to a raw or
      standardized mean difference. If only one value is supplied, it is
      applied to all. Default is `TRUE` for standardized mean
      differences. The standardization factor is the pooled standard
      deviation when `estimand = "ATE"`, the standard deviation of the
      treated group when `estimand = "ATT"`, and the standard deviation
      of the control group when `estimand = "ATC"` (the same as used in
      [`summary.matchit()`](https://kosukeimai.github.io/MatchIt/reference/summary.matchit.md)).

  `solver`

  :   the name of solver to use to solve the optimization problem.
      Available options include `"highs"`, `"glpk"`, `"symphony"`, and
      `"gurobi"` for HiGHS (implemented in the *highs* package), GLPK
      (implemented in the *Rglpk* package), SYMPHONY (implemented in the
      *Rsymphony* package), and Gurobi (implemented in the *gurobi*
      package), respectively. The differences between them are in speed
      and solving ability. HiGHS (the default) and GLPK are the easiest
      to install, but Gurobi is recommended as it consistently
      outperforms other solvers and can find solutions even when others
      can't, and in less time. Gurobi is proprietary but can be used
      with a free trial or academic license. SYMPHONY may not produce
      reproducible results, even with a seed set.

  `time`

  :   the maximum amount of time before the optimization routine aborts,
      in seconds. Default is 120 (2 minutes). For large problems, this
      should be set much higher.

  The arguments `distance` (and related arguments), `replace`,
  `m.order`, and `caliper` (and related arguments) are ignored with a
  warning.

## Details

### Cardinality and Profile Matching

Two types of matching are available with `method = "cardinality"`:
cardinality matching and profile matching.

**Cardinality matching** finds the largest matched set that satisfies
the balance constraints between treatment groups, with the additional
constraint that the ratio of the number of matched control to matched
treated units is equal to `ratio` (1 by default), mimicking k:1
matching. When not all treated units are included in the matched set,
the estimand no longer corresponds to the ATT, so cardinality matching
should be avoided if retaining the ATT is desired. To request
cardinality matching, `estimand` should be set to `"ATT"` or `"ATC"` and
`ratio` should be set to a positive integer. 1:1 cardinality matching is
the default method when no arguments are specified.

**Profile matching** finds the largest matched set that satisfies
balance constraints between each treatment group and a specified target
sample. When `estimand = "ATT"`, it will find the largest subset of the
control units that satisfies the balance constraints with respect to the
treated group, which is left intact. When `estimand = "ATE"`, it will
find the largest subsets of the treated group and of the control group
that are balanced to the overall sample. To request profile matching for
the ATT, `estimand` should be set to `"ATT"` and `ratio` to `NA`. To
request profile matching for the ATE, `estimand` should be set to
`"ATE"` and `ratio` can be set either to `NA` to maximize the size of
each sample independently or to a positive integer to ensure that the
ratio of matched control units to matched treated treats is fixed,
mimicking k:1 matching. Unlike cardinality matching, profile matching
retains the requested estimand if a solution is found.

Neither method involves creating pairs in the matched set, but it is
possible to perform an additional round of pairing within the matched
sample after cardinality matching or profile matching for the ATE with a
fixed whole number sample size ratio by supplying the desired pairing
variables to `mahvars`. Doing so will trigger [optimal
matching](https://kosukeimai.github.io/MatchIt/reference/method_optimal.md)
using
[`optmatch::pairmatch()`](https://rdrr.io/pkg/optmatch/man/pairmatch.html)
on the Mahalanobis distance computed using the variables supplied to
`mahvars`. The balance or composition of the matched sample will not
change, but additional precision and robustness can be gained by forming
the pairs.

The weights are scaled so that the sum of the weights in each group is
equal to the number of matched units in the smaller group when
cardinality matching or profile matching for the ATE, and scaled so that
the sum of the weights in the control group is equal to the number of
treated units when profile matching for the ATT. When the sample sizes
of the matched groups is the same (i.e., when `ratio = 1`), no scaling
is done. Robust standard errors should be used in effect estimation
after cardinality or profile matching (and cluster-robust standard
errors if additional pairing is done in the matched sample). See
[`vignette("estimating-effects")`](https://kosukeimai.github.io/MatchIt/articles/estimating-effects.md)
for more information.

### Specifying Balance Constraints

The balance constraints are on the (standardized) mean differences
between the matched treatment groups for each covariate. Balance
constraints should be set by supplying arguments to `tols` and
`std.tols`. For example, setting `tols = .1` and `std.tols = TRUE`
requests that all the mean differences in the matched sample should be
within .1 standard deviations for each covariate. Different tolerances
can be set for different variables; it might be beneficial to constrain
the mean differences for highly prognostic covariates more tightly than
for other variables. For example, one could specify
`tols = c(.001, .05), std.tols = c(TRUE, FALSE)` to request that the
standardized mean difference for the first covariate is less than .001
and the raw mean difference for the second covariate is less than .05.
The values should be specified in the order they appear in `formula`,
except when interactions are present. One can run the following code:

    MatchIt:::get_assign(model.matrix(~X1*X2 + X3, data = data))[-1]

which will output a vector of numbers and the variable to which each
number corresponds; the first entry in `tols` corresponds to the
variable labeled 1, the second to the variable labeled 2, etc.

### Dealing with Errors and Warnings

When the optimization cannot be solved at all, or at least within the
time frame specified in the argument to `time`, an error or warning will
appear. Unfortunately, it is hard to know exactly the cause of the
failure and what measures should be taken to rectify it.

A warning that says
`"The optimizer failed to find an optimal solution in the time alotted. The returned solution may not be optimal."`
usually means that an optimal solution may be possible to find with more
time, in which case `time` should be increased or a faster solver should
be used. Even with this warning, a potentially usable solution will be
returned, so don't automatically take it to mean the optimization
failed. Sometimes, when there are multiple solutions with the same
resulting sample size, the optimizers will stall at one of them, not
thinking it has found the optimum. The result should be checked to see
if it can be used as the solution.

An error that says `"The optimization problem may be infeasible."`
usually means that there is a issue with the optimization problem, i.e.,
that there is no possible way to satisfy the constraints. To rectify
this, one can try relaxing the constraints by increasing the value of
`tols` or use another solver. Sometimes Gurobi can solve problems that
the other solvers cannot.

## Outputs

Most outputs described in
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
are returned with `method = "cardinality"`. Unless `mahvars` is
specified, the `match.matrix` and `subclass` components are omitted
because no pairing or subclassification is done. When
`include.obj = TRUE` in the call to
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md),
the output of the optimization function will be included in the output.
When `exact` is specified, this will be a list of such objects, one for
each stratum of the exact variables.

## References

In a manuscript, you should reference the solver used in the
optimization. For example, a sentence might read:

*Cardinality matching was performed using the MatchIt package (Ho, Imai,
King, & Stuart, 2011) in R with the optimization performed by HiGHS
(Huangfu & Hall, 2018).*

See
[`vignette("matching-methods")`](https://kosukeimai.github.io/MatchIt/articles/matching-methods.md)
for more literature on cardinality matching.

## See also

[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
for a detailed explanation of the inputs and outputs of a call to
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md).

*[designmatch](https://CRAN.R-project.org/package=designmatch)*, which
performs cardinality and profile matching with many more options and
more flexibility. The implementations of cardinality matching differ
between *MatchIt* and *designmatch*, so their results might differ.

*[optweight](https://CRAN.R-project.org/package=optweight)*, which
offers similar functionality but in the context of weighting rather than
matching.

## Examples

``` r
data("lalonde")

#Choose your solver; "gurobi" is best, "highs" is free and
#easy to install
solver <- "highs"

m.out1 <- matchit(treat ~ age + educ + re74,
                  data = lalonde,
                  method = "cardinality",
                  estimand = "ATT",
                  ratio = 1,
                  tols = .2,
                  solver = solver)
m.out1
#> A `matchit` object
#>  - method: Cardinality matching
#>  - number of obs.: 614 (original), 370 (matched)
#>  - target estimand: ATT
#>  - covariates: age, educ, re74
summary(m.out1)
#> 
#> Call:
#> matchit(formula = treat ~ age + educ + re74, data = lalonde, 
#>     method = "cardinality", estimand = "ATT", ratio = 1, tols = 0.2, 
#>     solver = solver)
#> 
#> Summary of Balance for All Data:
#>      Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean eCDF Max
#> age        25.8162       28.0303         -0.3094     0.4400    0.0813   0.1577
#> educ       10.3459       10.2354          0.0550     0.4959    0.0347   0.1114
#> re74     2095.5737     5619.2365         -0.7211     0.5181    0.2248   0.4470
#> 
#> Summary of Balance for Matched Data:
#>      Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean eCDF Max
#> age        25.8162       27.2108         -0.1949     0.4858    0.0722   0.1459
#> educ       10.3459       10.3459          0.0000     0.5470    0.0284   0.0919
#> re74     2095.5737     3051.6162         -0.1956     1.2100    0.1026   0.3243
#> 
#> Sample Sizes:
#>           Control Treated
#> All           429     185
#> Matched       185     185
#> Unmatched     244       0
#> Discarded       0       0
#> 

# Profile matching for the ATT
m.out2 <- matchit(treat ~ age + educ + re74,
                  data = lalonde,
                  method = "cardinality",
                  estimand = "ATT",
                  ratio = NA,
                  tols = .2,
                  solver = solver)
m.out2
#> A `matchit` object
#>  - method: Cardinality matching
#>  - number of obs.: 614 (original), 540 (matched)
#>  - target estimand: ATT
#>  - covariates: age, educ, re74
summary(m.out2, un = FALSE)
#> 
#> Call:
#> matchit(formula = treat ~ age + educ + re74, data = lalonde, 
#>     method = "cardinality", estimand = "ATT", ratio = NA, tols = 0.2, 
#>     solver = solver)
#> 
#> Summary of Balance for Matched Data:
#>      Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean eCDF Max
#> age        25.8162       26.3662         -0.0769     0.4781    0.0695   0.1593
#> educ       10.3459       10.0761          0.1342     0.4756    0.0375   0.1106
#> re74     2095.5737     3069.7682         -0.1994     1.7132    0.1305   0.3926
#> 
#> Sample Sizes:
#>           Control Treated
#> All           429     185
#> Matched       355     185
#> Unmatched      74       0
#> Discarded       0       0
#> 

# Profile matching for the ATE
m.out3 <- matchit(treat ~ age + educ + re74,
                  data = lalonde,
                  method = "cardinality",
                  estimand = "ATE",
                  ratio = NA,
                  tols = .2,
                  solver = solver)
m.out3
#> A `matchit` object
#>  - method: Cardinality matching
#>  - number of obs.: 614 (original), 516 (matched)
#>  - target estimand: ATE
#>  - covariates: age, educ, re74
summary(m.out3, un = FALSE)
#> 
#> Call:
#> matchit(formula = treat ~ age + educ + re74, data = lalonde, 
#>     method = "cardinality", estimand = "ATE", ratio = NA, tols = 0.2, 
#>     solver = solver)
#> 
#> Summary of Balance for Matched Data:
#>      Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean eCDF Max
#> age        27.1753       27.7852         -0.0666     0.5914    0.0573   0.1154
#> educ       10.1856       10.1718          0.0056     0.6960    0.0212   0.1051
#> re74     3971.1184     5146.7881         -0.1988     1.0186    0.0831   0.2069
#> 
#> Sample Sizes:
#>           Control Treated
#> All           429     185
#> Matched       419      97
#> Unmatched      10      88
#> Discarded       0       0
#> 
# Pairing after 1:1 cardinality matching:
m.out1b <- matchit(treat ~ age + educ + re74,
                   data = lalonde,
                   method = "cardinality",
                   estimand = "ATT",
                   ratio = 1,
                   tols = .15,
                   solver = solver,
                   mahvars = ~ age + educ + re74)

# Note that balance doesn't change but pair distances
# are lower for the paired-upon variables
summary(m.out1b, un = FALSE)
#> 
#> Call:
#> matchit(formula = treat ~ age + educ + re74, data = lalonde, 
#>     method = "cardinality", estimand = "ATT", mahvars = ~age + 
#>         educ + re74, ratio = 1, tols = 0.15, solver = solver)
#> 
#> Summary of Balance for Matched Data:
#>      Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean eCDF Max
#> age        25.8162       26.8432         -0.1435     0.5018    0.0673   0.1351
#> educ       10.3459       10.3514         -0.0027     0.5453    0.0287   0.0919
#> re74     2095.5737     2827.7621         -0.1498     1.3037    0.0890   0.3189
#>      Std. Pair Dist.
#> age           0.4487
#> educ          0.3576
#> re74          0.2796
#> 
#> Sample Sizes:
#>           Control Treated
#> All           429     185
#> Matched       185     185
#> Unmatched     244       0
#> Discarded       0       0
#> 
summary(m.out1, un = FALSE)
#> 
#> Call:
#> matchit(formula = treat ~ age + educ + re74, data = lalonde, 
#>     method = "cardinality", estimand = "ATT", ratio = 1, tols = 0.2, 
#>     solver = solver)
#> 
#> Summary of Balance for Matched Data:
#>      Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean eCDF Max
#> age        25.8162       27.2108         -0.1949     0.4858    0.0722   0.1459
#> educ       10.3459       10.3459          0.0000     0.5470    0.0284   0.0919
#> re74     2095.5737     3051.6162         -0.1956     1.2100    0.1026   0.3243
#> 
#> Sample Sizes:
#>           Control Treated
#> All           429     185
#> Matched       185     185
#> Unmatched     244       0
#> Discarded       0       0
#> 

# In these examples, a high tol was used and
# few covariate matched on in order to not take too long;
# with real data, tols should be much lower and more
# covariates included if possible.
```
