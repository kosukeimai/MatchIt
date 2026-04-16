# Coarsened Exact Matching

In
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md),
setting `method = "cem"` performs coarsened exact matching. With
coarsened exact matching, covariates are coarsened into bins, and a
complete cross of the coarsened covariates is used to form subclasses
defined by each combination of the coarsened covariate levels. Any
subclass that doesn't contain both treated and control units is
discarded, leaving only subclasses containing treatment and control
units that are exactly equal on the coarsened covariates. The coarsening
process can be controlled by an algorithm or by manually specifying
cutpoints and groupings. The benefits of coarsened exact matching are
that the tradeoff between exact matching and approximate balancing can
be managed to prevent discarding too many units, which can otherwise
occur with exact matching.

This page details the allowable arguments with `method = "cem"`. See
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
for an explanation of what each argument means in a general context and
how it can be specified.

Below is how
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
is used for coarsened exact matching:

    matchit(formula,
            data = NULL,
            method = "cem",
            estimand = "ATT",
            s.weights = NULL,
            verbose = FALSE,
            ...) 

## Arguments

- formula:

  a two-sided [formula](https://rdrr.io/r/stats/formula.html) object
  containing the treatment and covariates to be used in creating the
  subclasses defined by a full cross of the coarsened covariate levels.

- data:

  a data frame containing the variables named in `formula`. If not found
  in `data`, the variables will be sought in the environment.

- method:

  set here to `"cem"`.

- estimand:

  a string containing the desired estimand. Allowable options include
  `"ATT"`, `"ATC"`, and `"ATE"`. The estimand controls how the weights
  are computed; see the Computing Weights section at
  [`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
  for details. When `k2k = TRUE` (see below), `estimand` also controls
  how the matching is done.

- s.weights:

  the variable containing sampling weights to be incorporated into
  balance statistics or the scaling factors when `k2k = TRUE` and
  certain methods are used.

- verbose:

  `logical`; whether information about the matching process should be
  printed to the console.

- ...:

  additional arguments to control the matching process.

  `grouping`

  :   a named list with an (optional) entry for each categorical
      variable to be matched on. Each element should itself be a list,
      and each entry of the sublist should be a vector containing levels
      of the variable that should be combined to form a single level.
      Any categorical variables not included in `grouping` will remain
      as they are in the data, which means exact matching, with no
      coarsening, will take place on these variables. See Details.

  `cutpoints`

  :   a named list with an (optional) entry for each numeric variable to
      be matched on. Each element describes a way of coarsening the
      corresponding variable. They can be a vector of cutpoints that
      demarcate bins, a single number giving the number of bins, or a
      string corresponding to a method of computing the number of bins.
      Allowable strings include `"sturges"`, `"scott"`, and `"fd"`,
      which use the functions
      [`grDevices::nclass.Sturges()`](https://rdrr.io/r/grDevices/nclass.html),
      [`grDevices::nclass.scott()`](https://rdrr.io/r/grDevices/nclass.html),
      and
      [`grDevices::nclass.FD()`](https://rdrr.io/r/grDevices/nclass.html),
      respectively. The default is `"sturges"` for variables that are
      not listed or if no argument is supplied. Can also be a single
      value to be applied to all numeric variables. See Details.

  `k2k`

  :   `logical`; whether 1:1 matching should occur within the matched
      strata. If `TRUE` nearest neighbor matching without replacement
      will take place within each stratum, and any unmatched units will
      be dropped (e.g., if there are more treated than control units in
      the stratum, the treated units without a match will be dropped).
      The `k2k.method` argument controls how the distance between units
      is calculated.

  `k2k.method`

  :   `character`; how the distance between units should be calculated
      if `k2k = TRUE`. Allowable arguments include `NULL` (for random
      matching), any argument to
      [`distance()`](https://kosukeimai.github.io/MatchIt/reference/distance.md)
      for computing a distance matrix from covariates (e.g.,
      `"mahalanobis"`), or any allowable argument to `method` in
      [`dist()`](https://rdrr.io/r/stats/dist.html). Matching will take
      place on the original (non-coarsened) variables. The default is
      `"mahalanobis"`.

  `mpower`

  :   if `k2k.method = "minkowski"`, the power used in creating the
      distance. This is passed to the `p` argument of
      [`dist()`](https://rdrr.io/r/stats/dist.html).

  `m.order`

  :   `character`; the order that the matching takes place when
      `k2k = TRUE`. Allowable options include `"closest"`, where
      matching takes place in ascending order of the smallest distance
      between units; `"farthest"`, where matching takes place in
      descending order of the smallest distance between units;
      `"random"`, where matching takes place in a random order; and
      `"data"` where matching takes place based on the order of units in
      the data. When `m.order = "random"`, results may differ across
      different runs of the same code unless a seed is set and specified
      with [`set.seed()`](https://rdrr.io/r/base/Random.html). The
      default of `NULL` corresponds to `"data"`. See
      [`method_nearest`](https://kosukeimai.github.io/MatchIt/reference/method_nearest.md)
      for more information.

  The arguments `distance` (and related arguments), `exact`, `mahvars`,
  `discard` (and related arguments), `replace`, `caliper` (and related
  arguments), and `ratio` are ignored with a warning.

## Details

If the coarsening is such that there are no exact matches with the
coarsened variables, the `grouping` and `cutpoints` arguments can be
used to modify the matching specification. Reducing the number of
cutpoints or grouping some variable values together can make it easier
to find matches. See Examples below. Removing variables can also help
(but they will likely not be balanced unless highly correlated with the
included variables). To take advantage of coarsened exact matching
without failing to find any matches, the covariates can be manually
coarsened outside of
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
and then supplied to the `exact` argument in a call to
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
with another matching method.

Setting `k2k = TRUE` is equivalent to first doing coarsened exact
matching with `k2k = FALSE` and then supplying stratum membership as an
exact matching variable (i.e., in `exact`) to another call to
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
with `method = "nearest"`. It is also equivalent to performing nearest
neighbor matching supplying coarsened versions of the variables to
`exact`, except that `method = "cem"` automatically coarsens the
continuous variables. The `estimand` argument supplied with
`method = "cem"` functions the same way it would in these alternate
matching calls, i.e., by determining the "focal" group that controls the
order of the matching.

### Grouping and Cutpoints

The `grouping` and `cutpoints` arguments allow one to fine-tune the
coarsening of the covariates. `grouping` is used for combining
categories of categorical covariates and `cutpoints` is used for binning
numeric covariates. The values supplied to these arguments should be
iteratively changed until a matching solution that balances covariate
balance and remaining sample size is obtained. The arguments are
described below.

#### `grouping`

The argument to `grouping` must be a list, where each component has the
name of a categorical variable, the levels of which are to be combined.
Each component must itself be a list; this list contains one or more
vectors of levels, where each vector corresponds to the levels that
should be combined into a single category. For example, if a variable
`amount` had levels `"none"`, `"some"`, and `"a lot"`, one could enter
`grouping = list(amount = list(c("none"), c("some", "a lot")))`, which
would group `"some"` and `"a lot"` into a single category and leave
`"none"` in its own category. Any levels left out of the list for each
variable will be left alone (so `c("none")` could have been omitted from
the previous code). Note that if a categorical variable does not appear
in `grouping`, it will not be coarsened, so exact matching will take
place on it. `grouping` should not be used for numeric variables with
more than a few values; use `cutpoints`, described below, instead.

#### `cutpoints`

The argument to `cutpoints` must also be a list, where each component
has the name of a numeric variables that is to be binned. (As a
shortcut, it can also be a single value that will be applied to all
numeric variables). Each component can take one of three forms: a vector
of cutpoints that separate the bins, a single number giving the number
of bins, or a string corresponding to an algorithm used to compute the
number of bins. Any values at a boundary will be placed into the higher
bin; e.g., if the cutpoints were `c(0, 5, 10)`, values of 5 would be
placed into the same bin as values of 6, 7, 8, or 9, and values of 10
would be placed into a different bin. Internally, values of `-Inf` and
`Inf` are appended to the beginning and end of the range. When given as
a single number defining the number of bins, the bin boundaries are the
maximum and minimum values of the variable with bin boundaries evenly
spaced between them, i.e., not quantiles. A value of 0 will not perform
any binning (equivalent to exact matching on the variable), and a value
of 1 will remove the variable from the exact matching variables but it
will be still used for pair matching when `k2k = TRUE`. The allowable
strings include `"sturges"`, `"scott"`, and `"fd"`, which use the
corresponding binning method, and `"q#"` where `#` is a number, which
splits the variable into `#` equally-sized bins (i.e., quantiles).

An example of a way to supply an argument to `cutpoints` would be the
following:

    cutpoints = list(X1 = 4,
                     X2 = c(1.7, 5.5, 10.2),
                     X3 = "scott",
                     X4 = "q5") 

This would split `X1` into 4 bins, `X2` into bins based on the provided
boundaries, `X3` into a number of bins determined by
[`grDevices::nclass.scott()`](https://rdrr.io/r/grDevices/nclass.html),
and `X4` into quintiles. All other numeric variables would be split into
a number of bins determined by
[`grDevices::nclass.Sturges()`](https://rdrr.io/r/grDevices/nclass.html),
the default.

## Note

This method does not rely on the *cem* package, instead using code
written for *MatchIt*, but its design is based on the original *cem*
functions. Versions of *MatchIt* prior to 4.1.0 did rely on *cem*, so
results may differ between versions. There are a few differences between
the ways *MatchIt* and *cem* (and older versions of *MatchIt*) differ in
executing coarsened exact matching, described below.

- In *MatchIt*, when a single number is supplied to `cutpoints`, it
  describes the number of bins; in *cem*, it describes the number of
  cutpoints separating bins. The *MatchIt* method is closer to how
  [`hist()`](https://rdrr.io/r/graphics/hist.html) processes breaks
  points to create bins.

- In *MatchIt*, values on the cutpoint boundaries will be placed into
  the higher bin; in *cem*, they are placed into the lower bin. To avoid
  consequences of this choice, ensure the bin boundaries do not coincide
  with observed values of the variables.

- When `cutpoints` are used, `"ss"` (for Shimazaki-Shinomoto's rule) can
  be used in *cem* but not in *MatchIt*.

- When `k2k = TRUE`, *MatchIt* matches on the original variables
  (scaled), whereas *cem* matches on the coarsened variables. Because
  the variables are already exactly matched on the coarsened variables,
  matching in *cem* is equivalent to random matching within strata.

- When `k2k = TRUE`, in *MatchIt* matched units are identified by pair
  membership, and the original stratum membership prior to 1:1 matching
  is discarded. In *cem*, pairs are not identified beyond the stratum
  the members are part of.

- When `k2k = TRUE`, `k2k.method = "mahalanobis"` can be requested in
  *MatchIt* but not in *cem*.

## Outputs

All outputs described in
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
are returned with `method = "cem"` except for `match.matrix`. When
`k2k = TRUE`, a `match.matrix` component with the matched pairs is also
included. `include.obj` is ignored.

## References

In a manuscript, you don't need to cite another package when using
`method = "cem"` because the matching is performed completely within
*MatchIt*. For example, a sentence might read:

*Coarsened exact matching was performed using the MatchIt package (Ho,
Imai, King, & Stuart, 2011) in R.*

It would be a good idea to cite the following article, which develops
the theory behind coarsened exact matching:

Iacus, S. M., King, G., & Porro, G. (2012). Causal Inference without
Balance Checking: Coarsened Exact Matching. *Political Analysis*, 20(1),
1–24. [doi:10.1093/pan/mpr013](https://doi.org/10.1093/pan/mpr013)

## See also

[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
for a detailed explanation of the inputs and outputs of a call to
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md).

The *cem* package, upon which this method is based and which provided
the workhorse in previous versions of *MatchIt*.

[`method_exact`](https://kosukeimai.github.io/MatchIt/reference/method_exact.md)
for exact matching, which performs exact matching on the covariates
without coarsening.

## Examples

``` r
data("lalonde")

# Coarsened exact matching on age, race, married, and educ with educ
# coarsened into 5 bins and race coarsened into 2 categories,
# grouping "white" and "hispan" together
cutpoints <- list(educ = 5)
grouping <- list(race = list(c("white", "hispan"),
                             c("black")))

m.out1 <- matchit(treat ~ age + race + married + educ,
                  data = lalonde,
                  method = "cem",
                  cutpoints = cutpoints,
                  grouping = grouping)
m.out1
#> A `matchit` object
#>  - method: Coarsened exact matching
#>  - number of obs.: 614 (original), 393 (matched)
#>  - target estimand: ATT
#>  - covariates: age, race, married, educ
summary(m.out1)
#> 
#> Call:
#> matchit(formula = treat ~ age + race + married + educ, data = lalonde, 
#>     method = "cem", cutpoints = cutpoints, grouping = grouping)
#> 
#> Summary of Balance for All Data:
#>            Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean
#> age              25.8162       28.0303         -0.3094     0.4400    0.0813
#> raceblack         0.8432        0.2028          1.7615          .    0.6404
#> racehispan        0.0595        0.1422         -0.3498          .    0.0827
#> racewhite         0.0973        0.6550         -1.8819          .    0.5577
#> married           0.1892        0.5128         -0.8263          .    0.3236
#> educ             10.3459       10.2354          0.0550     0.4959    0.0347
#>            eCDF Max
#> age          0.1577
#> raceblack    0.6404
#> racehispan   0.0827
#> racewhite    0.5577
#> married      0.3236
#> educ         0.1114
#> 
#> Summary of Balance for Matched Data:
#>            Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean
#> age              24.1007       23.7091          0.0547     0.8882    0.0140
#> raceblack         0.8188        0.8188          0.0000          .    0.0000
#> racehispan        0.0671        0.0360          0.1317          .    0.0312
#> racewhite         0.1141        0.1452         -0.1051          .    0.0312
#> married           0.1745        0.1745         -0.0000          .    0.0000
#> educ             10.5369       10.8649         -0.1631     0.7974    0.0177
#>            eCDF Max Std. Pair Dist.
#> age          0.0897          0.1555
#> raceblack    0.0000          0.0000
#> racehispan   0.0312          0.4829
#> racewhite    0.0312          0.3853
#> married      0.0000          0.0000
#> educ         0.1865          0.3962
#> 
#> Sample Sizes:
#>               Control Treated
#> All            429.       185
#> Matched (ESS)   52.84     149
#> Matched        244.       149
#> Unmatched      185.        36
#> Discarded        0.         0
#> 

# The same but requesting 1:1 Mahalanobis distance matching with
# the k2k and k2k.method argument. Note the remaining number of units
# is smaller than when retaining the full matched sample.
m.out2 <- matchit(treat ~ age + race + married + educ,
                  data = lalonde,
                  method = "cem",
                  cutpoints = cutpoints,
                  grouping = grouping,
                  k2k = TRUE,
                  k2k.method = "mahalanobis")
m.out2
#> A `matchit` object
#>  - method: Coarsened exact matching
#>  - number of obs.: 614 (original), 170 (matched)
#>  - target estimand: ATT
#>  - covariates: age, race, married, educ
summary(m.out2, un = FALSE)
#> 
#> Call:
#> matchit(formula = treat ~ age + race + married + educ, data = lalonde, 
#>     method = "cem", cutpoints = cutpoints, grouping = grouping, 
#>     k2k = TRUE, k2k.method = "mahalanobis")
#> 
#> Summary of Balance for Matched Data:
#>            Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean
#> age              23.6353       23.4118          0.0312     0.9022    0.0126
#> raceblack         0.6824        0.6824          0.0000          .    0.0000
#> racehispan        0.1176        0.0706          0.1990          .    0.0471
#> racewhite         0.2000        0.2471         -0.1588          .    0.0471
#> married           0.2118        0.2118          0.0000          .    0.0000
#> educ             10.2588       10.5176         -0.1287     0.8036    0.0149
#>            eCDF Max Std. Pair Dist.
#> age          0.0824          0.1036
#> raceblack    0.0000          0.0000
#> racehispan   0.0471          0.1990
#> racewhite    0.0471          0.1588
#> married      0.0000          0.0000
#> educ         0.1294          0.2106
#> 
#> Sample Sizes:
#>           Control Treated
#> All           429     185
#> Matched        85      85
#> Unmatched     344     100
#> Discarded       0       0
#> 
```
