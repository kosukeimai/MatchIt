# Compute a Distance Matrix

The functions compute a distance matrix, either for a single dataset
(i.e., the distances between all pairs of units) or for two groups
defined by a splitting variable (i.e., the distances between all units
in one group and all units in the other). These distance matrices
include the Mahalanobis distance, Euclidean distance, scaled Euclidean
distance, and robust (rank-based) Mahalanobis distance. These functions
can be used as inputs to the `distance` argument to
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
and are used to compute the corresponding distance matrices within
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
when named.

## Usage

``` r
mahalanobis_dist(
  formula = NULL,
  data = NULL,
  s.weights = NULL,
  var = NULL,
  discarded = NULL,
  ...
)

scaled_euclidean_dist(
  formula = NULL,
  data = NULL,
  s.weights = NULL,
  var = NULL,
  discarded = NULL,
  ...
)

robust_mahalanobis_dist(
  formula = NULL,
  data = NULL,
  s.weights = NULL,
  discarded = NULL,
  ...
)

euclidean_dist(formula = NULL, data = NULL, ...)
```

## Arguments

- formula:

  a formula with the treatment (i.e., splitting variable) on the left
  side and the covariates used to compute the distance matrix on the
  right side. If there is no left-hand-side variable, the distances will
  be computed between all pairs of units. If `NULL`, all the variables
  in `data` will be used as covariates.

- data:

  a data frame containing the variables named in `formula`. If `formula`
  is `NULL`, all variables in `data` will be used as covariates.

- s.weights:

  when `var = NULL`, an optional vector of sampling weights used to
  compute the variances used in the Mahalanobis, scaled Euclidean, and
  robust Mahalanobis distances.

- var:

  for `mahalanobis_dist()`, a covariance matrix used to scale the
  covariates. For `scaled_euclidean_dist()`, either a covariance matrix
  (from which only the diagonal elements will be used) or a vector of
  variances used to scale the covariates. If `NULL`, these values will
  be calculated using formulas described in Details.

- discarded:

  a `logical` vector denoting which units are to be discarded or not.
  This is used only when `var = NULL`. The scaling factors will be
  computed only using the non-discarded units, but the distance matrix
  will be computed for all units (discarded and non-discarded).

- ...:

  ignored. Included to make cycling through these functions easier
  without having to change the arguments supplied.

## Value

A numeric distance matrix. When `formula` has a left-hand-side
(treatment) variable, the matrix will have one row for each treated unit
and one column for each control unit. Otherwise, the matrix will have
one row and one column for each unit.

## Details

The **Euclidean distance** (computed using `euclidean_dist()`) is the
raw distance between units, computed as \$\$d\_{ij} = \sqrt{(x_i -
x_j)(x_i - x_j)'}\$\$ where \\x_i\\ and \\x_j\\ are vectors of
covariates for units \\i\\ and \\j\\, respectively. The Euclidean
distance is sensitive to the scales of the variables and their
redundancy (i.e., correlation). It should probably not be used for
matching unless all of the variables have been previously scaled
appropriately or are already on the same scale. It forms the basis of
the other distance measures.

The **scaled Euclidean distance** (computed using
`scaled_euclidean_dist()`) is the Euclidean distance computed on the
scaled covariates. Typically the covariates are scaled by dividing by
their standard deviations, but any scaling factor can be supplied using
the `var` argument. This leads to a distance measure computed as
\$\$d\_{ij} = \sqrt{(x_i - x_j)S_d^{-1}(x_i - x_j)'}\$\$ where \\S_d\\
is a diagonal matrix with the squared scaling factors on the diagonal.
Although this measure is not sensitive to the scales of the variables
(because they are all placed on the same scale), it is still sensitive
to redundancy among the variables. For example, if 5 variables measure
approximately the same construct (i.e., are highly correlated) and 1
variable measures another construct, the first construct will have 5
times as much influence on the distance between units as the second
construct. The Mahalanobis distance attempts to address this issue.

The **Mahalanobis distance** (computed using `mahalanobis_dist()`) is
computed as \$\$d\_{ij} = \sqrt{(x_i - x_j)S^{-1}(x_i - x_j)'}\$\$ where
\\S\\ is a scaling matrix, typically the covariance matrix of the
covariates. It is essentially equivalent to the Euclidean distance
computed on the scaled principal components of the covariates. This is
the most popular distance matrix for matching because it is not
sensitive to the scale of the covariates and accounts for redundancy
between them. The scaling matrix can also be supplied using the `var`
argument.

The Mahalanobis distance can be sensitive to outliers and long-tailed or
otherwise non-normally distributed covariates and may not perform well
with categorical variables due to prioritizing rare categories over
common ones. One solution is the rank-based **robust Mahalanobis
distance** (computed using `robust_mahalanobis_dist()`), which is
computed by first replacing the covariates with their ranks (using
average ranks for ties) and rescaling each ranked covariate by a
constant scaling factor before computing the usual Mahalanobis distance
on the rescaled ranks.

The Mahalanobis distance and its robust variant are computed internally
by transforming the covariates in such a way that the Euclidean distance
computed on the scaled covariates is equal to the requested distance.
For the Mahalanobis distance, this involves replacing the covariates
vector \\x_i\\ with \\x_iS^{-.5}\\, where \\S^{-.5}\\ is the Cholesky
decomposition of the (generalized) inverse of the covariance matrix
\\S\\.

When a left-hand-side splitting variable is present in `formula` and
`var = NULL` (i.e., so that the scaling matrix is computed internally),
the covariance matrix used is the "pooled" covariance matrix, which
essentially is a weighted average of the covariance matrices computed
separately within each level of the splitting variable to capture
within-group variation and reduce sensitivity to covariate imbalance.
This is also true of the scaling factors used in the scaled Euclidean
distance.

## References

Rosenbaum, P. R. (2010). *Design of observational studies*. Springer.

Rosenbaum, P. R., & Rubin, D. B. (1985). Constructing a Control Group
Using Multivariate Matched Sampling Methods That Incorporate the
Propensity Score. *The American Statistician*, 39(1), 33–38.
[doi:10.2307/2683903](https://doi.org/10.2307/2683903)

Rubin, D. B. (1980). Bias Reduction Using Mahalanobis-Metric Matching.
*Biometrics*, 36(2), 293–298.
[doi:10.2307/2529981](https://doi.org/10.2307/2529981)

## See also

[`distance`](https://kosukeimai.github.io/MatchIt/reference/distance.md),
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md),
[`dist()`](https://rdrr.io/r/stats/dist.html) (which is used internally
to compute some Euclidean distances)

[`optmatch::match_on()`](https://rdrr.io/pkg/optmatch/man/match_on-methods.html),
which provides similar functionality but with fewer options and a focus
on efficient storage of the output.

## Author

Noah Greifer

## Examples

``` r

data("lalonde")

# Computing the scaled Euclidean distance between all units:
d <- scaled_euclidean_dist(~ age + educ + race + married,
                           data = lalonde)

# Another interface using the data argument:
dat <- subset(lalonde, select = c(age, educ, race, married))
d <- scaled_euclidean_dist(data = dat)

# Computing the Mahalanobis distance between treated and
# control units:
d <- mahalanobis_dist(treat ~ age + educ + race + married,
                      data = lalonde)

# Supplying a covariance matrix or vector of variances (note:
# a bit more complicated with factor variables)
dat <- subset(lalonde, select = c(age, educ, married, re74))
vars <- sapply(dat, var)

d <- scaled_euclidean_dist(data = dat, var = vars)

# Same result:
d <- scaled_euclidean_dist(data = dat, var = diag(vars))

# Discard units:
discard <- sample(c(TRUE, FALSE), nrow(lalonde),
                  replace = TRUE, prob = c(.2, .8))

d <- mahalanobis_dist(treat ~ age + educ + race + married,
                      data = lalonde, discarded = discard)
dim(d) #all units present in distance matrix
#> [1] 185 429
table(lalonde$treat)
#> 
#>   0   1 
#> 429 185 
```
