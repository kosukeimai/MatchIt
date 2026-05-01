# Construct a matched dataset from a `matchit` object

`match_data()` and `get_matches()` create a data frame with additional
variables for the distance measure, matching weights, and subclasses
after matching. This dataset can be used to estimate treatment effects
after matching or subclassification. `get_matches()` is most useful
after matching with replacement; otherwise, `match_data()` is more
flexible. See Details below for the difference between them.

## Usage

``` r
match_data(
  object,
  group = "all",
  distance = "distance",
  weights = "weights",
  subclass = "subclass",
  data = NULL,
  include.s.weights = TRUE,
  drop.unmatched = TRUE
)

match.data(...)

get_matches(
  object,
  distance = "distance",
  weights = "weights",
  subclass = "subclass",
  id = "id",
  data = NULL,
  include.s.weights = TRUE
)
```

## Arguments

- object:

  a `matchit` object; the output of a call to
  [`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md).

- group:

  which group should comprise the matched dataset: `"all"` for all
  units, `"treated"` for just treated units, or `"control"` for just
  control units. Default is `"all"`.

- distance:

  a string containing the name that should be given to the variable
  containing the distance measure in the data frame output. Default is
  `"distance"`, but `"prop.score"` or similar might be a good
  alternative if propensity scores were used in matching. Ignored if a
  distance measure was not supplied or estimated in the call to
  [`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md).

- weights:

  a string containing the name that should be given to the variable
  containing the matching weights in the data frame output. Default is
  `"weights"`.

- subclass:

  a string containing the name that should be given to the variable
  containing the subclasses or matched pair membership in the data frame
  output. Default is `"subclass"`.

- data:

  a data frame containing the original dataset to which the computed
  output variables (`distance`, `weights`, and/or `subclass`) should be
  appended. If empty, `match_data()` and `get_matches()` will attempt to
  find the dataset using the environment of the `matchit` object, which
  can be unreliable; see Notes.

- include.s.weights:

  `logical`; whether to multiply the estimated weights by the sampling
  weights supplied to
  [`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md),
  if any. Default is `TRUE`. If `FALSE`, the weights in the
  `match_data()` or `get_matches()` output should be multiplied by the
  sampling weights before being supplied to the function estimating the
  treatment effect in the matched data.

- drop.unmatched:

  `logical`; whether the returned data frame should contain all units
  (`FALSE`) or only units that were matched (i.e., have a matching
  weight greater than zero) (`TRUE`). Default is `TRUE` to drop
  unmatched units.

- ...:

  arguments passed to `match_data()`.

- id:

  a string containing the name that should be given to the variable
  containing the unit IDs in the data frame output. Default is `"id"`.
  Only used with `get_matches()`; for `match_data()`, the units IDs are
  stored in the row names of the returned data frame.

## Value

A data frame containing the data supplied in the `data` argument or in
the original call to
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
with the computed output variables appended as additional columns, named
according the arguments above. For `match_data()`, the `group` and
`drop.unmatched` arguments control whether only subsets of the data are
returned. See Details above for how `match_data()` and `get_matches()`
differ. Note that `get_matches` sorts the data by subclass and treatment
status, unlike `match_data()`, which uses the order of the data.

The returned data frame will contain the variables in the original data
set or dataset supplied to `data` and the following columns:

- distance:

  The propensity score, if estimated or supplied to the `distance`
  argument in
  [`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
  as a vector.

- weights:

  The computed matching weights. These must be used in effect estimation
  to correctly incorporate the matching.

- subclass:

  Matching strata membership. Units with the same value are in the same
  stratum.

- id:

  The ID of each unit, corresponding to the row names in the original
  data or dataset supplied to `data`. Only included in `get_matches`
  output. This column can be used to identify which rows belong to the
  same unit since the same unit may appear multiple times if reused in
  matching with replacement.

These columns will take on the name supplied to the corresponding
arguments in the call to `match_data()` or `get_matches()`. See Examples
for an example of rename the `distance` column to `"prop.score"`.

If `data` or the original dataset supplied to
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
was a `data.table` or `tbl`, the `match_data()` output will have the
same class, but the `get_matches()` output will always be a base R
`data.frame`.

In addition to their base class (e.g., `data.frame` or `tbl`), returned
objects have the class `matchdata` or `getmatches`. This class is
important when using
[`rbind()`](https://kosukeimai.github.io/MatchIt/reference/rbind.matchdata.md)
to append matched datasets.

## Details

`match_data()` creates a dataset with one row per unit. It will be
identical to the dataset supplied except that several new columns will
be added containing information related to the matching. When
`drop.unmatched = TRUE`, the default, units with weights of zero, which
are those units that were discarded by common support or the caliper or
were simply not matched, will be dropped from the dataset, leaving only
the subset of matched units. The idea is for the output of
`match_data()` to be used as the dataset input in calls to
[`glm()`](https://rdrr.io/r/stats/glm.html) or similar to estimate
treatment effects in the matched sample. It is important to include the
weights in the estimation of the effect and its standard error. The
subclass column, when created, contains pair or subclass membership and
should be used to estimate the effect and its standard error. Subclasses
will only be included if there is a `subclass` component in the
`matchit` object, which does not occur with matching with replacement,
in which case `get_matches()` should be used. See
[`vignette("estimating-effects")`](https://kosukeimai.github.io/MatchIt/articles/estimating-effects.md)
for information on how to use `match_data()` output to estimate effects.
`match.data()` is an alias for `match_data()`.

`get_matches()` is similar to `match_data()`; the primary difference
occurs when matching is performed with replacement, i.e., when units do
not belong to a single matched pair. In this case, the output of
`get_matches()` will be a dataset that contains one row per unit for
each pair they are a part of. For example, if matching was performed
with replacement and a control unit was matched to two treated units,
that control unit will have two rows in the output dataset, one for each
pair it is a part of. Weights are computed for each row, and, for
control units, are equal to the inverse of the number of control units
in each control unit's subclass; treated units get a weight of 1.
Unmatched units are dropped. An additional column with unit IDs will be
created (named using the `id` argument) to identify when the same unit
is present in multiple rows. This dataset structure allows for the
inclusion of both subclass membership and repeated use of units, unlike
the output of `match_data()`, which lacks subclass membership when
matching is done with replacement. A `match.matrix` component of the
`matchit` object must be present to use `get_matches()`; in some forms
of matching, it is absent, in which case `match_data()` should be used
instead. See
[`vignette("estimating-effects")`](https://kosukeimai.github.io/MatchIt/articles/estimating-effects.md)
for information on how to use `get_matches()` output to estimate effects
after matching with replacement.

## Note

The most common way to use `match_data()` and `get_matches()` is by
supplying just the `matchit` object, e.g., as `match_data(m.out)`. A
data set will first be searched in the environment of the `matchit`
formula, then in the calling environment of `match_data()` or
`get_matches()`, and finally in the `model` component of the `matchit`
object if a propensity score was estimated.

When called from an environment different from the one in which
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
was originally called and a propensity score was not estimated (or was
but with `discard` not `"none"` and `reestimate = TRUE`), this syntax
may not work because the original dataset used to construct the matched
dataset will not be found. This can occur when
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
was run within an [`lapply()`](https://rdrr.io/r/base/lapply.html) or
[`purrr::map()`](https://purrr.tidyverse.org/reference/map.html) call.
The solution, which is recommended in all cases, is simply to supply the
original dataset to the `data` argument of `match_data()`, e.g., as
`match_data(m.out, data = original_data)`, as demonstrated in the
Examples.

## See also

[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md);
[`rbind.matchdata()`](https://kosukeimai.github.io/MatchIt/reference/rbind.matchdata.md)

[`vignette("estimating-effects")`](https://kosukeimai.github.io/MatchIt/articles/estimating-effects.md)
for uses of `match_data()` and `get_matches()` in estimating treatment
effects.

## Examples

``` r

data("lalonde")

# 4:1 matching w/replacement
m.out1 <- matchit(treat ~ age + educ + married +
                    race + nodegree + re74 + re75,
                  data = lalonde,
                  replace = TRUE,
                  caliper = .05,
                  ratio = 4)

m.data1 <- match_data(m.out1,
                      data = lalonde,
                      distance = "prop.score")
dim(m.data1) #one row per matched unit
#> [1] 347  11
head(m.data1, 10)
#>       treat age educ   race married nodegree re74 re75       re78 prop.score
#> NSW1      1  37   11  black       1        1    0    0  9930.0460 0.63876993
#> NSW2      1  22    9 hispan       0        1    0    0  3595.8940 0.22463424
#> NSW3      1  30   12  black       0        0    0    0 24909.4500 0.67824388
#> NSW4      1  27   11  black       0        1    0    0  7506.1460 0.77632408
#> NSW5      1  33    8  black       0        1    0    0   289.7899 0.70163874
#> NSW6      1  22    9  black       0        1    0    0  4056.4940 0.69906990
#> NSW7      1  23   12  black       0        0    0    0     0.0000 0.65368426
#> NSW8      1  32   11  black       0        1    0    0  8472.1580 0.78972311
#> NSW9      1  22   16  black       0        0    0    0  2164.0220 0.77983825
#> NSW10     1  33   12  white       1        0    0    0 12418.0700 0.04292461
#>       weights
#> NSW1        1
#> NSW2        1
#> NSW3        1
#> NSW4        1
#> NSW5        1
#> NSW6        1
#> NSW7        1
#> NSW8        1
#> NSW9        1
#> NSW10       1

g.matches1 <- get_matches(m.out1,
                          data = lalonde,
                          distance = "prop.score")
dim(g.matches1) #multiple rows per matched unit
#> [1] 820  13
head(g.matches1, 10)
#>         id subclass weights treat age educ   race married nodegree        re74
#> 1     NSW1        1    1.00     1  37   11  black       1        1     0.00000
#> 2   PSID69        1    0.25     0  30   17  black       0        0 17827.37000
#> 3  PSID387        1    0.25     0  55    4  black       0        1     0.00000
#> 4  PSID373        1    0.25     0  20   12  black       0        0     0.00000
#> 5  PSID386        1    0.25     0  20   12  black       0        0     0.00000
#> 6     NSW2        2    1.00     1  22    9 hispan       0        1     0.00000
#> 7  PSID111        2    0.25     0  51   11  white       0        1    48.98167
#> 8   PSID66        2    0.25     0  26    8 hispan       0        1  3168.13400
#> 9  PSID339        2    0.25     0  26    9 hispan       0        1  1563.49500
#> 10 PSID150        2    0.25     0  22   11 hispan       0        1  7341.37300
#>        re75      re78 prop.score
#> 1     0.000  9930.046  0.6387699
#> 2  5546.419 14421.130  0.6385542
#> 3     0.000     0.000  0.6357539
#> 4     0.000     0.000  0.6428929
#> 5     0.000 11594.240  0.6428929
#> 6     0.000  3595.894  0.2246342
#> 7  3813.387  1525.014  0.2240794
#> 8  5872.258 11136.150  0.2225951
#> 9     0.000  2862.356  0.2161957
#> 10 2535.097 14187.650  0.2128728
```
