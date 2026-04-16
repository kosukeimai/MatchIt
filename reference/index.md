# Package index

## Matching

- [`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
  : Matching for Causal Inference

- [`method_cardinality`](https://kosukeimai.github.io/MatchIt/reference/method_cardinality.md)
  : Cardinality Matching

- [`method_cem`](https://kosukeimai.github.io/MatchIt/reference/method_cem.md)
  : Coarsened Exact Matching

- [`method_exact`](https://kosukeimai.github.io/MatchIt/reference/method_exact.md)
  : Exact Matching

- [`method_full`](https://kosukeimai.github.io/MatchIt/reference/method_full.md)
  : Optimal Full Matching

- [`method_genetic`](https://kosukeimai.github.io/MatchIt/reference/method_genetic.md)
  : Genetic Matching

- [`method_nearest`](https://kosukeimai.github.io/MatchIt/reference/method_nearest.md)
  : Nearest Neighbor Matching

- [`method_optimal`](https://kosukeimai.github.io/MatchIt/reference/method_optimal.md)
  : Optimal Pair Matching

- [`method_quick`](https://kosukeimai.github.io/MatchIt/reference/method_quick.md)
  : Fast Generalized Full Matching

- [`method_subclass`](https://kosukeimai.github.io/MatchIt/reference/method_subclass.md)
  : Subclassification

- [`distance`](https://kosukeimai.github.io/MatchIt/reference/distance.md)
  : Propensity scores and other distance measures

- [`mahalanobis_dist()`](https://kosukeimai.github.io/MatchIt/reference/mahalanobis_dist.md)
  [`scaled_euclidean_dist()`](https://kosukeimai.github.io/MatchIt/reference/mahalanobis_dist.md)
  [`robust_mahalanobis_dist()`](https://kosukeimai.github.io/MatchIt/reference/mahalanobis_dist.md)
  [`euclidean_dist()`](https://kosukeimai.github.io/MatchIt/reference/mahalanobis_dist.md)
  : Compute a Distance Matrix

- [`add_s.weights()`](https://kosukeimai.github.io/MatchIt/reference/add_s.weights.md)
  :

  Add sampling weights to a `matchit` object

## Assessing Balance

- [`summary(`*`<matchit>`*`)`](https://kosukeimai.github.io/MatchIt/reference/summary.matchit.md)
  [`summary(`*`<matchit.subclass>`*`)`](https://kosukeimai.github.io/MatchIt/reference/summary.matchit.md)
  [`print(`*`<summary.matchit>`*`)`](https://kosukeimai.github.io/MatchIt/reference/summary.matchit.md)
  :

  View a balance summary of a `matchit` object

- [`plot(`*`<summary.matchit>`*`)`](https://kosukeimai.github.io/MatchIt/reference/plot.summary.matchit.md)
  : Generate a Love Plot of Standardized Mean Differences

- [`plot(`*`<matchit>`*`)`](https://kosukeimai.github.io/MatchIt/reference/plot.matchit.md)
  [`plot(`*`<matchit.subclass>`*`)`](https://kosukeimai.github.io/MatchIt/reference/plot.matchit.md)
  : Generate Balance Plots after Matching and Subclassification

## Extracting Matched Data

- [`match_data()`](https://kosukeimai.github.io/MatchIt/reference/match_data.md)
  [`match.data()`](https://kosukeimai.github.io/MatchIt/reference/match_data.md)
  [`get_matches()`](https://kosukeimai.github.io/MatchIt/reference/match_data.md)
  :

  Construct a matched dataset from a `matchit` object

- [`rbind(`*`<matchdata>`*`)`](https://kosukeimai.github.io/MatchIt/reference/rbind.matchdata.md)
  [`rbind(`*`<getmatches>`*`)`](https://kosukeimai.github.io/MatchIt/reference/rbind.matchdata.md)
  : Append matched datasets together

## Datasets

- [`lalonde`](https://kosukeimai.github.io/MatchIt/reference/lalonde.md)
  : Data from National Supported Work Demonstration and PSID, as
  analyzed by Dehejia and Wahba (1999).
