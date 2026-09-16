# Upcoming features

### Multi-category treatments

`cardinality_matchit()` and `cem_matchit()` already ready. Rest in R/MatchItMulti.


### Weights no longer scaled by default

Matching weights will come directly from formulas and not divide by total of weights. More transparent given questions often asked, and in most cases should keep group sizes equal (correctly). Could add option to scale weights to recover old behavior.

### eCDF Mean removed from `summary()`

No recommendations, hard to interpret, not sure if valid. Possible renaming of eCDF max to KS statistic.

### Matching w/ replacement for ATE

Simple to match for ATT and then to ATC. Need to adjust `match.matrix`. 

### Missingness indicator approach for missing data

Maybe to function similar to `WeightIt`. Possibly can use other missingness methods for different propensity scores, e.g., SAEM or surrogate.

### Hand-coded genetic matching

Exists in MatchItMulti

### Sampling weights in stratficiation

PS should be computed using weighted sums of units / weighted regression of treatment on subclass membership

### Different optimization algorithms for optimal matching

RcppHungarian already used in MatchItMulti; other used in couplr