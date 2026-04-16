# Matching Methods

## Introduction

`MatchIt` implements several matching methods with a variety of options.
Though the help pages for the individual methods describe each method
and how they can be used, this vignette provides a broad overview of the
available matching methods and their associated options. The choice of
matching method depends on the goals of the analysis (e.g., the
estimand, whether low bias or high precision is important) and the
unique qualities of each dataset to be analyzed, so there is no single
optimal choice for any given analysis. A benefit of nonparametric
preprocessing through matching is that a number of matching methods can
be tried and their quality assessed without consulting the outcome,
reducing the possibility of capitalizing on chance while allowing for
the benefits of an exploratory analysis in the design phase ([Ho et al.
2007](#ref-ho2007)).

This vignette describes each matching method available in `MatchIt` and
the various options that are allowed with matching methods and the
consequences of their use. For a brief introduction to the use of
`MatchIt` functions, see
[`vignette("MatchIt")`](https://kosukeimai.github.io/MatchIt/articles/MatchIt.md).
For details on how to assess and report covariate balance, see
[`vignette("assessing-balance")`](https://kosukeimai.github.io/MatchIt/articles/assessing-balance.md).
For details on how to estimate treatment effects and standard errors
after matching, see
[`vignette("estimating-effects")`](https://kosukeimai.github.io/MatchIt/articles/estimating-effects.md).

## Matching

Matching as implemented in `MatchIt` is a form of *subset selection*,
that is, the pruning and weighting of units to arrive at a (weighted)
subset of the units from the original dataset. Ideally, and if done
successfully, subset selection produces a new sample where the treatment
is unassociated with the covariates so that a comparison of the outcomes
treatment and control groups is not confounded by the measured and
balanced covariates. Although statistical estimation methods like
regression can also be used to remove confounding due to measured
covariates, Ho et al. ([2007](#ref-ho2007)) argue that fitting
regression models in matched samples reduces the dependence of the
validity of the estimated treatment effect on the correct specification
of the model.

Matching is nonparametric in the sense that the estimated weights and
pruning of the sample are not direct functions of estimated model
parameters but rather depend on the organization of discrete units in
the sample; this is in contrast to propensity score weighting (also
known as inverse probability weighting), where the weights come more
directly from the estimated propensity score model and therefore are
more sensitive to its correct specification. These advantages, as well
as the intuitive understanding of matching by the public compared to
regression or weighting, make it a robust and effective way to estimate
treatment effects.

It is important to note that this implementation of matching differs
from the methods described by Abadie and Imbens
([2006](#ref-abadie2006), [2016](#ref-abadie2016)) and implemented in
the `Matching` R package and `teffects` routine in Stata. That form of
matching is *matching imputation*, where the missing potential outcomes
for each unit are imputed using the observed outcomes of paired units.
This is a critical distinction because matching imputation is a specific
estimation method with its own effect and standard error estimators, in
contrast to subset selection, which is a preprocessing method that does
not require specific estimators and is broadly compatible with other
parametric and nonparametric analyses. The benefits of matching
imputation are that its theoretical properties (i.e., the rate of
convergence and asymptotic variance of the estimator) are well
understood, it can be used in a straightforward way to estimate not just
the average treatment effect in the treated (ATT) but also the average
treatment effect in the population (ATE), and additional effective
matching methods can be used in the imputation (e.g., kernel matching).
The benefits of matching as nonparametric preprocessing are that it is
far more flexible with respect to the types of effects that can be
estimated because it does not involve any specific estimator, its
empirical and finite-sample performance has been examined in depth and
is generally well understood, and it aligns well with the design of
experiments, which are more familiar to non-technical audiences.

In addition to subset selection, matching often (though not always)
involves a form of *stratification*, the assignment of units to pairs or
strata containing multiple units. The distinction between subset
selection and stratification is described by Zubizarreta, Paredes, and
Rosenbaum ([2014](#ref-zubizarretaMatchingBalancePairing2014)), who
separate them into two separate steps. In `MatchIt`, with almost all
matching methods, subset selection is performed by stratification; for
example, treated units are paired with control units, and unpaired units
are then dropped from the matched sample. With some methods, subclasses
are used to assign matching or stratification weights to individual
units, which increase or decrease each unit’s leverage in a subsequent
analysis. There has been some debate about the importance of
stratification after subset selection; while some authors have argued
that, with some forms of matching, pair membership is incidental
([Stuart 2008](#ref-stuart2008); [Schafer and Kang
2008](#ref-schafer2008)), others have argued that correctly
incorporating pair membership into effect estimation can improve the
quality of inferences ([Austin and Small 2014](#ref-austin2014a); [Wan
2019](#ref-wan2019)). For methods that allow it, `MatchIt` includes
stratum membership as an additional output of each matching
specification. How these strata can be used is detailed in
[`vignette("estimating-effects")`](https://kosukeimai.github.io/MatchIt/articles/estimating-effects.md).

At the heart of `MatchIt` are three classes of methods: distance
matching, stratum matching, and pure subset selection. *Distance
matching* involves considering a focal group (usually the treated group)
and selecting members of the non-focal group (i.e., the control group)
to pair with each member of the focal group based on the *distance*
between units, which can be computed in one of several ways. Members of
either group that are not paired are dropped from the sample. Nearest
neighbor matching (`method = "nearest"`), optimal pair matching
(`method = "optimal"`), optimal full matching (`method = "full"`),
generalized full matching (`method = "quick"`), and genetic matching
(`method = "genetic"`) are the methods of distance matching implemented
in `MatchIt`. Typically, only the average treatment in the treated (ATT)
or average treatment in the control (ATC), if the control group is the
focal group, can be estimated after distance matching in `MatchIt` (full
matching is an exception, described later).

*Stratum matching* involves creating strata based on unique values of
the covariates and assigning units with those covariate values into
those strata. Any units that are in strata that lack either treated or
control units are then dropped from the sample. Strata can be formed
using the raw covariates (`method = "exact"`), coarsened versions of the
covariates (`method = "cem"`), or coarsened versions of the propensity
score (`method = "subclass"`). When no units are discarded, either the
ATT, ATC, or ATE can be estimated after stratum matching, though often
some units are discarded, especially with exact and coarsened exact
matching, making the estimand less clear. For use in estimating marginal
treatment effects after exact matching, stratification weights are
computed for the matched units first by computing a new “stratum
propensity score” for each unit, which is the proportion of treated
units in its stratum. The formulas for computing inverse probability
weights from standard propensity scores are then applied to the new
stratum propensity scores to form the new weights.

Pure subset selection involves selecting a subset of units form the
original sample without considering the distance between individual
units or strata that units might fall into. Subsets are selected to
optimize a criterion subject to constraint on balance and remaining
sample size. Cardinality and profile matching (`method = "cardinality"`)
are the methods of pure subset selection implemented in `MatchIt`. Both
methods allow the user to specify the largest imbalance allowed in the
resulting matched sample, and an optimization routine attempts to find
the largest matched sample that satisfies those balance constraints.
While cardinality matching does not target a specific estimand, profile
matching can be used to target the ATT, ATC, or ATE.

Below, we describe each of the matching methods implemented in
`MatchIt`.

## Matching Methods

### Nearest Neighbor Matching (`method = "nearest"`)

Nearest neighbor matching is also known as greedy matching. It involves
running through the list of treated units and selecting the closest
eligible control unit to be paired with each treated unit. It is greedy
in the sense that each pairing occurs without reference to how other
units will be or have been paired, and therefore does not aim to
optimize any criterion. Nearest neighbor matching is the most common
form of matching used ([Thoemmes and Kim 2011](#ref-thoemmes2011);
[Zakrison, Austin, and McCredie 2018](#ref-zakrison2018)) and has been
extensively studied through simulations. See
[`?method_nearest`](https://kosukeimai.github.io/MatchIt/reference/method_nearest.md)
for the documentation for
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
with `method = "nearest"`.

Nearest neighbor matching requires the specification of a distance
measure to define which control unit is closest to each treated unit.
The default and most common distance is the *propensity score
difference*, which is the difference between the propensity scores of
each treated and control unit ([Stuart 2010](#ref-stuart2010)). Another
popular distance is the Mahalanobis distance, described in the section
“Mahalanobis distance matching” below. The order in which the treated
units are to be paired must also be specified and has the potential to
change the quality of the matches ([Austin 2013](#ref-austin2013b);
[Rubin 1973](#ref-rubin1973)); this is specified by the `m.order`
argument. With propensity score matching, the default is to go in
descending order from the highest propensity score; doing so allows the
units that would have the hardest time finding close matches to be
matched first ([Rubin 1973](#ref-rubin1973)). Other orderings are
possible, including random ordering, which can be tried multiple times
until an adequate matched sample is found. When matching with
replacement (i.e., where each control unit can be reused to be matched
with any number of treated units), the matching order doesn’t matter.

When using a matching ratio greater than 1 (i.e., when more than 1
control units are requested to be matched to each treated unit),
matching occurs in a cycle, where each treated unit is first paired with
one control unit, and then each treated unit is paired with a second
control unit, etc. Ties are broken deterministically based on the order
of the units in the dataset to ensure that multiple runs of the same
specification yield the same result (unless the matching order is
requested to be random).

Nearest neighbor matching is implemented in `MatchIt` using internal C++
code through `Rcpp`. When matching on a propensity score, this makes
matching extremely fast, even for large datasets. Using a caliper on the
propensity score (described below) makes it even faster. Run times may
be a bit longer when matching on other distance measures (e.g., the
Mahalanobis distance). In contrast to optimal pair matching (described
below), nearest neighbor matching does not require computing the full
distance matrix between units, which makes it more applicable to large
datasets.

### Optimal Pair Matching (`method = "optimal"`)

Optimal pair matching (often just called optimal matching) is very
similar to nearest neighbor matching in that it attempts to pair each
treated unit with one or more control units. Unlike nearest neighbor
matching, however, it is “optimal” rather than greedy; it is optimal in
the sense that it attempts to choose matches that collectively optimize
an overall criterion ([Hansen and Klopfer 2006](#ref-hansen2006); [Gu
and Rosenbaum 1993](#ref-gu1993)). The criterion used is the sum of the
absolute pair distances in the matched sample. See
[`?method_optimal`](https://kosukeimai.github.io/MatchIt/reference/method_optimal.md)
for the documentation for
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
with `method = "optimal"`. Optimal pair matching in `MatchIt` depends on
the `fullmatch()` function in the `optmatch` package ([Hansen and
Klopfer 2006](#ref-hansen2006)).

Like nearest neighbor matching, optimal pair matching requires the
specification of a distance measure between units. Optimal pair matching
can be thought of simply as an alternative to selecting the order of the
matching for nearest neighbor matching. Optimal pair matching and
nearest neighbor matching often yield the same or very similar matched
samples; indeed, some research has indicated that optimal pair matching
is not much better than nearest neighbor matching at yielding balanced
matched samples ([Austin 2013](#ref-austin2013b)).

The `tol` argument in `fullmatch()` can be supplied to
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
with `method = "optimal"`; this controls the numerical tolerance used to
determine whether the optimal solution has been found. The default is
fairly high and, for smaller problems, should be set much lower (e.g.,
by setting `tol = 1e-7`).

### Optimal Full Matching (`method = "full"`)

Optimal full matching (often just called full matching) assigns every
treated and control unit in the sample to one subclass each ([Hansen
2004](#ref-hansen2004); [Stuart and Green 2008](#ref-stuart2008a)). Each
subclass contains one treated unit and one or more control units or one
control units and one or more treated units. It is optimal in the sense
that the chosen number of subclasses and the assignment of units to
subclasses minimize the sum of the absolute within-subclass distances in
the matched sample. Weights are computed based on subclass membership,
and these weights then function like propensity score weights and can be
used to estimate a weighted treatment effect, ideally free of
confounding by the measured covariates. See
[`?method_full`](https://kosukeimai.github.io/MatchIt/reference/method_full.md)
for the documentation for
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
with `method = "full"`. Optimal full matching in `MatchIt` depends on
the `fullmatch()` function in the `optmatch` package ([Hansen and
Klopfer 2006](#ref-hansen2006)).

Like the other distance matching methods, optimal full matching requires
the specification of a distance measure between units. It can be seen a
combination of distance matching and stratum matching: subclasses are
formed with varying numbers of treated and control units, as with
stratum matching, but the subclasses are formed based on minimizing
within-pair distances and do not involve forming strata based on any
specific variable, similar to distance matching. Unlike other distance
matching methods, full matching can be used to estimate the ATE. Full
matching can also be seen as a form of propensity score weighting that
is less sensitive to the form of the propensity score model because the
original propensity scores are used just to create the subclasses, not
to form the weights directly ([Austin and Stuart
2015a](#ref-austin2015a)). In addition, full matching does not have to
rely on estimated propensity scores to form the subclasses and weights;
other distance measures are allowed as well.

Although full matching uses all available units, there is a loss in
precision due to the weights. Units may be weighted in such a way that
they contribute less to the sample than would unweighted units, so the
effective sample size (ESS) of the full matching weighted sample may be
lower than even that of 1:1 pair matching. Balance is often far better
after full matching than it is with 1:k matching, making full matching a
good option to consider especially when 1:k matching is not effective or
when the ATE is the target estimand.

The specification of the full matching optimization problem can be
customized by supplying additional arguments that are passed to
[`optmatch::fullmatch()`](https://rdrr.io/pkg/optmatch/man/fullmatch.html),
such as `min.controls`, `max.controls`, `mean.controls`, and
`omit.fraction`. As with optimal pair matching, the numerical tolerance
value can be set much lower than the default with small problems by
setting, e.g., `tol = 1e-7`.

### Generalized Full Matching (`method = "quick"`)

Generalized full matching is a variant of full matching that uses a
special fast clustering algorithm to dramatically speed up the matching,
even for large datasets ([Fredrik Sävje, Higgins, and Sekhon
2021](#ref-savjeGeneralizedFullMatching2021)). Like with optimal full
matching, generalized full matching assigns every unit to a subclass.
What makes generalized full match “generalized” is that the user can
customize the matching in a number of ways, such as by specifying an
arbitrary minimum number of units from each treatment group or total
number of units per subclass, or by allowing not all units from a
treatment group to have to be matched. Generalized full matching
minimizes the largest within-subclass distances in the matched sample,
but it does so in a way that is not completely optimal (though the
solution is often very close to the optimal solution). Matching weights
are computed based on subclass membership, and these weights then
function like propensity score weights and can be used to estimate a
weighted treatment effect, ideally free of confounding by the measured
covariates. See
[`?method_quick`](https://kosukeimai.github.io/MatchIt/reference/method_quick.md)
for the documentation for
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
with `method = "quick"`. Generalized full matching in `MatchIt` depends
on the `quickmatch()` function in the `quickmatch` package ([Fredrik
Sävje, Sekhon, and Higgins
2018](#ref-savjeQuickmatchQuickGeneralized2018)).

Generalized full matching includes different options for customization
than optimal full matching. The user cannot supply their own distance
matrix, but propensity scores and distance metrics that are computed
from the supplied covariates (e.g., Mahalanobis distance) are allowed.
Calipers can only be placed on the propensity score, if supplied. As
with optimal full matching, generalized full matching can target the
ATE. Matching performance tends to be similar between the two methods,
but generalized full matching will be much quicker and can accommodate
larger datasets, making it a good substitute. Generalized full matching
is often faster than even nearest neighbor matching, especially for
large datasets.

### Genetic Matching (`method = "genetic"`)

Genetic matching is less a specific form of matching and more a way of
specifying a distance measure for another form of matching. In practice,
though, the form of matching used is nearest neighbor pair matching.
Genetic matching uses a genetic algorithm, which is an optimization
routine used for non-differentiable objective functions, to find scaling
factors for each variable in a generalized Mahalanobis distance formula
([Diamond and Sekhon 2013](#ref-diamond2013)). The criterion optimized
by the algorithm is one based on covariate balance. Once the scaling
factors have been found, nearest neighbor matching is performed on the
scaled generalized Mahalanobis distance. See
[`?method_genetic`](https://kosukeimai.github.io/MatchIt/reference/method_genetic.md)
for the documentation for
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
with `method = "genetic"`. Genetic matching in `MatchIt` depends on the
`GenMatch()` function in the `Matching` package ([Sekhon
2011](#ref-sekhon2011)) to perform the genetic search and uses the
`Match()` function to perform the nearest neighbor match using the
scaled generalized Mahalanobis distance.

Genetic matching considers the generalized Mahalanobis distance between
a treated unit $i$ and a control unit $j$ as
$$\delta_{GMD}\left( \mathbf{x}_{i},\mathbf{x}_{j},\mathbf{W} \right) = \sqrt{\left( \mathbf{x}_{i} - \mathbf{x}_{j} \right)\prime\left( \mathbf{S}^{- 1/2} \right)\prime\mathbf{W}\left( \mathbf{S}^{- 1/2} \right)\left( \mathbf{x}_{i} - \mathbf{x}_{j} \right)}$$
where $\mathbf{x}$ is a $p \times 1$ vector containing the value of each
of the $p$ included covariates for that unit, $\mathbf{S}^{- 1/2}$ is
the Cholesky decomposition of the covariance matrix $\mathbf{S}$ of the
covariates, and $\mathbf{W}$ is a diagonal matrix with scaling factors
$w$ on the diagonal: $$\mathbf{W} = \begin{bmatrix}
w_{1} & & & \\
 & w_{2} & & \\
 & & \ddots & \\
 & & & w_{p} \\
 & & & 
\end{bmatrix}$$

When $w_{k} = 1$ for all covariates $k$, the computed distance is the
standard Mahalanobis distance between units. Genetic matching estimates
the optimal values of the $w_{k}$s, where a user-specified criterion is
used to define what is optimal. The default is to maximize the smallest
p-value among balance tests for the covariates in the matched sample
(both Kolmogorov-Smirnov tests and t-tests for each covariate).

In `MatchIt`, if a propensity score is specified, the default is to
include the propensity score and the covariates in $\mathbf{x}$ and to
optimize balance on the covariates. When `distance = "mahalanobis"` or
the `mahvars` argument is specified, the propensity score is left out of
$\mathbf{x}$.

In all other respects, genetic matching functions just like nearest
neighbor matching except that the matching itself is carried out by
[`Matching::Match()`](https://rdrr.io/pkg/Matching/man/Match.html)
instead of by `MatchIt`. When using `method = "genetic"` in `MatchIt`,
additional arguments passed to
[`Matching::GenMatch()`](https://rdrr.io/pkg/Matching/man/GenMatch.html)
to control the genetic search process should be specified; in
particular, the `pop.size` argument should be increased from its default
of 100 to a much higher value. Doing so will make the algorithm take
more time to finish but will generally improve the quality of the
resulting matches. Different functions can be supplied to be used as the
objective in the optimization using the `fit.func` argument.

### Exact Matching (`method = "exact"`)

Exact matching is a form of stratum matching that involves creating
subclasses based on unique combinations of covariate values and
assigning each unit into their corresponding subclass so that only units
with identical covariate values are placed into the same subclass. Any
units that are in subclasses lacking either treated or control units
will be dropped. Exact matching is the most powerful matching method in
that no functional form assumptions are required on either the treatment
or outcome model for the method to remove confounding due to the
measured covariates; the covariate distributions are exactly balanced.
The problem with exact matching is that in general, few if any units
will remain after matching, so the estimated effect will only generalize
to a very limited population and can lack precision. Exact matching is
particularly ineffective with continuous covariates, for which it might
be that no two units have the same value, and with many covariates, for
which it might be the case that no two units have the same combination
of all covariates; this latter problem is known as the “curse of
dimensionality”. See
[`?method_exact`](https://kosukeimai.github.io/MatchIt/reference/method_exact.md)
for the documentation for
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
with `method = "exact"`.

It is possible to use exact matching on some covariates and another form
of matching on the rest. This makes it possible to have exact balance on
some covariates (typically categorical) and approximate balance on
others, thereby gaining the benefits of both exact matching and the
other matching method used. To do so, the other matching method should
be specified in the `method` argument to
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
and the `exact` argument should be specified to contain the variables on
which exact matching is to be done.

### Coarsened Exact Matching (`method = "cem"`)

Coarsened exact matching (CEM) is a form of stratum matching that
involves first coarsening the covariates by creating bins and then
performing exact matching on the new coarsened versions of the
covariates ([Iacus, King, and Porro 2012](#ref-iacus2012)). The degree
and method of coarsening can be controlled by the user to manage the
trade-off between exact and approximate balancing. For example,
coarsening a covariate to two bins will mean that units that differ
greatly on the covariate might be placed into the same subclass, while
coarsening a variable to five bins may require units to be dropped due
to not finding matches. Like exact matching, CEM is susceptible to the
curse of dimensionality, making it a less viable solution with many
covariates, especially with few units. Dropping units can also change
the target population of the estimated effect. See
[`?method_cem`](https://kosukeimai.github.io/MatchIt/reference/method_cem.md)
for the documentation for
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
with `method = "cem"`. CEM in `MatchIt` does not depend on any other
package to perform the coarsening and matching, though it used to rely
on the `cem` package.

### Subclassification (`method = "subclass"`)

Propensity score subclassification can be thought of as a form of
coarsened exact matching with the propensity score as the sole covariate
to be coarsened and matched on. The bins are usually based on specified
quantiles of the propensity score distribution either in the treated
group, control group, or overall, depending on the desired estimand.
Propensity score subclassification is an old and well-studied method,
though it can perform poorly compared to other, more modern propensity
score methods such as full matching and weighting ([Austin
2010a](#ref-austin2010)). See
[`?method_subclass`](https://kosukeimai.github.io/MatchIt/reference/method_subclass.md)
for the documentation for
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
with `method = "subclass"`.

The binning of the propensity scores is typically based on dividing the
distribution of covariates into approximately equally sized bins. The
user specifies the number of subclasses using the `subclass` argument
and which group should be used to compute the boundaries of the bins
using the `estimand` argument. Sometimes, subclasses can end up with no
units from one of the treatment groups; by default,
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
moves a unit from an adjacent subclass into the lacking one to ensure
that each subclass has at least one unit from each treatment group. The
minimum number of units required in each subclass can be chosen by the
`min.n` argument to
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md).
If set to 0, an error will be thrown if any subclass lacks units from
one of the treatment groups. Moving units from one subclass to another
generally worsens the balance in the subclasses but can increase
precision.

The default number of subclasses is 6, which is arbitrary and should not
be taken as a recommended value. Although early theory has recommended
the use of 5 subclasses, in general there is an optimal number of
subclasses that is typically much larger than 5 but that varies among
datasets ([Orihara and Hamada 2021](#ref-orihara2021)). Rather than
trying to figure this out for oneself, one can use optimal full matching
(i.e., with `method = "full"`) or generalized full matching
(`method = "quick"`) to optimally create subclasses that optimize a
within-subclass distance criterion.

The output of propensity score subclassification includes the assigned
subclasses and the subclassification weights. Effects can be estimated
either within each subclass and then averaged across them, or a single
marginal effect can be estimated using the subclassification weights.
This latter method has been called marginal mean weighting through
subclassification \[MMWS; Hong ([2010](#ref-hong2010))\] and fine
stratification weighting ([Desai et al. 2017](#ref-desai2017)). It is
also implemented in the `WeightIt` package.

### Cardinality and Profile Matching (`method = "cardinality"`)

Cardinality and profile matching are pure subset selection methods that
involve selecting a subset of the original sample without considering
the distance between individual units or assigning units to pairs or
subclasses. They can be thought of as a weighting method where the
weights are restricted to be zero or one. Cardinality matching involves
finding the largest sample that satisfies user-supplied balance
constraints and constraints on the ratio of matched treated to matched
control units ([Zubizarreta, Paredes, and Rosenbaum
2014](#ref-zubizarretaMatchingBalancePairing2014)). It does not consider
a specific estimand and can be a useful alternative to matching with a
caliper for handling data with little overlap ([Visconti and Zubizarreta
2018](#ref-visconti2018)). Profile matching involves identifying a
target distribution (e.g., the full sample for the ATE or the treated
units for the ATT) and finding the largest subset of the treated and
control groups that satisfy user-supplied balance constraints with
respect to that target ([Cohn and Zubizarreta
2022](#ref-cohnProfileMatchingGeneralization2021)). See
[`?method_cardinality`](https://kosukeimai.github.io/MatchIt/reference/method_cardinality.md)
for the documentation for using
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
with `method = "cardinality"`, including which inputs are required to
request either cardinality matching or profile matching.

Subset selection is performed by solving a mixed integer programming
optimization problem with linear constraints. The problem involves
maximizing the size of the matched sample subject to constraints on
balance and sample size. For cardinality matching, the balance
constraints refer to the mean difference for each covariate between the
matched treated and control groups, and the sample size constraints
require the matched treated and control groups to be the same size (or
differ by a user-supplied factor). For profile matching, the balance
constraints refer to the mean difference for each covariate between each
treatment group and the target distribution; for the ATE, this requires
the mean of each covariate in each treatment group to be within a given
tolerance of the mean of the covariate in the full sample, and for the
ATT, this requires the mean of each covariate in the control group to be
within a given tolerance of the mean of the covariate in the treated
group, which is left intact. The balance tolerances are controlled by
the `tols` and `std.tols` arguments. One can also create pairs in the
matched sample by using the `mahvars` argument, which requests that
optimal Mahalanobis matching be done after subset selection; doing so
can add additional precision and robustness ([Zubizarreta, Paredes, and
Rosenbaum 2014](#ref-zubizarretaMatchingBalancePairing2014)).

The optimization problem requires a special solver to solve. Currently,
the available options in `MatchIt` are the HiGHS solver (through the
`highs` package), the GLPK solver (through the `Rglpk` package), the
SYMPHONY solver (through the `Rsymphony` package), and the Gurobi solver
(through the `gurobi` package). The differences among the solvers are in
performance; Gurobi is by far the best (fastest, least likely to fail to
find a solution), but it is proprietary (though has a free trial and
academic license) and is a bit more complicated to install. HiGHS is the
default due to being open source, easily installed, and with performance
comparable to Gurobi. The `designmatch` package also provides an
implementation of cardinality matching with more options than `MatchIt`
offers.

## Customizing the Matching Specification

In addition to the specific matching method, other options are available
for many of the matching methods to further customize the matching
specification. These include different specifications of the distance
measure, methods to perform alternate forms of matching in addition to
the main method, prune units far from other units prior to matching,
restrict possible matches, etc. Not all options are compatible with all
matching methods.

### Specifying the propensity score or other distance measure (`distance`)

The distance measure is used to define how close two units are. In
nearest neighbor matching, this is used to choose the nearest control
unit to each treated unit. In optimal matching, this is used in the
criterion that is optimized. By default, the distance measure is the
propensity score difference, and the argument supplied to `distance`
corresponds to the method of estimating the propensity score. In
`MatchIt`, propensity scores are often labeled as “distance” values,
even though the propensity score itself is not a distance measure. This
is to reflect that the propensity score is used in creating the distance
value, but other scores could be used, such as prognostic scores for
prognostic score matching ([Hansen 2008](#ref-hansen2008a)). The
propensity score is more like a “position” value, in that it reflects
the position of each unit in the matching space, and the difference
between positions is the distance between them. If the argument to
`distance` is one of the allowed methods for estimating propensity
scores (see
[`?distance`](https://kosukeimai.github.io/MatchIt/reference/distance.md)
for these values) or is a numeric vector with one value per unit, the
distance between units will be computed as the pairwise difference
between propensity scores or the supplied values. Propensity scores are
also used in propensity score subclassification and can optionally be
used in genetic matching as a component of the generalized Mahalanobis
distance. For exact, coarsened exact, and cardinality matching, the
`distance` argument is ignored.

The default `distance` argument is `"glm"`, which estimates propensity
scores using logistic regression or another generalized linear model.
The `link` and `distance.options` arguments can be supplied to further
specify the options for the propensity score models, including whether
to use the raw propensity score or a linearized version of it (e.g., the
logit of a logistic regression propensity score, which has been commonly
referred to and recommended in the propensity score literature ([Austin
2011](#ref-austin2011a); [Stuart 2010](#ref-stuart2010))). Allowable
options for the propensity score model include parametric and machine
learning-based models, each of which have their strengths and
limitations and may perform differently depending on the unique
qualities of each dataset. We recommend multiple types of models be
tried to find one that yields the best balance, as there is no way to
make a single recommendation that will work for all cases.

The `distance` argument can also be specified as a method of computing
pairwise distances from the covariates directly (i.e., without
estimating propensity scores). The options include `"mahalanobis"`,
`"robust_mahalanobis"`, `"euclidean"`, and `"scaled_euclidean"`. These
methods compute a distance metric for a treated unit $i$ and a control
unit $j$ as
$$\delta\left( \mathbf{x}_{i},\mathbf{x}_{j} \right) = \sqrt{\left( \mathbf{x}_{i} - \mathbf{x}_{j} \right)\prime S^{- 1}\left( \mathbf{x}_{i} - \mathbf{x}_{j} \right)}$$

where $\mathbf{x}$ is a $p \times 1$ vector containing the value of each
of the $p$ included covariates for that unit, $S$ is a scaling matrix,
and $S^{- 1}$ is the (generalized) inverse of $S$. For Mahalanobis
distance matching, $S$ is the pooled covariance matrix of the covariates
([Rubin 1980](#ref-rubinBiasReductionUsing1980)); for Euclidean distance
matching, $S$ is the identity matrix (i.e., no scaling); and for scaled
Euclidean distance matching, $S$ is the diagonal of the pooled
covariance matrix (containing just the variances). The robust
Mahalanobis distance is computed not on the covariates directly but
rather on their ranks and uses a correction for ties (see Rosenbaum
([2010](#ref-rosenbaumDesignObservationalStudies2010)), ch 8). For
creating close pairs, matching with these distance measures tends work
better than propensity score matching because paired units will have
close values on all of the covariates, whereas propensity score-paired
units may be close on the propensity score but not on any of the
covariates themselves. This feature was the basis of King and Nielsen’s
([2019](#ref-king2019)) warning against using propensity scores for
matching. That said, they do not always outperform propensity score
matching ([Ripollone et al.
2018](#ref-ripolloneImplicationsPropensityScore2018)).

`distance` can also be supplied as a matrix of distance values between
units. This makes it possible to use handcrafted distance matrices or
distances created outside `MatchIt`. Only nearest neighbor, optimal
pair, and optimal full matching allow this specification.

The propensity score can have uses other than as the basis for matching.
It can be used to define a region of common support, outside which units
are dropped prior to matching; this is implemented by the `discard`
option. It can also be used to define a caliper, the maximum distance
two units can be before they are prohibited from being paired with each
other; this is implemented by the `caliper` argument. To estimate or
supply a propensity score for one of these purposes but not use it as
the distance measure for matching (i.e., to perform Mahalanobis distance
matching instead), the `mahvars` argument can be specified. These
options are described below.

### Implementing common support restrictions (`discard`)

The region of *common support* is the region of overlap between
treatment groups. A common support restriction discards units that fall
outside of the region of common support, preventing them from being
matched to other units and included in the matched sample. This can
reduce the potential for extrapolation and help the matching algorithms
to avoid overly distant matches from occurring. In `MatchIt`, the
`discard` option implements a common support restriction based on the
propensity score. The argument can be supplied as `"treated"`,
`"control"`, or `"both"`, which discards units in the corresponding
group that fall outside the region of common support for the propensity
score. The `reestimate` argument can be supplied to choose whether to
re-estimate the propensity score in the remaining units. **If units from
the treated group are discarded based on a common support restriction,
the estimand no longer corresponds to the ATT.**

### Caliper matching (`caliper`)

A *caliper* can be though of as a ring around each unit that limits to
which other units that unit can be paired. Calipers are based on the
propensity score or other covariates. Two units whose distance on a
calipered covariate is larger than the caliper width for that covariate
are not allowed to be matched to each other. Any units for which there
are no available matches within the caliper are dropped from the matched
sample. Calipers ensure paired units are close to each other on the
calipered covariates, which can ensure good balance in the matched
sample. Multiple variables can be supplied to `caliper` to enforce
calipers on all of them simultaneously. Using calipers can be a good
alternative to exact or coarsened exact matching to ensure only similar
units are paired with each other. The `std.caliper` argument controls
whether the provided calipers are in raw units or standard deviation
units. When negative calipers are supplied, this forces units whose
distance on the calipered covariate is *smaller* than the absolute
caliper width for that covariate to be disallowed from being matched to
each other. **If units from the treated group are left unmatched due to
a caliper, the estimand no longer corresponds to the ATT.**

### Mahalanobis distance matching (`mahvars`)

To perform Mahalanobis distance matching without the need to estimate or
use a propensity score, the `distance` argument can be set to
`"mahalanobis"`. If a propensity score is to be estimated or used for a
different purpose, such as in a common support restriction or a caliper,
but you still want to perform Mahalanobis distance matching, variables
should be supplied to the `mahvars` argument. The propensity scores will
be generated using the `distance` specification, and matching will occur
not on the covariates supplied to the main formula of
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.md)
but rather on the covariates supplied to `mahvars`. To perform
Mahalanobis distance matching within a propensity score caliper, for
example, the `distance` argument should be set to the method of
estimating the propensity score (e.g., `"glm"` for logistic regression),
the `caliper` argument should be specified to the desired caliper width,
and `mahvars` should be specified to perform Mahalanobis distance
matching on the desired covariates within the caliper. `mahvars` has a
special meaning for genetic matching and cardinality matching; see their
respective help pages for details.

### Exact matching (`exact`)

To perform exact matching on all supplied covariates, the `method`
argument can be set to `"exact"`. To perform exact matching only on some
covariates and some other form of matching within exact matching strata
on other covariates, the `exact` argument can be used. Covariates
supplied to the `exact` argument will be matched exactly, and the form
of matching specified by `method` (e.g., `"nearest"` for nearest
neighbor matching) will take place within each exact matching stratum.
This can be a good way to gain some of the benefits of exact matching
without completely succumbing to the curse of dimensionality. As with
exact matching performed with `method = "exact"`, any units in strata
lacking members of one of the treatment groups will be left unmatched.
Note that although matching occurs within each exact matching stratum,
propensity score estimation and computation of the Mahalanobis or other
distance matrix occur in the full sample. **If units from the treated
group are unmatched due to an exact matching restriction, the estimand
no longer corresponds to the ATT.**

### Anti-exact matching (`antiexact`)

Anti-exact matching adds a restriction such that a treated and control
unit with same values of any of the specified anti-exact matching
variables cannot be paired. This can be useful when finding comparison
units outside of a unit’s group, such as when matching units in one
group to units in another when units within the same group might
otherwise be close matches. See examples
[here](https://stackoverflow.com/q/66526115/6348551) and
[here](https://stackoverflow.com/q/61120201/6348551). A similar effect
can be implemented by supplying negative caliper values.

### Matching with replacement (`replace`)

Nearest neighbor matching and genetic matching have the option of
matching with or without replacement, and this is controlled by the
`replace` argument. Matching without replacement means that each control
unit is matched to only one treated unit, while matching with
replacement means that control units can be reused and matched to
multiple treated units. Matching without replacement carries certain
statistical benefits in that weights for each unit can be omitted or are
more straightforward to include and dependence between units depends
only on pair membership. However, it is not asymptotically consistent
unless the propensity scores for all treated units are below .5 and
there are many more control units than treated units ([F. Sävje
2022](#ref-savjeInconsistencyMatchingReplacement2022)). Special standard
error estimators are sometimes required for estimating effects after
matching with replacement ([Austin and Cafri 2020](#ref-austin2020a)),
and methods for accounting for uncertainty are not well understood for
non-continuous outcomes. Matching with replacement will tend to yield
better balance though, because the problem of “running out” of close
control units to match to treated units is avoided, though the reuse of
control units will decrease the effect sample size, thereby worsening
precision ([Austin 2013](#ref-austin2013b)). (This problem occurs in the
Lalonde dataset used in
[`vignette("MatchIt")`](https://kosukeimai.github.io/MatchIt/articles/MatchIt.md),
which is why nearest neighbor matching without replacement is not very
effective there.) After matching with replacement, control units are
assigned to more than one subclass, so the
[`get_matches()`](https://kosukeimai.github.io/MatchIt/reference/match_data.md)
function should be used instead of
[`match_data()`](https://kosukeimai.github.io/MatchIt/reference/match_data.md)
after matching with replacement if subclasses are to be used in
follow-up analyses; see
[`vignette("estimating-effects")`](https://kosukeimai.github.io/MatchIt/articles/estimating-effects.md)
for details.

The `reuse.max` argument can also be used with `method = "nearest"` to
control how many times each control unit can be reused as a match.
Setting `reuse.max = 1` is equivalent to requiring matching without
replacement (i.e., because each control can be used only once). Other
values allow control units to be matched more than once, though only up
to the specified number of times. Higher values will tend to improve
balance at the cost of precision.

### $k$:1 matching (`ratio`)

The most common form of matching, 1:1 matching, involves pairing one
control unit with each treated unit. To perform $k$:1 matching (e.g.,
2:1 or 3:1), which pairs (up to) $k$ control units with each treated
unit, the `ratio` argument can be specified. Performing $k$:1 matching
can preserve precision by preventing too many control units from being
unmatched and dropped from the matched sample, though the gain in
precision by increasing $k$ diminishes rapidly after 4 ([Rosenbaum
2020](#ref-rosenbaum2020)). Importantly, for $k > 1$, the matches after
the first match will generally be worse than the first match in terms of
closeness to the treated unit, so increasing $k$ can also worsen balance
([Rassen et al. 2012](#ref-rassenOnetomanyPropensityScore2012)). Austin
([2010b](#ref-austin2010a)) found that 1:1 or 1:2 matching generally
performed best in terms of mean squared error. In general, it makes
sense to use higher values of $k$ while ensuring that balance is
satisfactory.

With nearest neighbor and optimal pair matching, variable $k$:1
matching, in which the number of controls matched to each treated unit
varies, can also be used; this can have improved performance over
“fixed” $k$:1 matching ([Ming and Rosenbaum 2000](#ref-ming2000);
[Rassen et al. 2012](#ref-rassenOnetomanyPropensityScore2012)). See
[`?method_nearest`](https://kosukeimai.github.io/MatchIt/reference/method_nearest.md)
and
[`?method_optimal`](https://kosukeimai.github.io/MatchIt/reference/method_optimal.md)
for information on implementing variable $k$:1 matching.

### Matching order (`m.order`)

For nearest neighbor matching (including genetic matching), units are
matched in an order, and that order can affect the quality of individual
matches and of the resulting matched sample. With `method = "nearest"`,
the allowable options to `m.order` to control the matching order are
`"largest"`, `"smallest"`, `"closest"`, `"farthest"`, `"random"`, and
`"data"`. With `method = "genetic"`, all but `"closest"` and
`"farthest"` can be used. Requesting `"largest"` means that treated
units with the largest propensity scores, i.e., those least like the
control units, will be matched first, which prevents them from having
bad matches after all the close control units have been used up.
`"smallest"` means that treated units with the smallest propensity
scores are matched first. `"closest"` means that potential pairs with
the smallest distance between units will be matched first, which ensures
that the best possible matches are included in the matched sample but
can yield poor matches for units whose best match is far from them; this
makes it particularly useful when matching with a caliper. `"farthest"`
means that closest pairs with the largest distance between them will be
matched first, which ensures the hardest units to match are given the
best chance to find matches. `"random"` matches in a random order, and
`"data"` matches in order of the data. A propensity score is required
for `"largest"` and `"smallest"` but not for the other options.

Rubin ([1973](#ref-rubin1973)) recommends using `"largest"` or
`"random"`, though Austin ([2013](#ref-austin2013b)) recommends against
`"largest"` and instead favors `"closest"` or `"random"`. `"closest"`
and `"smallest"` are best for prioritizing the best possible matches,
while `"farthest"` and `"largest"` are best for preventing extreme
pairwise distances between matched units.

## Choosing a Matching Method

Choosing the best matching method for one’s data depends on the unique
characteristics of the dataset as well as the goals of the analysis. For
example, because different matching methods can target different
estimands, when certain estimands are desired, specific methods must be
used. On the other hand, some methods may be more effective than others
when retaining the target estimand is less important. Below we provide
some guidance on choosing a matching method. Remember that multiple
methods can (and should) be tried as long as the treatment effect is not
estimated until a method has been settled on.

The criteria on which a matching specification should be judged are
balance and remaining (effective) sample size after matching. Assessing
balance is described in
[`vignette("assessing-balance")`](https://kosukeimai.github.io/MatchIt/articles/assessing-balance.md).
A typical workflow is similar to that demonstrated in
[`vignette("MatchIt")`](https://kosukeimai.github.io/MatchIt/articles/MatchIt.md):
try a matching method, and if it yields poor balance or an unacceptably
low remaining sample size, try another, until a satisfactory
specification has been found. It is important to assess balance broadly
(i.e., beyond comparing the means of the covariates in the treated and
control groups), and the search for a matching specification should not
stop when a threshold is reached, but should attempt to come as close as
possible to perfect balance ([Ho et al. 2007](#ref-ho2007)). Even if the
first matching specification appears successful at reducing imbalance,
there may be another specification that could reduce it even further,
thereby increasing the robustness of the inference and the plausibility
of an unbiased effect estimate.

If the target of inference is the ATE, optimal or generalized full
matching, subclassification, or profile matching can be used. If the
target of inference is the ATT or ATC, any matching method may be used.
When retaining the target estimand is not so important, additional
options become available that involve discarding units in such a way
that the original estimand is distorted. These include matching with a
caliper, matching within a region of common support, cardinality
matching, or exact or coarsened exact matching, perhaps on a subset of
the covariates.

Because exact and coarsened exact matching aim to balance the entire
joint distribution of covariates, they are the most powerful methods. If
it is possible to perform exact matching, this method should be used. If
continuous covariates are present, coarsened exact matching can be
tried. Care should be taken with retaining the target population and
ensuring enough matched units remain; unless the control pool is much
larger than the treated pool, it is likely some (or many) treated units
will be discarded, thereby changing the estimand and possibly
dramatically reducing precision. These methods are typically only
available in the most optimistic of circumstances, but they should be
used first when those circumstances arise. It may also be useful to
combine exact or coarsened exact matching on some covariates with
another form of matching on the others (i.e., by using the `exact`
argument).

When estimating the ATE, either subclassification, full matching, or
profile matching can be used. Optimal and generalized full matching can
be effective because they optimize a balance criterion, often leading to
better balance. With full matching, it’s also possible to exact match on
some variables and match using the Mahalanobis distance, eliminating the
need to estimate propensity scores. Profile matching also ensures good
balance, but because units are only given weights of zero or one, a
solution may not be feasible and many units may have to be discarded.
For large datasets, neither optimal full matching nor profile matching
may be possible, in which case generalized full matching and
subclassification are faster solutions. When using subclassification,
the number of subclasses should be varied. With large samples, higher
numbers of subclasses tend to yield better performance; one should not
immediately settle for the default (6) or the often-cited recommendation
of 5 without trying several other numbers. The documentation for
[`cobalt::bal.compute()`](https://ngreifer.github.io/cobalt/reference/bal.compute.html)
contains an example of using balance to select the optimal number of
subclasses.

When estimating the ATT, a variety of methods can be tried. Genetic
matching can perform well at achieving good balance because it directly
optimizes covariate balance. With larger datasets, it may take a long
time to reach a good solution (though that solution will tend to be good
as well). Profile matching also will achieve good balance if a solution
is feasible because balance is controlled by the user. Optimal pair
matching and nearest neighbor matching without replacement tend to
perform similarly to each other; nearest neighbor matching may be
preferable for large datasets that cannot be handled by optimal
matching. Nearest neighbor, optimal, and genetic matching allow some
customizations like including covariates on which to exactly match,
using the Mahalanobis distance instead of a propensity score difference,
and performing $k$:1 matching with $k > 1$. Nearest neighbor matching
with replacement, full matching, and subclassification all involve
weighting the control units with nonuniform weights, which often allows
for improved balancing capabilities but can be accompanied by a loss in
effective sample size, even when all units are retained. There is no
reason not to try many of these methods, varying parameters here and
there, in search of good balance and high remaining sample size. As
previously mentioned, no single method can be recommended above all
others because the optimal specification depends on the unique qualities
of each dataset.

When the target population is less important, for example, when engaging
in treatment effect discovery or when the sampled population is not of
particular interest (e.g., it corresponds to an arbitrarily chosen
hospital or school; see Mao, Li, and Greene ([2018](#ref-mao2018)) for
these and other reasons why retaining the target population may not be
important), other methods that do not retain the characteristics of the
original sample become available. These include matching with a caliper
(on the propensity score or on the covariates themselves), cardinality
matching, and more restrictive forms of matching like exact and
coarsened exact matching, either on all covariates or just a subset,
that are prone to discard units from the sample in such a way that the
target population is changed. Austin ([2013](#ref-austin2013b)) and
Austin and Stuart ([2015b](#ref-austin2015c), [2015a](#ref-austin2015a))
have found that caliper matching can be a particularly effective
modification to nearest neighbor matching for eliminating imbalance and
reducing bias when the target population is less relevant, but when
inference to a specific target population is desired, using calipers can
induce bias due to incomplete matching ([Rosenbaum and Rubin
1985a](#ref-rosenbaum1985); [Wang 2020](#ref-wang2020)). Cardinality
matching can be particularly effective in data with little overlap
between the treatment groups ([Visconti and Zubizarreta
2018](#ref-visconti2018)) and can perform better than caliper matching
([de los Angeles Resa and Zubizarreta
2020](#ref-delosangelesresaDirectStableWeight2020)).

It is important not to rely excessively on theoretical or
simulation-based findings or specific recommendations when making
choices about the best matching method to use. For example, although
nearest neighbor matching without replacement balance covariates better
than did subclassification with five or ten subclasses in Austin’s
([2009](#ref-austin2009c)) simulation, this does not imply it will be
superior in all datasets. Likewise, though Rosenbaum and Rubin
([1985b](#ref-rosenbaum1985a)) and Austin ([2011](#ref-austin2011a))
both recommend using a caliper of .2 standard deviations of the logit of
the propensity score, this does not imply that caliper will be optimal
in all scenarios, and other widths should be tried, though it should be
noted that tightening the caliper on the propensity score can sometimes
degrade performance ([King and Nielsen 2019](#ref-king2019)).

For large datasets (i.e., in 10,000s to millions), some matching methods
will be too slow to be used at scale. Instead, users should consider
generalized full matching, subclassification, or coarsened exact
matching, which are all very fast and designed to work with large
datasets. Nearest neighbor matching on the propensity score has been
optimized to run quickly for large datasets as well.

## Reporting the Matching Specification

When reporting the results of a matching analysis, it is important to
include the relevant details of the final matching specification and the
process of arriving at it. Using
[`print()`](https://rdrr.io/r/base/print.html) on the `matchit` object
synthesizes information on how the above arguments were used to provide
a description of the matching specification. It is best to be as
specific as possible to ensure the analysis is replicable and to allow
audiences to assess its validity. Although citations recommending
specific matching methods can be used to help justify a choice, the only
sufficient justification is adequate balance and remaining sample size,
regardless of published recommendations for specific methods. See
[`vignette("assessing-balance")`](https://kosukeimai.github.io/MatchIt/articles/assessing-balance.md)
for instructions on how to assess and report the quality of a matching
specification. After matching and estimating an effect, details of the
effect estimation must be included as well; see
[`vignette("estimating-effects")`](https://kosukeimai.github.io/MatchIt/articles/estimating-effects.md)
for instructions on how to perform and report on the analysis of a
matched dataset.

## References

Abadie, Alberto, and Guido W. Imbens. 2006. “Large Sample Properties of
Matching Estimators for Average Treatment Effects.” *Econometrica* 74
(1): 235–67. <https://doi.org/10.1111/j.1468-0262.2006.00655.x>.

———. 2016. “Matching on the Estimated Propensity Score.” *Econometrica*
84 (2): 781–807. <https://doi.org/10.3982/ECTA11293>.

Austin, Peter C. 2009. “The Relative Ability of Different Propensity
Score Methods to Balance Measured Covariates Between Treated and
Untreated Subjects in Observational Studies.” *Medical Decision Making*
29 (6): 661–77. <https://doi.org/10.1177/0272989x09341755>.

———. 2010a. “The Performance of Different Propensity-Score Methods for
Estimating Differences in Proportions (Risk Differences or Absolute Risk
Reductions) in Observational Studies.” *Statistics in Medicine* 29 (20):
2137–48. <https://doi.org/10.1002/sim.3854>.

———. 2010b. “Statistical Criteria for Selecting the Optimal Number of
Untreated Subjects Matched to Each Treated Subject When Using
Many-to-One Matching on the Propensity Score.” *American Journal of
Epidemiology* 172 (9): 1092–97. <https://doi.org/10.1093/aje/kwq224>.

———. 2011. “Optimal Caliper Widths for Propensity-Score Matching When
Estimating Differences in Means and Differences in Proportions in
Observational Studies.” *Pharmaceutical Statistics* 10 (2): 150–61.
<https://doi.org/10.1002/pst.433>.

———. 2013. “A Comparison of 12 Algorithms for Matching on the Propensity
Score.” *Statistics in Medicine* 33 (6): 1057–69.
<https://doi.org/10.1002/sim.6004>.

Austin, Peter C., and Guy Cafri. 2020. “Variance Estimation When Using
Propensity-Score Matching with Replacement with Survival or
Time-to-Event Outcomes.” *Statistics in Medicine* 39 (11): 1623–40.
<https://doi.org/10.1002/sim.8502>.

Austin, Peter C., and Dylan S. Small. 2014. “The Use of Bootstrapping
When Using Propensity-Score Matching Without Replacement: A Simulation
Study.” *Statistics in Medicine* 33 (24): 4306–19.
<https://doi.org/10.1002/sim.6276>.

Austin, Peter C., and Elizabeth A. Stuart. 2015a. “The Performance of
Inverse Probability of Treatment Weighting and Full Matching on the
Propensity Score in the Presence of Model Misspecification When
Estimating the Effect of Treatment on Survival Outcomes.” *Statistical
Methods in Medical Research* 26 (4): 1654–70.
<https://doi.org/10.1177/0962280215584401>.

———. 2015b. “Estimating the Effect of Treatment on Binary Outcomes Using
Full Matching on the Propensity Score.” *Statistical Methods in Medical
Research* 26 (6): 2505–25. <https://doi.org/10.1177/0962280215601134>.

Cohn, Eric R., and José R. Zubizarreta. 2022. “Profile Matching for the
Generalization and Personalization of Causal Inferences.” *Epidemiology*
33 (5): 678. <https://doi.org/10.1097/EDE.0000000000001517>.

de los Angeles Resa, María, and José R. Zubizarreta. 2020. “Direct and
Stable Weight Adjustment in Non-Experimental Studies with Multivalued
Treatments: Analysis of the Effect of an Earthquake on Post-Traumatic
Stress.” *Journal of the Royal Statistical Society: Series A (Statistics
in Society)* n/a (n/a). <https://doi.org/10.1111/rssa.12561>.

Desai, Rishi J., Kenneth J. Rothman, Brian T. Bateman, Sonia
Hernandez-Diaz, and Krista F. Huybrechts. 2017. “A
Propensity-Score-Based Fine Stratification Approach for Confounding
Adjustment When Exposure Is Infrequent:” *Epidemiology* 28 (2): 249–57.
<https://doi.org/10.1097/EDE.0000000000000595>.

Diamond, Alexis, and Jasjeet S. Sekhon. 2013. “Genetic Matching for
Estimating Causal Effects: A General Multivariate Matching Method for
Achieving Balance in Observational Studies.” *Review of Economics and
Statistics* 95 (3): 932945. <https://doi.org/10.1162/REST_a_00318>.

Gu, Xing Sam, and Paul R. Rosenbaum. 1993. “Comparison of Multivariate
Matching Methods: Structures, Distances, and Algorithms.” *Journal of
Computational and Graphical Statistics* 2 (4): 405.
<https://doi.org/10.2307/1390693>.

Hansen, Ben B. 2004. “Full Matching in an Observational Study of
Coaching for the SAT.” *Journal of the American Statistical Association*
99 (467): 609–18. <https://doi.org/10.1198/016214504000000647>.

———. 2008. “The Prognostic Analogue of the Propensity Score.”
*Biometrika* 95 (2): 481–88. <https://doi.org/10.1093/biomet/asn004>.

Hansen, Ben B., and Stephanie O. Klopfer. 2006. “Optimal Full Matching
and Related Designs via Network Flows.” *Journal of Computational and
Graphical Statistics* 15 (3): 609–27.
<https://doi.org/10.1198/106186006X137047>.

Ho, Daniel E., Kosuke Imai, Gary King, and Elizabeth A. Stuart. 2007.
“Matching as Nonparametric Preprocessing for Reducing Model Dependence
in Parametric Causal Inference.” *Political Analysis* 15 (3): 199–236.
<https://doi.org/10.1093/pan/mpl013>.

Hong, Guanglei. 2010. “Marginal Mean Weighting Through Stratification:
Adjustment for Selection Bias in Multilevel Data.” *Journal of
Educational and Behavioral Statistics* 35 (5): 499–531.
<https://doi.org/10.3102/1076998609359785>.

Iacus, Stefano M., Gary King, and Giuseppe Porro. 2012. “Causal
Inference Without Balance Checking: Coarsened Exact Matching.”
*Political Analysis* 20 (1): 1–24. <https://doi.org/10.1093/pan/mpr013>.

King, Gary, and Richard Nielsen. 2019. “Why Propensity Scores Should Not
Be Used for Matching.” *Political Analysis*, May, 1–20.
<https://doi.org/10.1017/pan.2019.11>.

Mao, Huzhang, Liang Li, and Tom Greene. 2018. “Propensity Score
Weighting Analysis and Treatment Effect Discovery.” *Statistical Methods
in Medical Research*, June, 096228021878117.
<https://doi.org/10.1177/0962280218781171>.

Ming, Kewei, and Paul R. Rosenbaum. 2000. “Substantial Gains in Bias
Reduction from Matching with a Variable Number of Controls.”
*Biometrics* 56 (1): 118–24.
<https://doi.org/10.1111/j.0006-341X.2000.00118.x>.

Orihara, Shunichiro, and Etsuo Hamada. 2021. “Determination of the
Optimal Number of Strata for Propensity Score Subclassification.”
*Statistics & Probability Letters* 168 (January): 108951.
<https://doi.org/10.1016/j.spl.2020.108951>.

Rassen, Jeremy A., Abhi A. Shelat, Jessica Myers, Robert J. Glynn,
Kenneth J. Rothman, and Sebastian Schneeweiss. 2012. “One-to-Many
Propensity Score Matching in Cohort Studies.” *Pharmacoepidemiology and
Drug Safety* 21 (S2): 69–80. <https://doi.org/10.1002/pds.3263>.

Ripollone, John E., Krista F. Huybrechts, Kenneth J. Rothman, Ryan E.
Ferguson, and Jessica M. Franklin. 2018. “Implications of the Propensity
Score Matching Paradox in Pharmacoepidemiology.” *American Journal of
Epidemiology* 187 (9): 1951–61. <https://doi.org/10.1093/aje/kwy078>.

Rosenbaum, Paul R. 2010. *Design of Observational Studies*. Springer
Series in Statistics. New York: Springer.

———. 2020. “Modern Algorithms for Matching in Observational Studies.”
*Annual Review of Statistics and Its Application* 7 (1): 143–76.
<https://doi.org/10.1146/annurev-statistics-031219-041058>.

Rosenbaum, Paul R., and Donald B. Rubin. 1985a. “The Bias Due to
Incomplete Matching.” *Biometrics* 41 (1): 103–16.
<https://doi.org/10.2307/2530647>.

———. 1985b. “Constructing a Control Group Using Multivariate Matched
Sampling Methods That Incorporate the Propensity Score.” *The American
Statistician* 39 (1): 33. <https://doi.org/10.2307/2683903>.

Rubin, Donald B. 1973. “Matching to Remove Bias in Observational
Studies.” *Biometrics* 29 (1): 159. <https://doi.org/10.2307/2529684>.

———. 1980. “Bias Reduction Using Mahalanobis-Metric Matching.”
*Biometrics* 36 (2): 293–98. <https://doi.org/10.2307/2529981>.

Sävje, F. 2022. “On the Inconsistency of Matching Without Replacement.”
*Biometrika* 109 (2): 551–58. <https://doi.org/10.1093/biomet/asab035>.

Sävje, Fredrik, Michael J. Higgins, and Jasjeet S. Sekhon. 2021.
“Generalized Full Matching.” *Political Analysis* 29 (4): 423–47.
<https://doi.org/10.1017/pan.2020.32>.

Sävje, Fredrik, Jasjeet Sekhon, and Michael Higgins. 2018. *Quickmatch:
Quick Generalized Full Matching*.
<https://CRAN.R-project.org/package=quickmatch>.

Schafer, Joseph L., and Joseph Kang. 2008. “Average Causal Effects from
Nonrandomized Studies: A Practical Guide and Simulated Example.”
*Psychological Methods* 13 (4): 279–313.
<https://doi.org/10.1037/a0014268>.

Sekhon, Jasjeet S. 2011. “Multivariate and Propensity Score Matching
Software with Automated Balance Optimization: The Matching Package for
R.” *Journal of Statistical Software* 42 (1): 1–52.
<https://doi.org/10.18637/jss.v042.i07>.

Stuart, Elizabeth A. 2008. “Developing Practical Recommendations for the
Use of Propensity Scores: Discussion of ‘A Critical Appraisal of
Propensity Score Matching in the Medical Literature Between 1996 and
2003’ by Peter Austin,Statistics in Medicine.” *Statistics in Medicine*
27 (12): 2062–65. <https://doi.org/10.1002/sim.3207>.

———. 2010. “Matching Methods for Causal Inference: A Review and a Look
Forward.” *Statistical Science* 25 (1): 1–21.
<https://doi.org/10.1214/09-STS313>.

Stuart, Elizabeth A., and Kerry M. Green. 2008. “Using Full Matching to
Estimate Causal Effects in Nonexperimental Studies: Examining the
Relationship Between Adolescent Marijuana Use and Adult Outcomes.”
*Developmental Psychology* 44 (2): 395–406.
<https://doi.org/10.1037/0012-1649.44.2.395>.

Thoemmes, Felix J., and Eun Sook Kim. 2011. “A Systematic Review of
Propensity Score Methods in the Social Sciences.” *Multivariate
Behavioral Research* 46 (1): 90–118.
<https://doi.org/10.1080/00273171.2011.540475>.

Visconti, Giancarlo, and José R. Zubizarreta. 2018. “Handling Limited
Overlap in Observational Studies with Cardinality Matching.”
*Observational Studies* 4 (1): 217–49.
<https://doi.org/10.1353/obs.2018.0012>.

Wan, Fei. 2019. “Matched or Unmatched Analyses with
Propensity-Scorematched Data?” *Statistics in Medicine* 38 (2): 289–300.
<https://doi.org/10.1002/sim.7976>.

Wang, Jixian. 2020. “To Use or Not to Use Propensity Score Matching?”
*Pharmaceutical Statistics*, August. <https://doi.org/10.1002/pst.2051>.

Zakrison, T. L., Peter C. Austin, and V. A. McCredie. 2018. “A
Systematic Review of Propensity Score Methods in the Acute Care Surgery
Literature: Avoiding the Pitfalls and Proposing a Set of Reporting
Guidelines.” *European Journal of Trauma and Emergency Surgery* 44 (3):
385–95. <https://doi.org/10.1007/s00068-017-0786-6>.

Zubizarreta, José R., Ricardo D. Paredes, and Paul R. Rosenbaum. 2014.
“Matching for Balance, Pairing for Heterogeneity in an Observational
Study of the Effectiveness of for-Profit and Not-for-Profit High Schools
in Chile.” *The Annals of Applied Statistics* 8 (1): 204–31.
<https://doi.org/10.1214/13-AOAS713>.
