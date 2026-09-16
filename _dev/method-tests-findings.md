# Findings from adding tests for the untested matching methods

Behavior-recording tests were added for `method = "exact"`, `"subclass"`, `"optimal"`,
`"quick"`, `"cardinality"`, `"cem"`, `"full"`, and `NULL`, and for `add_s.weights()`,
`summary()`, and `match_data()` (342 tests across eleven files), plus test helpers
taken from WeightIt and cobalt. Everything below was found while writing them and verified by running code.

Environment: R 4.6.1; optmatch 0.10.8, quickmatch 0.2.3, highs 1.14.0.2,
Rglpk 0.6.5.1, cobalt 5.0.0.

Items 1, 2, 4, 6, 7, 8, 9, and 10 from the original version of this document, plus the
`add_s.weights()` bug and the CRAN-safety problem raised later, have been resolved and
removed from the list below; the verification for each is summarized at the bottom. What
follows is what is still open.

## Open issues

### 1. RESOLVED: `method = "cardinality"` and non-constant `s.weights`

**Cause.** The `cardinality` and `profile_ate` branches constrained the *weighted* group
sizes to be exactly equal (up to `ratio`) with `==` rows on a continuous slack. With
constant weights -- all exactly 1 after `.make_sum_to_n()` -- that is an equality between
integer counts, and HiGHS solves at the root in 1 node. With non-constant weights it is an
exact equality between two real-valued weighted subset sums, which almost no selection
satisfies: 14595 nodes and a 33% gap (incumbent 176.75, dual bound 235.71) after 20 s. It
was erratic because difficulty is not monotone in how much the weights vary -- 0.25 s at
+/-1e-4, timeout at +/-1e-3, 3.4 s at +/-1e-2, timeout at +/-5e-2. `profile_att` was never
affected, having no cross-group equality.

**Resolution**, following the decision that sampling weights only make sense against a
fixed target:

1. `s.weights` is rejected for cardinality matching, which balances the groups against
   each other and so has no fixed target population to generalize to. Profile matching
   (`estimand = "ATE"`, or `ratio = NA`) does have one and keeps the feature.
2. In `profile_ate`, the objective and the `ratio` constraint now count units rather than
   summing weights, so `s.weights` enters only the balance rows where it weights the
   covariate means. The per-group weighted-size slacks are kept: they are what make the
   balance rows an exact linearization of the weighted mean, and because they are
   *definitions* of each group's weighted size rather than equalities linking the groups,
   they cost nothing.

Verified: every `s.weights` case now solves to optimality inside the default `time`
(2-30 s, erratic but reliable), all group deviations land within `tols/2` as the
formulation promises, the `ratio` constraint holds exactly on the unweighted counts, and
unweighted results are unchanged -- every existing snapshot still passes, which follows
because with unit weights `rep(1, n)` and `s.weights` are the same vector.

**Tried and rejected**, in case it comes up again:

- Count-based sizing while keeping group-vs-group balance, i.e. retaining cardinality
  matching and just swapping the denominators: fast (1 node) but silently loosens the
  tolerance. The balance row `sum_1 sw X z - sum_0 sw X z <= tol * n1` is the weighted
  *mean* difference only when `n1` is the weighted size of both groups; with a count
  denominator the error scales with each covariate's mean-to-SD ratio (5.15 for `educ` on
  lalonde), giving a realized SMD of 0.070-0.090 against a 0.05 request.
- Centering the covariates to kill that error term: destroys the LP structure, 18097 nodes
  even *without* `s.weights`.
- Bounding the two weighted sizes to within a band of one another: needs a ~0.1% band to
  control a mean/SD of 5, which is the knife-edge equality again.
- Two slacks with per-group weighted denominators and group-vs-group balance: not linear,
  since that is a difference of two ratios with different variable denominators.
- Tightening the bounds on the weighted-size and count slacks: valid, but made no
  measurable difference -- presolve already derives them. Kept anyway, being correct.

The reason no formulation is exact, group-vs-group, and easy at the same time: linearity
requires either one shared denominator, which forces the exact cross-group weighted
equality, or constant denominators, which means a fixed target and therefore profile
matching. There is no third option, which is why the resolution is a scope decision rather
than a reformulation.

Still outstanding, independent of all this: the highs path *errors* on a timeout where
`glpk` and `gurobi` warn and return the feasible incumbent, and its message advises
increasing `time`, which was misleading under the old formulation.

### 2. `method = "optimal"` inherits optmatch's loose default `tol`

Optimal matching's advantage over greedy matching is that it minimizes the total
within-pair distance. On lalonde it still loses:

| | total within-pair \|PS\| distance |
|---|---|
| `method = "nearest"` | 39.69279872 |
| `method = "optimal"`, default `tol = 1e-3` | 39.69374188 |
| `method = "optimal"`, `tol = 1e-6` | 39.69279872 |
| `method = "optimal"`, `tol = 1e-9` | 39.69279872 |

MatchIt passes `tol` through to optmatch but does not override its default of `1e-3`,
which on a problem this small is loose enough to return a measurably suboptimal
solution — one a user comparing methods would notice. Tightening the default costs
nothing measurable in runtime here. Worth considering, since the current default makes
the method fail its own headline property on the package's example dataset.

## `summary()` and `match_data()`: findings from writing their tests

23 tests for `summary.matchit()` and `summary.matchit.subclass()`, 22 for
`match_data()`/`match.data()`/`get_matches()`. Two real bugs, both now fixed.

### FIXED: `match_data()` silently ignored an invalid `data` argument

`match_data()` looked for the dataset in four places and took the first candidate that was
two-dimensional with as many rows as the `matchit` object. The user's own `data` argument
was only the *first* candidate, so when it failed that test it was discarded and the loop
moved on to `eval(object$call$data, envir = environment(object$formula))` -- the original
dataset. Both documented errors were therefore unreachable whenever the original data was
still in scope, and the user got a result computed on a dataset they did not pass, with no
warning. Removing the original from the calling frame did not help, because the formula
environment still reached it.

Fixed by splitting the two cases: a supplied `data` is now authoritative and is validated
(two-dimensional, coerced from a matrix if need be, and with exactly as many rows as there
were units, the error reporting both counts), and the four-place search runs only when
`data` is missing. `get_matches()` inherits the fix, since it resolves its data through
`match_data()`. The documentation for `data` now says it must be the original dataset and
that it is usually best left unspecified, and the Note says what goes wrong if something
else is passed.

A related but separate point, recorded rather than flagged and still true: the rows of a
supplied `data` are assumed to line up with the original and this is not checked, so a
reordered frame of the right size is accepted and silently mismatched. That cannot be
detected in general, so the assumption is just worth stating; it is pinned as such.

The same defect was present, in weaker forms, everywhere else a dataset is accepted
alongside a `matchit` object, so all four now share one validator, `.check_supplied_data()`
in `aux_functions.R`:

- `summary()` ran the same four-place search and replaced a wrongly sized `data` with the
  recovered original, which surfaced only as an error about `addlvariables` and only when
  the resulting variables happened to be the wrong length.
- `plot()` validated `data` only inside the branch that consumes `which.xs`, so
  `plot(m, data = anything)` ignored the argument entirely.
- `add_s.weights()` already validated it; only its message changed.

The validator reports both row counts and coerces a matrix. `match_data()` phrases it as
"must be the original dataset"; the other three, where `data` is a source of additional
variables rather than the dataset itself, phrase it as "must have one row for each of the
N units".

### FIXED: `addlvariables` as a matrix

`.process_X()` explicitly admitted a matrix -- `!is.matrix(x) && !is.data.frame(x)` was
the rejection test -- but a matrix then failed several frames down inside
`get_covs_matrix()` with base R's `'data' must be a data.frame, not a matrix or an array`,
which names the wrong argument and describes a restriction the documentation does not
state. Fixed by coercing a matrix to a data frame at the point where it is already
accepted, so matrices work as the guard intended; the roxygen for `addlvariables` now
mentions matrices. Every remaining error on that path names `addlvariables`.

While there: the string branch tested `is.character(addlvariables)`, which is also true of
a *character matrix*, so a matrix of covariate values was interpreted as a vector of
variable names. It now requires `is_null(dim(addlvariables))`, and a character matrix is
treated as covariates, its columns becoming factors exactly as character columns of a data
frame do.

### Non-findings worth recording

- **`improvement` controls whether `reduction` is computed at all**, not just whether it
  prints: `summary(m)$reduction` is `NULL` and `summary(m, improvement = TRUE)$reduction`
  is a matrix. I had assumed it was print-only and checked.
- **`Std. Pair Dist.` comes from the strata, not from `match.matrix`.** It is reported for
  `cem`, `exact`, and `full`, none of which produce a `match.matrix` in these
  specifications, and is `NA` only for `cardinality`, which produces neither strata nor
  pairs. Under `exact` every within-stratum distance is exactly zero, since the units are
  identical on the covariates; that is now a test.
- **`which.subclass` is not a `summary()` argument.** The user-facing argument is
  `subclass`, which takes `TRUE`, `FALSE`, or a vector of indices; `which.subclass` is the
  internal name it is assigned to, and also a real argument to `plot.matchit()`. Passing
  it to `summary()` is silently absorbed by `...`.
- **The `nn` table reconciles.** Matched, unmatched, and discarded counts sum to the
  sample size, each row matches a direct count from the object, and the ESS rows equal the
  raw counts exactly when there are no sampling weights and matching is 1:1.
- **`get_matches()` differs from `match_data()` in the way it should.** With replacement a
  reused control gets one row per pair, sharing an `id` and distinguished by `subclass`;
  every subclass holds exactly one treated unit; and it errors informatively for methods
  with no `match.matrix`.

The balance matrices are pinned rather than the printed tables, since printed output
depends on `getOption("digits")` and `OutDec`. Printing is checked structurally instead,
by asserting which blocks appear under `un`, `improvement`, and `subclass`.


## CEM and full matching: findings from writing their tests

`method = "cem"` has the largest option surface of any method here, and now has 43 tests
covering it; `method = "full"` has 28. Nothing found was serious, but four things are
worth knowing.

### 1. `k2k.method = NULL` is not random

The documentation describes `NULL` as "for random matching". It is not: with the default
`m.order = "data"` it sets every distance to zero and matches in data order, giving the
same answer on every run and under every seed. Randomness comes from
`m.order = "random"`, not from the absence of a distance. Verified both ways -- `NULL`
alone is seed-invariant, `NULL` plus `m.order = "random"` is not. The doc wording is worth
changing to something like "no distance (matching in the order given by `m.order`)".

### 2. The `grouping` error message named the wrong argument -- FIXED

`matchit2cem.R:357` read `{.arg groupings}` where the argument is `grouping`. Fixed as a
one-line change; no test depended on the misspelling.

### 3. `cutpoints = NA` works but is undocumented

A `cutpoints` entry of `NA` reaches the same branch as `0`: the variable is not binned, so
exact matching happens on it. Verified identical to `cutpoints = 0` on the same variable.
The documented forms are a cutpoint vector, a bin count, or an algorithm string, so `NA`
is an accident that happens to behave sensibly. Either document it or reject it.

### 4. `k2k.method`, `m.order`, and `mpower` are silently ignored without `k2k = TRUE`

Every argument that belongs to another method (`distance`, `exact`, `mahvars`, `discard`,
`replace`, `caliper`, `ratio`) warns when supplied to `method = "cem"`. The three that
belong to `cem` but only apply when `k2k = TRUE` are dropped without comment, so
`matchit(..., method = "cem", k2k.method = "manhattan")` quietly does plain CEM. Recorded
in a test rather than treated as a bug, since it is a consistency question rather than a
wrong result.

### Non-findings worth recording

- **Every `k2k.method` value does something.** All nine named methods give pairwise
  distinct pairings, and `minkowski` with the default `mpower = 2` correctly coincides
  with `euclidean`. A test asserts the pairwise distinctness so that an option cannot
  silently become a no-op.
- **`mpower` plateaus, which is expected.** Nine values from 0.5 to 50 give six distinct
  pairings: 3, 4, and 6 agree, and 10 and 50 agree. The within-stratum ranking of
  candidate controls stabilizes as the power grows. `mpower = 50` is still not identical
  to `k2k.method = "maximum"`, which it approaches only in the limit. The test therefore
  compares well-separated powers rather than adjacent ones.
- **`s.weights` reaches the k2k distance only for the scaled methods.** It changes the
  pairing for `mahalanobis`, `robust_mahalanobis`, and `scaled_euclidean`, and not for
  `euclidean` or `manhattan`, which have no scaling factor to weight. That matches the
  documentation exactly and is now pinned in both directions.
- **Documented equivalences hold.** `cutpoints = 0` on every numeric variable reproduces
  `method = "exact"` on the same formula, weights included; values on a bin boundary go
  into the higher bin; a length-1 numeric is a bin count while a longer vector is a set of
  cut points; `"q4"` gives exactly equal-sized bins; levels omitted from `grouping` keep
  their own category; and every k2k pair comes from a single coarsened stratum.
- **Everything is deterministic** except `m.order = "random"`, for both methods, so all of
  it is safe to pin.

For `full`, the constraint arguments were checked against their meaning rather than only
for running: `min.controls` forces at least that many controls per treated unit in every
stratum, `max.controls` caps it, tightening either produces more strata, `mean.controls`
and `omit.fraction` each drop units and cannot be combined (optmatch rejects it), and a
caliper with `min.controls` is rejected by MatchIt. One cross-method invariant is also
pinned: mean within-stratum distance under full matching is no worse than under 1:1
optimal matching, which is the property that motivates the method.


## `print()`, `plot()`, and `rbind()`: findings from writing their tests

10 tests for `print.matchit()`, 16 for `plot.matchit()`/`plot.matchit.subclass()`/
`plot.summary.matchit()`, 11 for `rbind.matchdata()`/`rbind.getmatches()`. Four bugs, all
fixed; none of them touch a number the package reports.

### FIXED: `plot()` on a subclassification object hung in non-interactive sessions

`plot.matchit.subclass()` entered its subclass-selection menu whenever `interactive = TRUE`
(the default) and `subclass` was missing, without checking whether the session was actually
interactive. The menu loop is `while (!ans %in% 0:k) { message; readline() }`, and
`readline()` returns `""` immediately in a non-interactive session, so the loop never
terminated: 26,000 lines of output in 12 seconds and no way out but a kill signal. This is
worse than an error -- it hangs `R CMD check`, knitr, and any script that plots a
subclassification fit without naming a subclass.

Fixed by gating the menu on `interactive()`, so a non-interactive session falls through to
the aggregate plot, which is what the argument documentation already describes. There is a
test that would hang rather than fail if the guard were removed, which is worth knowing
before anyone "simplifies" it.

### FIXED: `print()` left the distance line unterminated

The newline that ended the ` - distance:` line was emitted inside the branch for
propensity-score-like distances, so matching on the Mahalanobis distance or on a
user-supplied one produced

    - distance: Robust Mahalanobis - number of obs.: 614 (original), 370 (matched)

And because the "estimated with" line began with its own newline, a distance that had both a
bracketed annotation and an estimating method produced a blank line in the middle of the
block. Both come from the same cause -- line termination distributed across the branches
that print the pieces -- and are fixed by terminating the line once, where it ends.

### FIXED: the bracketed annotation never said "matching"

`print.matchit()` computed `nm <- is_null(x[["method"]])`, but `matchit` objects have no
`method` component; the method is in `x$info$method`. `nm` was therefore always `TRUE`, and
the `matching` and `subclassification` labels in the bracketed list of what the distance
measure was used for were dead code. A propensity score used for matching with a caliper on
itself printed `Propensity score [caliper]` instead of `Propensity score [matching,
caliper]`. Fixed by reading `info$method`. The annotation only ever appears when there is a
caliper or a common-support restriction, so this changes the output of those cases and no
others.

### FIXED: a `var.order` error named the wrong value

`plot.summary.matchit()` rejects `var.order = "unmatched"` when the summary holds no
unmatched statistics, but said "if `un = TRUE` in the call to `summary()`" when it is
`un = FALSE` that causes it.

### Non-findings worth recording

- **A character subclass index is caught, but not where you would expect.**
  `plot(m, subclass = "2")` passes `all(subclass %in% seq_along(subclasses))` by coercion,
  then indexes the level vector by name, gets `NA`, and fails one level down with "the
  argument supplied to `subclass` is not the index of any subclass". Informative enough
  that it is pinned rather than flagged.
- **`type = "jitter"` with `interactive = TRUE` is harmless non-interactively.** It prints
  its "use first mouse button" instruction and `identify()` returns immediately on a
  non-interactive device. Only the subclass menu had the blocking problem.
- **`rbind()` renames by position, not by precedence.** The output column names are taken
  from the first input that has each one, so `rbind(a, b)` and `rbind(b, a)` can disagree on
  whether the distance column is called `distance` or `prop.score`. Documented behavior,
  and pinned in both directions.
- **`rbind()` output can be stacked again**, and subclass labels accumulate prefixes
  (`1_1_1`). Uniqueness is preserved, which is all that matters.
- **An input lacking a column gets `NA` for it**, so a CEM dataset with no distance can be
  stacked with a propensity-score one; the column exists and is `NA` for the CEM rows.
- **`par()` is restored** by `plot.summary.matchit()`, which sets several parameters via
  `dotchart()`. Checked by comparing `par(no.readonly = TRUE)` across the call.

## Verified as resolved

Checked against the working-tree diff by re-running each reproduction:

- **`subclass = 0.5` segfault** — fixed. The rewritten input check in
  `matchit2subclass()` reads a lone value in (0, 1) as a quantile, so `subclass = 0.5`
  now yields 2 subclasses identical to `subclass = c(0, 0.5, 1)`. Duplicated quantiles
  are collapsed. `has_n_unique()` is additionally passed `as.integer(round(subclass))`,
  and the R-side branch guarantees `subclass >= 2` at that call, so the non-integer `n`
  path is unreachable. Now covered by two tests.
  - Residual, not reachable from `matchit()`: `has_n_unique()` itself is still
    unguarded. `Vector<RTYPE> seen(n); seen[0] = x[0]` writes out of bounds for any
    `n < 1`, and `Rcpp::Vector::operator[]` does not bounds-check. All six call sites
    are now safe (five pass `2L`), so this is latent robustness only — worth a guard
    whenever `has_n_unique.cpp` is next touched.
- **`subclass = 1` silently accepted** — fixed; it now errors, along with `0`,
  `c(0, 1)`, `2.5`, `-1`, `NA`, and `"a"`. Covered by a test over all seven.
- **SYMPHONY non-determinism** — resolved by removing the solver. Verified gone from
  `match_arg()`, `check_installed()`, `dispatch_optimizer()`,
  `cardinality_error_report()`, the solution-extraction switch, `DESCRIPTION` Suggests,
  the roxygen block, and `man/method_cardinality.Rd`. `solver = "symphony"` now errors;
  a test pins that, and the vignette and CI workflows were cleaned up too.
- **"alotted"** — fixed; no occurrences remain in `R/` or `man/`.
- **`s.weights` inert for stratification methods** — fixed, and the implementation is
  correct. `get_weights_from_subclass()` now takes `s.weights` and builds stratum masses
  from weighted rather than raw counts. Verified two ways: the within-stratum ratio of
  weighted control mass to weighted treated mass is constant across strata for `exact`
  (ATT and ATC), `subclass`, `full`, `quick`, and `cem` — which is the property that
  must hold, since `matchit()`'s per-group rescaling multiplies that ratio by one global
  constant but cannot make it vary between strata — and constant `s.weights` (1 and 3)
  reproduce the previous unweighted weights exactly for all five methods. Covered by
  three new tests in `test-method_exact.R`, including a reusable
  `expect_stratum_mass_balanced()`. The `add_s.weights()` half of the feature needed a
  separate fix; see the next entry.
- **`add_s.weights()` never recomputed the matching weights** — fixed. The guard now
  reads `m$info$method`, and the recomputed weights go through the same normalization
  `matchit()` applies, gated on the `normalize` value from the original call. `normalize`
  is now recorded in `info` rather than recovered from `m$call`, which would fail
  whenever it was passed as a variable (`matchit(..., normalize = nrm)` records the
  symbol `nrm`, unresolvable later). Objects predating the `info` field default to
  `TRUE`. Verified: `add_s.weights()` reproduces `matchit(s.weights = )` exactly for
  `exact` and `cem`, whose strata do not depend on the propensity score; for `subclass`,
  `full`, and `quick` the two routes differ only because `s.weights` changes the score
  and therefore the strata, and given identical strata they agree exactly, under both
  `normalize = TRUE` and `FALSE`. The weighted stratum masses balance in every case, and
  the non-stratification methods are untouched. Covered by a new
  `test-add_s.weights.R` (16 tests), which is also the first coverage this exported
  function has had.
- **`subclass = Inf`, `NaN`, and `NA_real_`** -- fixed. A finiteness check was added to
  the count branch, and `all(is.finite(subclass))` to the quantile branch, which
  short-circuits before the `all(subclass >= 0 & subclass <= 1)` that returned `NA` for
  `NaN` and `NA_real_` and left the `if` with nothing to branch on. The padding was also
  simplified to `sort(unique(c(0, 1, subclass)))`. Verified across 27 input forms: every
  invalid one now produces the intended message and none reaches base R; `0.5` still
  equals `c(0, 0.5, 1)`, duplicates and unsorted input still collapse correctly, and the
  integer-count path is unchanged.
- **Stale SYMPHONY reference in the vignette** -- fixed, and also removed from the
  pkgdown and rhub workflows, which still passed `Rsymphony=?ignore`.
- **Snapshots not CRAN-safe** — fixed by adding `skip_on_cran()`, keeping the tests in
  the build as intended. The skip lives inside `expect_matchit_snapshot()` rather than at
  the top of each test, so the structural checks preceding it still execute on CRAN;
  verified that a failure recorded before a skip is still reported as a failure, so this
  does not hide regressions. The one bare `expect_snapshot_value()` call not routed
  through that helper ("estimated propensity scores are stable across links") got an
  explicit `skip_on_cran()`. With `NOT_CRAN` unset: 148 tests skip, 0 fail, and 5922
  expectations still run, including all of `test-add_s.weights.R`. No `.new.md` files are
  written in that mode. `waldo` added to Suggests, since `helpers.R` calls
  `waldo::compare()` and the tests now ship.
  - Residual: the snapshot files (1.4 MB) still go into the tarball without ever being
    read on CRAN. Adding `^tests/testthat/_snaps$` to `.Rbuildignore` would drop the dead
    weight while keeping the tests shipping. Safe only because every snapshot test skips
    on CRAN — if one ever ran without its `_snaps` file it would silently record rather
    than compare.
- **`method = "optimal"` + `antiexact` internal error** — fixed on the MatchIt side. The
  `allNA(pair)` guard now gives "No matches were found." instead of letting an empty
  result reach `subclass2mmC()`. optmatch still judges this particular restriction
  impossible while `method = "full"` matches the same specification successfully, so if
  that verdict is worth chasing it is a separate question from the crash. Test updated
  to assert the informative error.
- **`glpk` slower than `highs`** — accepted as inherent to the solver; `highs` is
  already the default.
- **Literal `tc[1L]` in a warning** — fixed earlier as a one-line change; the corrected
  message is pinned in `test-method_optimal.R`.

## Test helpers added

Ported from `WeightIt/tests/testthat/helpers.R` and adapted:

- `expect_no_unexpected_warning(expr, known)` — asserts no warning other than the
  listed ones. For specifications that warn only on some samples or solver versions,
  where `expect_warning()` would be wrong in both directions.
- `expect_not_equal()` — `expect_equal()`'s complement, via `waldo`, for asserting
  something actually changed. Used for the estimand and `s.weights` comparisons.
- `expect_balance_improved(m)` — the weakest useful check that matching did something,
  via `cobalt::col_w_smd()`. Catches sign errors and mixed-up weights that both
  structural checks and snapshots pass over.
- `inject_missingness(data, cols, prop, seed)` — sets a fixed proportion of columns to
  `NA` without disturbing the calling test's RNG stream.

WeightIt's `expect_ATT_weights_okay()` and `expect_M_parts_okay()` were not ported;
both are specific to weighting output that `matchit` objects do not have.

Taken from cobalt's test helpers, replacing an interim `expect_matchit_condition()` and
the `.w()` whitespace-regex helper, both now removed:

- `expect_err()`, `expect_wrn()`, `expect_msg()`, and `squish()` -- for the conditions
  `arg::err()`, `arg::wrn()`, and `arg::msg()` raise. They collapse whitespace in the
  *observed* message and then match a literal substring, which sidesteps both cli
  problems at once: the hard wrapping that defeats a `fixed = TRUE` match against a long
  message, and the metacharacters cli introduces that a regex would have to escape.
  Matching is case-insensitive, which matters more than it looks: cli capitalizes the
  first letter, so the lowercase source text in `arg::err("covariates must be ...")`
  would not otherwise match the rendered "Covariates must be ...". Warnings and messages
  are caught with `withCallingHandlers()` so the expression runs to completion, which
  removed the duplicate `suppressWarnings()` re-runs several tests needed in order to
  reach the fitted object afterward, and output goes through `capture.output()` so a
  `print()` under test does not bury the reporter.

One limitation worth knowing: `expect_wrn()` keeps only the *first* warning, so a second
one cannot be asserted through it. The five tests that pin two warnings from a single
call therefore keep testthat's native nested `expect_warning()`, which re-signals what
it does not match. Everything else uses the cobalt helpers.

Also taken from cobalt, when the plotting tests needed it:

- `local_null_device()` -- sends plots to a null `pdf()` device for the rest of the
  calling test, registering `dev.off()` as an `on.exit()` expression in the caller's
  frame so the device closes even when an expectation fails. Verbatim from cobalt, so
  the two packages' copies stay interchangeable.

Not ported from cobalt: `local_cobalt_options()`, `local_bal_tab_snapshot()`, and
`skip_if_slow()`, which are specific to cobalt's options and snapshot surface.

## Testing notes

- Solver-backed methods pin values that are a joint property of MatchIt and the
  installed solver, so a snapshot failure after an optmatch, quickmatch, or highs
  upgrade is expected rather than a MatchIt regression. Noted in each file's header.
  See item 4 for why this now matters more than it did.
- An error raised inside `test_that()` is reported as an error, not a failure. Summaries
  that only count `failed` will call a broken suite green — check `error` too.
- Full suite: 342 tests, 7533 passing expectations, 0 failures, 0 errors. With
  `NOT_CRAN` unset: 187 skipped, 6682 expectations still run, 0 failures.

## Still unpinned

- `method = "genetic"`, deliberately skipped.
- The rendered content of any plot. The plotting tests establish that each type and
  argument combination runs on a null device and that nothing prompts for input; what is
  actually drawn is not checked, and would need vdiffr or an equivalent to be.
