# `src/` cleanup: what changed and what it bought

Companion to `_dev/cpp-cleanup-bench.R`. Every change below is behavior-preserving;
the 381-test suite passes with **no snapshot churn**, which is the evidence for that
claim, since the snapshots pin `match.matrix`, `weights`, and `subclass` for every
method. `devtools::check()` is clean (0 errors, 0 warnings, 0 notes).

## Measurement trap, recorded because it cost an hour

`pkgbuild::compile_dll()` compiles at **`-O0`** by default. The first before/after
comparison was nonsense: the "before" numbers came from the committed `.so` (built
at `-O2` by `R CMD INSTALL`) and the "after" numbers from a debug build, so
everything looked 30-60% *slower*. Worse, `compile_dll()` only recompiles changed
files, so a tree can end up with a mix of `-O0` and `-O2` objects and produce timings
that move when an unrelated file is touched.

Whenever these benchmarks are re-run: `pkgbuild::clean_dll()` then
`pkgbuild::compile_dll(debug = FALSE)`, on **both** sides. `_dev/cpp-cleanup-bench.R`
takes `PKG_DIR` from the environment so a baseline can be benchmarked from a
`git worktree` of the commit being compared against.

## Timings (both sides clean `-O2`, min of 3 runs, n = 20,000 unless noted)

| case | before | after | change |
|---|---|---|---|
| `nn_matchC_mahcovs`, lalonde | 0.014 | 0.005 | **-64%** |
| `nn_matchC_mahcovs` | 1.539 | 0.958 | **-38%** |
| `nn_matchC_mahcovs_closest` | 2.086 | 1.094 | **-48%** |
| `reuse.max` | 0.032 | 0.029 | -9% |
| everything else (nearest, optimal, cem, full, quick, `pair.dist`) | | | within ±3% |

The Mahalanobis win is one change: `euc_dist_sq()` took two `NumericVector`s, so
`euc_dist_sq(mah_covs.row(i), mah_covs.row(j))` allocated two R vectors for every
candidate control considered. It now takes the matrix and two row indices.

A raw-pointer variant (`&x[0]` plus manual `i + k * nrow`) was tried and was another
~2% faster, i.e. inside the noise, so the readable `x(i, k)` version was kept. Worth
knowing that `Rcpp::Matrix::operator()` bounds-checks and throws, so it is not free;
if a future hot loop needs it, the pointer version is the escape hatch.

The other allocation removals (the `mm(row, sum(!is_na(mm(row, _))))` idiom, the
`na_omit()` in `mm_okay()`, the row copies in `antiexact_okay()`, the per-unit
`max(as<IntegerVector>(n_eligible[g_c]))`) did not move the needle measurably: they
run once per *match written* or per *treated unit*, not per candidate, and the
candidate scan dominates. They are still worth having, and the `antiexact_okay()` one
does run per candidate whenever `antiexact` is used.

## Alternatives to `base::order()`: tested, rejected

The matching code calls `base::order()` from C++ in seven places. Four candidate
replacements were written and checked against `order()` on 17 adversarial vectors
(all-equal, two-valued, presorted, reverse-sorted, `NA`/`NaN` mixes, +/-`Inf`, `-0`,
denormals, length 0 and 1) and 300 random fuzz cases weighted toward ties and
missingness, ascending and descending. **All five implementations agree with
`base::order()` byte for byte**, so the choice is purely about speed. See
`_dev/order-alternatives.{cpp,R}`, which re-runs both halves.

Microseconds per call:

| | rcall (current) | R_orderVector1 | std::sort on indices | std::stable_sort | std::sort on pairs |
|---|---|---|---|---|---|
| n=614 continuous | 12.5 | **5.0** | 15.5 | 23.5 | 17.5 |
| n=20,000 continuous | **350** | 1150 | 1180 | 1700 | 1320 |
| n=20,000 discrete | **140** | 1440 | 1260 | 1230 | 1310 |
| n=20,000 presorted | **10** | 80 | 70 | 200 | 100 |
| n=100,000 continuous | **1333** | 8833 | 7500 | 10867 | 8100 |
| n=8,000 (heap_ord shape) | **150** | 275 | 410 | 790 | 465 |

**`base::order()` stays.** It wins by 3-8x everywhere the sort is a measurable part of
the work, because `order()` on a double vector defaults to a **radix** sort: measured
directly, `order(x, method = "radix")` is 320 us at n = 20,000 against 1160 us for
`method = "shell"`. `R_orderVector1()`, the C API entry point, *is* the shell sort, and
its 1150 us matches. The comparison sorts lose for the same reason, and sorting
`(value, index)` pairs for locality does not rescue them: the `ISNAN` checks in the
comparator, not the moves, are the cost.

The one place a C++ version wins is n = 614, where `R_orderVector1()` saves 7.5 us by
skipping the closure call. A whole lalonde-sized `matchit()` takes ~4 ms, so that is
0.2% of one call, against carrying a second code path and a hand-written NA ordering.
Not worth it.

What did come out of this: the comments at all seven call sites now say *why*
`base::order()` is there, with the measurement referenced, so the next person to look
does not spend the afternoon proving it again. Writing our own radix sort would be the
only way to beat it, and getting the `NA`/`NaN`/`-0` tie order byte-identical by hand
is exactly the kind of thing R has already got right.

## Bugs fixed

1. **`Function o("order")` resolved in the global environment.** `Rcpp::Function(name)`
   looks up with `Rf_findFun` from `R_GlobalEnv`, so a user function named `order`
   would be called instead of `base::order()`. Seven sites, all now bound with
   `Environment::base_env()["order"]`. This is the only fix here a user could observe,
   so it is the one with a NEWS entry.
2. **`std::isfinite()` on an integer NA** in `subclass2mmC()`: `NA_INTEGER` is
   `INT_MIN`, which is finite as a double, so the guard never fired. Masked by the
   `si != ss[s]` test on the next line, so no wrong results — but the check was dead.
3. **`preprocess_matchC()` reserved `n1 * n0`** ints up front (~100 MB at 5,000 x
   5,000) in the function whose purpose is to avoid the dense product.
4. **`has_n_unique_()`** read `x[0]` and wrote `seen[0]` without checking that `x` is
   non-empty or that `n >= 1`. Unreachable from current R code; now guarded.
5. **`double` loop variables used as matrix indices** in `eucdistC_N1xN0()`.
6. **`weights_matrixC()` compared a raw `focal`** against a `treat` recoded to
   `0..g-1`. The identity for the 0/1 coding `matchit()` produces, so latent only;
   `focal` is now recoded the same way, as `nn_matchC_vec()` already did.

## Cleanups

- `internal.cpp` now includes `internal.h`, so declaration/definition drift is a
  compile error rather than a link error. Doing this immediately caught a duplicated
  default argument, which is exactly the class of thing it guards against. Defaults
  now live only in the header, and `<numeric>`/`<vector>`, previously arriving only
  by way of `Rcpp.h`, are included explicitly.
- Scalars are passed by value rather than as `const int&`/`const double&`/`const bool&`
  throughout `internal.{h,cpp}` and the small standalone files. The six exported
  `nn_matchC_*` signatures were left alone deliberately, to keep the diff in the
  matching code small.
- `using namespace Rcpp;` is gone from `internal.h` (it leaked into every translation
  unit); the header's own declarations are `Rcpp::`-qualified.
- The duplicated ~50-line tail of `find_control_vec()` and `find_control_mat()` is now
  one `take_closest()` helper. It keeps the `partial_sort()`/`sort()` split as it was:
  the two order ties differently, and collapsing them to one call would silently
  change which control lands in which column of `match.matrix` for tied distances.
- `subclass2mmC()` went from O(n * n1) to O(n + n1) by building a subclass-to-row
  lookup once instead of scanning every focal unit for every control unit. The lookup
  keeps the first matching row, which is what the old `break` did.
- `get_splitsC()` accumulated splits with `NumericVector::push_back()`, which
  reallocates and copies the whole vector each time; it now builds a `std::vector`.
- Sixteen `// [[Rcpp::interfaces(cpp)]]` comments (a file-scoped attribute applied
  per function) were generating `inst/include/MatchIt*.h`, a `registerCCallable`
  block in `RcppExports.cpp`, and a `methods::setLoadAction()` call in
  `RcppExports.R` — all inert, since the exported-signature table was empty and no
  package declares `LinkingTo: MatchIt`. Note that the `setLoadAction()` call
  referenced *methods*, which is not in `Imports`; that latent check problem is gone
  with it.
- `src/Makevars` and `src/Makevars.win` passed OpenMP, BLAS, LAPACK, and FLIBS flags
  with no `omp` pragma or BLAS call anywhere in `src/`. Removed, which also drops a
  pointless gfortran/OpenMP link from the shared object.
- The ten `// [[Rcpp::plugins(cpp11)]]` comments are honored only by `sourceCpp()`,
  never in a package build. Removed. The standard in force is R's default, which is
  at least C++14 for R >= 4.1 (the `DESCRIPTION` floor); everything written here is
  valid C++11 regardless.

## Warnings

`-Wall -Wextra -pedantic` over all of `src/`, before and after: **0** warnings from
the hand-written files in both cases. In the *generated* `RcppExports.cpp` the count
fell from 20 to 17, the three lost being the unused `validateSignature` and two
function-pointer casts that the removed interface required. The remaining 17
(`-Wcast-function-type-mismatch`, `-Wunused-parameter`) are inherent to Rcpp's
generated registration table.

## Not done

- Restructuring the nearest-neighbour search: passing `mm` plus a row index instead of
  `mm.row()` (an allocation per treated unit per group), the two-sided scan itself,
  and `find_control_mat()`'s prune against the running maximum rather than the
  ratio-th best. Out of scope by instruction, and each is a design change.
- Replacing `base::order()` for the primary sorts or for `heap_ord`. Radix sort is
  fast and its tie order is what the snapshots encode. Only the second `order()` call
  in four files was replaced, with the inverse permutation it provably computes.
- Restructuring `subclass_scootC()`. It relies on the invariant "`min(subtab) <= 0`
  implies some `subtab[s] == 0`", which holds; the search for that subclass now has a
  defensive `break` for `s == nsub`, so if the invariant is ever broken the function
  bails out of the while loop instead of assigning an out-of-range subclass and then
  reading `unique_sub` past its end.
