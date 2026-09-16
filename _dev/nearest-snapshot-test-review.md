# Review of the `method = "nearest"` snapshot tests

Audit and expansion of `tests/testthat/test-nearest_snapshots.R`, done to make it a usable regression baseline ahead of a major update. Everything below was verified by running code, not by reading it.

## Summary

|   | before | after |
|----|----|----|
| tests | 57 | 102 |
| snapshot values per test | 1 (`match.matrix`) | 3 (`match.matrix`, `weights`, `subclass`) |
| `_snaps/nearest_snapshots.md` | 269 KB | 928 KB |
| `nn_matchC_*()` entry points reached | 6 / 6 | 6 / 6 |
| (engine, constraint) pairs with zero coverage | 20 | 4, all structurally impossible |

All 57 original `match.matrix` snapshots still reproduce **byte for byte** — the only line removed from the snapshot file is one renamed test header. Nothing in the existing baseline was invalidated. Full suite: 108 tests, 0 failures, \~22 s.

## What the original tests got right

- All six C++ entry points (`nn_matchC_vec`, `_vec_closest`, `_mahcovs`, `_mahcovs_closest`, `_distmat`, `_distmat_closest`) were genuinely reached. Verified by wrapping each binding in the namespace and logging calls.
- No two tests produced identical snapshots except one intentional pair (below), so the specs were not accidentally redundant.
- `set.seed()` before each call is only load-bearing for `m.order = "random"`, but is harmless and worth keeping.
- Wrapping the calls that warn in `expect_warning()` pins the warning behavior too, and the short regexes avoid breaking on *cli*'s line wrapping.

## Problems found

### 1. `weights` and `subclass` were not pinned at all

`match.matrix` was the only thing recorded. When `reuse.max > 1` (i.e. `replace = TRUE` or a finite `reuse.max`), `subclass` is `NULL` and the weights come from `get_weights_from_mm()` — and the structural check in `expect_good_matchit()` only cross-checks the two weight routes when *both* `subclass` and `match.matrix` exist. So for every matching-with-replacement test, the weights — the output users actually consume — were unconstrained. One such case had 16 distinct weight values, none of them pinned.

Fixed by a new `expect_matchit_snapshot()` helper in `helpers.R` that records `match.matrix`, then `weights` (rounded to 8 digits, unnamed), then `subclass` (as unnamed integer codes). `match.matrix` stays first so the existing snapshots keep matching.

### 2. Two tests claimed to exercise `m.order` under `replace = TRUE`, which is a no-op

`matchit2nearest()` forces `m.order <- "data"` when `reuse.max >= n1`. Confirmed empirically: with `replace = TRUE`, all six `m.order` values give an identical `match.matrix`. So `"distmat + m.order='closest' + replace=TRUE + ratio=2"` was running `nn_matchC_distmat`, not `nn_matchC_distmat_closest` — a misleading name over a redundant spec.

Renamed to `"distmat + replace=TRUE + ratio=2"` with the `m.order` argument dropped (dropping it changed nothing, which the unchanged snapshot content confirms), the intended coverage restored by a new `reuse.max = 2` variant, and the invariant itself pinned by a new test, `"m.order is ignored when matching with replacement"`.

### 3. The `_closest` engines had almost no constraint coverage

Instrumenting `nn_matchC_dispatch()` across the whole file showed that `mahcovs_closest` and `distmat_closest` were reached only in unconstrained form — zero tests crossed either with `exact`, `antiexact`, a caliper, `unit.id`, `discard`, or variable ratio. These are the engines where a refactor is most likely to break something quietly, since they re-find and re-rank matches as controls are used up and so apply every constraint by a different mechanism than their non-`closest` counterparts. 14 tests added.

### 4. The propensity score itself was never pinned

The estimated PS is an input to most tests in the file, but no expectation recorded its values, so a change in how `glm()` is called would have surfaced only as an unexplained `match.matrix` diff. Demonstrated: shifting every PS by `+1e-6` changes no `match.matrix`, no weights, and no subclass anywhere in the suite (a uniform shift preserves all pairwise differences and the ordering) — it was completely invisible. New test `"estimated propensity scores are stable across links"` pins the PS vector for `logit`, `probit`, `cloglog`, and `linear.logit`.

### 5. Two `unit.id` tests were near-degenerate

`unit.id = ~ age` leaves the 429 control units sharing \~40 IDs, so 167 of 185 treated units go unmatched and the snapshot is mostly `NA`. Kept (it is a legitimate high-contention case, now flagged as such in a comment) and supplemented with a `lalonde_clust` fixture giving three observations per ID.

### 6. Untested arguments

Reachable specifications with no coverage at all, now added: `s.weights` (which changes both the PS fit and the Mahalanobis scaling), `discard` as a string plus `reestimate = TRUE` (only a logical vector was tested), `link` other than the default, `distance = "robust_mahalanobis"` / `"scaled_euclidean"` / `"euclidean"`, `mahvars` containing a factor, `mahvars` with a covariate caliper (the `get_splitsC()` caliper-splitting path, which only the `distmat` tests reached), `mahvars` with a negative covariate caliper (the branch that skips splitting), `estimand = "ATC"` with a supplied distance matrix (the transpose branch), variable ratio with `mahvars`, a distance matrix supplied as full *n* × *n* rather than *n*₁ × *n*₀, and full Mahalanobis on a single covariate (which collapses to the `vec` engine).

### 7. One intentional duplicate

`"baseline: Mahalanobis, m.order='data'"` and `"full Mahalanobis (distance='mahalanobis')"` record identical snapshots. This is correct — they match on the same Mahalanobis distance by different routes through `matchit2nearest()` — so both were kept, and the equivalence is now asserted directly by `"mahvars and distance='mahalanobis' agree on the same covariates"` rather than left implicit in two identical blobs.

## Are the snapshots actually informative?

Mutation-tested by patching package source, running the suite, and reverting:

| mutation                                | failing tests                     |
|-----------------------------------------|-----------------------------------|
| baseline                                | 0 / 102                           |
| `weights_matrixC()` result × (1 + 1e-6) | 86 / 102                          |
| `mm2subclass()` levels reversed         | 86 / 102                          |
| estimated PS + 1e-6                     | 1 / 102 — exactly the new PS test |
| estimated PS × (1 + 1e-12)              | 0 / 102                           |

The last row matters as much as the others: the suite does **not** trip on last-bit numerical noise, so it should not produce false alarms from a different BLAS or platform. Rounding weights and propensity scores to 8 digits before recording is what buys that.

## Remaining gaps

The four uncovered (engine, constraint) pairs are unreachable by construction: a distance caliper with `distance` supplied as a matrix (there is no propensity score to apply it to) and variable ratio matching with `distance` supplied as a matrix (`matchit()` errors). Both are noted in the header comment of the test file.

**The significant one is out of scope of this file:** `method = "nearest"` is the only method with any result pinning. `optimal`, `full`, `quick`, `cem`, `exact`, `subclass`, `cardinality`, and `genetic` have structural tests at best and no pinned results, as do `summary.matchit()` and `match_data()`. If the goal is to retain results across a major update, those need the same treatment.

## Repository cleanup

- Removed the leftover `worktree-rippling-dreaming-kazoo` worktree and branch (its commit, `9576da7`, was already an ancestor of `master`; no unique work, no uncommitted files).
- Removed two stray gitlinks committed by mistake in `e352e52`: `.claude/worktrees/nearest-snapshot-tests` and `.claude/worktrees/rippling-dreaming-kazoo` were recorded as mode 160000 submodule entries with no `.gitmodules`, which breaks a fresh clone.
- Added `.claude/worktrees/` to `.gitignore` so it cannot recur.
- Left the `testing` branch alone: it is three commits from 2023 that diverged at `ca8a039`, has a counterpart at `origin/testing`, and is not agent debris.
