# Snapshot tests for method = "nearest".
# These pin the exact match.matrix, weights, and subclass to detect changes in the
# C++ matching algorithms and in the R code that feeds and post-processes them.
# Complement test-method_nearest.R, which tests structural correctness.
#
# Every test is deterministic: set.seed() is only load-bearing for
# m.order = "random", but is set throughout so that adding randomness to a code
# path cannot silently make a snapshot unreproducible.
#
# Between them these tests reach all six nn_matchC_*() entry points and cross each
# of exact, antiexact, distance calipers, covariate calipers, unit.id, discard,
# ratio > 1, variable ratio, and bounded reuse with every engine that can reach it.
# Two combinations are unreachable by construction and so are absent: a distance
# caliper with `distance` supplied as a matrix (there is no propensity score to
# apply it to), and variable ratio matching with `distance` supplied as a matrix
# (matchit() errors).

data("lalonde", package = "MatchIt")

# Fixed subset for distance matrix tests (50 treated + 50 control)
lalonde_sub <- lalonde[c(1:50, 186:235), ]

# Fixed sampling weights, for tests of s.weights
lalonde_sw <- seq(0.5, 2, length.out = nrow(lalonde))

# Three consecutive observations per unit ID, for tests of unit.id
lalonde_clust <- lalonde
lalonde_clust$clust <- rep(seq_len(ceiling(nrow(lalonde) / 3)),
                           each = 3L)[seq_len(nrow(lalonde))]

lalonde_sub_clust <- lalonde_sub
lalonde_sub_clust$clust <- rep(seq_len(ceiling(nrow(lalonde_sub) / 3)),
                               each = 3L)[seq_len(nrow(lalonde_sub))]

# ===== Baseline tests: one per C++ code path =====

test_that("baseline: PS vector, m.order='largest' (default)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest")
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

test_that("baseline: Mahalanobis, m.order='data'", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               mahvars = ~ age + educ + re74 + re75,
               m.order = "data")
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

test_that("baseline: distance matrix, m.order='data'", {
  set.seed(12345)
  d <- scaled_euclidean_dist(treat ~ age + educ + re74 + re75, data = lalonde_sub)
  m <- matchit(treat ~ age + educ + re74 + re75, data = lalonde_sub,
               distance = d, m.order = "data")
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

test_that("baseline: PS vector, m.order='closest'", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               m.order = "closest")
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

test_that("baseline: Mahalanobis, m.order='closest'", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               mahvars = ~ age + educ + re74 + re75,
               m.order = "closest")
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

test_that("baseline: distance matrix, m.order='closest'", {
  set.seed(12345)
  d <- scaled_euclidean_dist(treat ~ age + educ + re74 + re75, data = lalonde_sub)
  m <- matchit(treat ~ age + educ + re74 + re75, data = lalonde_sub,
               distance = d, m.order = "closest")
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

# ===== m.order variants =====

test_that("PS vector, m.order='data'", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               m.order = "data")
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

test_that("PS vector, m.order='random'", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               m.order = "random")
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

test_that("PS vector, m.order='farthest' (close=FALSE)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               m.order = "farthest")
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

test_that("Mahalanobis, m.order='random' (mahcovs path)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               mahvars = ~ age + educ + re74 + re75,
               m.order = "random")
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

test_that("Mahalanobis, m.order='farthest' (mahcovs_closest close=FALSE)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               mahvars = ~ age + educ + re74 + re75,
               m.order = "farthest")
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

test_that("distance matrix, m.order='random'", {
  set.seed(12345)
  d <- scaled_euclidean_dist(treat ~ age + educ + re74 + re75, data = lalonde_sub)
  m <- matchit(treat ~ age + educ + re74 + re75, data = lalonde_sub,
               distance = d, m.order = "random")
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

test_that("distance matrix, m.order='farthest' (distmat_closest close=FALSE)", {
  set.seed(12345)
  d <- scaled_euclidean_dist(treat ~ age + educ + re74 + re75, data = lalonde_sub)
  m <- matchit(treat ~ age + educ + re74 + re75, data = lalonde_sub,
               distance = d, m.order = "farthest")
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

# ===== ratio + replacement =====

test_that("ratio=3, replace=FALSE (pool depletion)", {
  set.seed(12345)
  expect_wrn(
    m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               ratio = 3, replace = FALSE),
    "Not all treated units will get 3 matches"
  )
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 3L)
  expect_matchit_snapshot(m)
})

test_that("ratio=3, replace=TRUE (no pool depletion)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               ratio = 3, replace = TRUE)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = FALSE, ratio = 3L, replace = TRUE)
  expect_matchit_snapshot(m)
})

# ===== ratio + caliper =====

test_that("ratio=2, positive caliper (pool restriction)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               ratio = 2, caliper = 0.1, std.caliper = FALSE)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 2L)
  expect_matchit_snapshot(m)
})

test_that("ratio=2, negative caliper (anti-caliper)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               ratio = 2, caliper = -0.05, std.caliper = FALSE)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 2L)
  expect_matchit_snapshot(m)
})

# ===== caliper + replacement =====

test_that("caliper + replace=TRUE + ratio=2 (reuse within caliper)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               caliper = 0.1, std.caliper = FALSE,
               replace = TRUE, ratio = 2)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = FALSE, ratio = 2L, replace = TRUE)
  expect_matchit_snapshot(m)
})

# ===== m.order + caliper + no replacement =====

test_that("m.order='largest' + caliper + replace=FALSE", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               m.order = "largest", caliper = 0.1, std.caliper = FALSE,
               replace = FALSE)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

test_that("m.order='smallest' + caliper + replace=FALSE", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               m.order = "smallest", caliper = 0.1, std.caliper = FALSE,
               replace = FALSE)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

# ===== exact + other constraints =====

test_that("exact + ratio=2 (within-stratum depletion)", {
  set.seed(12345)
  expect_warning(
    expect_warning(
      m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
                   data = lalonde, method = "nearest",
                   exact = ~ race, ratio = 2),
      "Fewer control units than treated units"
    ),
    "Not all treated units will get 2 matches"
  )
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 2L)
  expect_matchit_snapshot(m)
})

test_that("exact + caliper (double constraint)", {
  set.seed(12345)
  expect_wrn(
    m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               exact = ~ race, caliper = 0.2, std.caliper = FALSE),
    "Fewer control units than treated units"
  )
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

test_that("exact + antiexact (inclusion + exclusion)", {
  set.seed(12345)
  expect_wrn(
    m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               exact = ~ race, antiexact = ~ married),
    "Fewer control units than treated units"
  )
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

test_that("exact + replace=TRUE + ratio=2 (reuse within strata)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               exact = ~ race, replace = TRUE, ratio = 2)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = FALSE, ratio = 2L, replace = TRUE)
  expect_matchit_snapshot(m)
})

# ===== mahvars + caliper =====

test_that("mahvars + distance caliper (Mahalanobis match with PS caliper)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               mahvars = ~ age + educ + re74 + re75,
               caliper = 0.2, std.caliper = FALSE)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

test_that("mahvars + exact + m.order='closest' (three-way)", {
  set.seed(12345)
  expect_wrn(
    m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               mahvars = ~ age + educ + re74 + re75,
               exact = ~ race, m.order = "closest"),
    "Fewer control units than treated units"
  )
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

# ===== mahvars + replace =====

test_that("mahvars + replace=TRUE (mahcovs path, reuse allowed)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               mahvars = ~ age + educ + re74 + re75,
               replace = TRUE)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = FALSE, ratio = 1L, replace = TRUE)
  expect_matchit_snapshot(m)
})

test_that("mahvars + replace=TRUE + ratio=2 (mahcovs, reuse, multi-match)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               mahvars = ~ age + educ + re74 + re75,
               replace = TRUE, ratio = 2)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = FALSE, ratio = 2L, replace = TRUE)
  expect_matchit_snapshot(m)
})

# ===== mahvars + antiexact =====

test_that("mahvars + antiexact (mahcovs path with antiexact)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               mahvars = ~ age + educ + re74 + re75,
               antiexact = ~ married)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

# ===== reuse.max =====

test_that("reuse.max=3 + ratio=2 (bounded replacement)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               ratio = 2, reuse.max = 3)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = FALSE, ratio = 2L,
                      replace = structure(TRUE, reuse.max = 3))
  expect_matchit_snapshot(m)
})

test_that("reuse.max=2 + ratio=2 + caliper (bounded replacement + caliper)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               ratio = 2, reuse.max = 2,
               caliper = 0.2, std.caliper = FALSE)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = FALSE, ratio = 2L,
                      replace = structure(TRUE, reuse.max = 2))
  expect_matchit_snapshot(m)
})

test_that("mahvars + reuse.max=3 + ratio=2 (mahcovs bounded replacement)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               mahvars = ~ age + educ + re74 + re75,
               ratio = 2, reuse.max = 3)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = FALSE, ratio = 2L,
                      replace = structure(TRUE, reuse.max = 3))
  expect_matchit_snapshot(m)
})

test_that("mahvars + reuse.max=2 + m.order='closest' (mahcovs_closest bounded)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               mahvars = ~ age + educ + re74 + re75,
               ratio = 2, reuse.max = 2, m.order = "closest")
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = FALSE, ratio = 2L,
                      replace = structure(TRUE, reuse.max = 2))
  expect_matchit_snapshot(m)
})

test_that("distmat + reuse.max=3 + ratio=2 (distmat bounded replacement)", {
  set.seed(12345)
  d <- scaled_euclidean_dist(treat ~ age + educ + re74 + re75, data = lalonde_sub)
  m <- matchit(treat ~ age + educ + re74 + re75, data = lalonde_sub,
               distance = d, ratio = 2, reuse.max = 3)
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = TRUE,
                      expect_subclass = FALSE, ratio = 2L,
                      replace = structure(TRUE, reuse.max = 3))
  expect_matchit_snapshot(m)
})

# ===== variable ratio =====

test_that("variable ratio + caliper (min/max with restriction)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               ratio = 2, min.controls = 1, max.controls = 4,
               caliper = 0.2, std.caliper = FALSE)
  ratio_attr <- structure(2L, min.controls = 1, max.controls = 4)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = ratio_attr)
  expect_matchit_snapshot(m)
})

test_that("variable ratio + exact (within strata)", {
  set.seed(12345)
  expect_warning(
    expect_warning(
      m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
                   data = lalonde, method = "nearest",
                   ratio = 2, min.controls = 1, max.controls = 4,
                   exact = ~ race),
      "Fewer control units than treated units"
    ),
    "Not enough control units"
  )
  ratio_attr <- structure(2L, min.controls = 1, max.controls = 4)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = ratio_attr)
  expect_matchit_snapshot(m)
})

test_that("variable ratio baseline (min/max.controls, PS vector)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               ratio = 2, min.controls = 1, max.controls = 4)
  ratio_attr <- structure(2L, min.controls = 1, max.controls = 4)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = ratio_attr)
  expect_matchit_snapshot(m)
})

# ===== estimand = "ATC" =====

test_that("estimand='ATC' baseline (PS vector, flipped focal)", {
  set.seed(12345)
  expect_wrn(
    m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               estimand = "ATC"),
    "Fewer treated units than control units"
  )
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

test_that("estimand='ATC' + mahvars (mahcovs path, ATC focal)", {
  set.seed(12345)
  expect_wrn(
    m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               mahvars = ~ age + educ + re74 + re75,
               estimand = "ATC"),
    "Fewer treated units than control units"
  )
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

test_that("estimand='ATC' + ratio=2 + exact (flipped focal)", {
  set.seed(12345)
  expect_warning(
    expect_warning(
      m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
                   data = lalonde, method = "nearest",
                   estimand = "ATC", ratio = 2, exact = ~ race),
      "Fewer treated units than control units"
    ),
    "Not all control units will get 2 matches"
  )
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 2L)
  expect_matchit_snapshot(m)
})

# ===== multiple calipers =====

test_that("covariate caliper + distance caliper (both simultaneously)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               caliper = c(.1, age = 2), std.caliper = FALSE)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

test_that("covariate caliper + antiexact + m.order='closest'", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               caliper = c(age = 5), std.caliper = FALSE,
               antiexact = ~ married, m.order = "closest")
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

# ===== std.caliper = TRUE =====

test_that("standardized caliper (std.caliper=TRUE, PS vector)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               caliper = 0.25, std.caliper = TRUE)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

test_that("standardized covariate caliper (std.caliper=TRUE on age)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               caliper = c(age = 1), std.caliper = TRUE)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

# ===== antiexact without exact =====

test_that("antiexact alone, PS vector", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               antiexact = ~ married)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

test_that("antiexact + m.order='closest' (vec_closest path)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               antiexact = ~ married, m.order = "closest")
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

# ===== discard =====

test_that("discard logical vector, PS vector path", {
  set.seed(12345)
  dis <- lalonde$re74 > 15000
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               discard = dis)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

test_that("discard + ratio=2 + m.order='closest' (vec_closest path with discards)", {
  set.seed(12345)
  dis <- lalonde$re74 > 10000
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               discard = dis, ratio = 2, m.order = "closest")
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 2L, replace = FALSE)
  expect_matchit_snapshot(m)
})

test_that("discard + mahvars (mahcovs path with discards)", {
  set.seed(12345)
  dis <- lalonde$re74 > 15000
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               mahvars = ~ age + educ + re74 + re75,
               discard = dis)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

# ===== unit.id =====

#`unit.id = ~ age` is a high-contention case: the 429 control units share only ~40
#IDs, so most treated units go unmatched. Realistic multi-observation clusters are
#covered by the `lalonde_clust` tests below.

test_that("unit.id with replacement=FALSE (clustered units, PS vector)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               unit.id = ~ age)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

test_that("unit.id + m.order='closest' (vec_closest path with unit.id)", {
  set.seed(12345)
  expect_wrn(
    m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               unit.id = ~ age, m.order = "closest", ratio = 2),
    "Not all treated units will get 2 matches"
  )
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 2L)
  expect_matchit_snapshot(m)
})

# ===== distance matrix + constraints =====

test_that("distmat + exact + ratio=2 (stratum loop + multi-match)", {
  set.seed(12345)
  d <- scaled_euclidean_dist(treat ~ age + educ + re74 + re75, data = lalonde_sub)
  expect_warning(
    expect_warning(
      m <- matchit(treat ~ age + educ + re74 + re75, data = lalonde_sub,
                   distance = d, exact = ~ race, ratio = 2),
      "Fewer control units than treated units"
    ),
    "Not all treated units will get 2 matches"
  )
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 2L)
  expect_matchit_snapshot(m)
})

#`m.order` is deliberately not set here: with `replace = TRUE` the matching is
#order-invariant and `matchit2nearest()` forces `m.order = "data"`, so the
#`distmat_closest` algorithm cannot be reached this way. It is covered instead by
#"distmat + reuse.max=2 + m.order='closest'" below.
test_that("distmat + replace=TRUE + ratio=2", {
  set.seed(12345)
  d <- scaled_euclidean_dist(treat ~ age + educ + re74 + re75, data = lalonde_sub)
  m <- matchit(treat ~ age + educ + re74 + re75, data = lalonde_sub,
               distance = d, replace = TRUE, ratio = 2)
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = TRUE,
                      expect_subclass = FALSE, ratio = 2L, replace = TRUE)
  expect_matchit_snapshot(m)
})

test_that("distmat + caliper (distmat path with caliper constraint)", {
  set.seed(12345)
  d <- scaled_euclidean_dist(treat ~ age + educ + re74 + re75, data = lalonde_sub)
  m <- matchit(treat ~ age + educ + re74 + re75, data = lalonde_sub,
               distance = d, caliper = c(age = 3), std.caliper = FALSE)
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

test_that("distmat + antiexact (distmat path with antiexact constraint)", {
  set.seed(12345)
  d <- scaled_euclidean_dist(treat ~ age + educ + re74 + re75, data = lalonde_sub)
  m <- matchit(treat ~ age + educ + re74 + re75, data = lalonde_sub,
               distance = d, antiexact = ~ married)
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

# ===== full Mahalanobis distance =====

test_that("full Mahalanobis (distance='mahalanobis')", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + re74 + re75,
               data = lalonde, method = "nearest",
               distance = "mahalanobis")
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

test_that("full Mahalanobis + m.order='closest'", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + re74 + re75,
               data = lalonde, method = "nearest",
               distance = "mahalanobis", m.order = "closest")
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

# ===== other full-distance transforms =====

test_that("full robust Mahalanobis (distance='robust_mahalanobis')", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + re74 + re75,
               data = lalonde, method = "nearest",
               distance = "robust_mahalanobis")
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

test_that("full scaled Euclidean (distance='scaled_euclidean')", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + re74 + re75,
               data = lalonde, method = "nearest",
               distance = "scaled_euclidean")
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

test_that("full Euclidean (distance='euclidean')", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + re74 + re75,
               data = lalonde, method = "nearest",
               distance = "euclidean")
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

#A single Mahalanobis covariate is collapsed to a distance vector and matched by
#`nn_matchC_vec()` rather than `nn_matchC_mahcovs()`.
test_that("full Mahalanobis, single covariate (collapses to vec path)", {
  set.seed(12345)
  m <- matchit(treat ~ age,
               data = lalonde, method = "nearest",
               distance = "mahalanobis")
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

# ===== link functions and propensity scores =====

test_that("link='probit'", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               link = "probit")
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

test_that("link='linear.logit' + caliper (caliper on the linear predictor)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               link = "linear.logit", caliper = 0.25)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

#The estimated propensity score is an input to almost every other test in this
#file, but nothing else pins its values, so a change in how `glm()` is called would
#surface only as an unexplained change in `match.matrix`.
test_that("estimated propensity scores are stable across links", {
  skip_on_cran()

  for (link in c("logit", "probit", "cloglog", "linear.logit")) {
    m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
                 data = lalonde, method = "nearest",
                 link = link)
    expect_snapshot_value(unname(round(m$distance, 8L)), style = "json2")
  }
})

# ===== sampling weights =====

test_that("s.weights + PS (weighted propensity score estimation)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               s.weights = lalonde_sw)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

test_that("s.weights + mahvars (weighted Mahalanobis scaling)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               mahvars = ~ age + educ + re74 + re75,
               s.weights = lalonde_sw)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

test_that("s.weights + full Mahalanobis", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + re74 + re75,
               data = lalonde, method = "nearest",
               distance = "mahalanobis", s.weights = lalonde_sw)
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

# ===== discard and reestimate =====

test_that("discard='both' (common support on the PS)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               discard = "both")
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

test_that("discard='both' + reestimate=TRUE (PS refit after discarding)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               discard = "both", reestimate = TRUE)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

test_that("discard='control'", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               discard = "control")
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

# ===== mahvars: factors and covariate calipers =====

test_that("mahvars with a factor covariate", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               mahvars = ~ age + educ + race)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

#A positive covariate caliper combined with `mahvars` is implemented by splitting
#the caliper variable into exact-matching strata via `get_splitsC()`.
test_that("mahvars + covariate caliper (caliper splitting)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               mahvars = ~ age + educ + re74 + re75,
               caliper = c(age = 2), std.caliper = FALSE)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

test_that("mahvars + exact + covariate caliper (split strata crossed with exact)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               mahvars = ~ age + educ + re74 + re75,
               exact = ~ married,
               caliper = c(age = 2), std.caliper = FALSE)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

#Negative calipers cannot be turned into strata, so this takes the branch that
#skips the split and applies the caliper directly.
test_that("mahvars + negative covariate caliper (no caliper splitting)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               mahvars = ~ age + educ + re74 + re75,
               caliper = c(age = -2), std.caliper = FALSE)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

test_that("mahvars + m.order='smallest'", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               mahvars = ~ age + educ + re74 + re75,
               m.order = "smallest")
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

# ===== estimand = "ATC" with non-PS distances =====

#With `estimand = "ATC"` a supplied distance matrix is transposed before matching.
test_that("estimand='ATC' + distance matrix (transposed)", {
  set.seed(12345)
  d <- scaled_euclidean_dist(treat ~ age + educ + re74 + re75, data = lalonde_sub)
  m <- matchit(treat ~ age + educ + re74 + re75, data = lalonde_sub,
               distance = d, estimand = "ATC")
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

test_that("estimand='ATC' + full Mahalanobis", {
  set.seed(12345)
  expect_wrn(
    m <- matchit(treat ~ age + educ + re74 + re75,
               data = lalonde, method = "nearest",
               distance = "mahalanobis", estimand = "ATC"),
    "Fewer treated units than control units"
  )
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

# ===== unit.id with multi-observation clusters =====

test_that("unit.id (multi-obs clusters) + ratio=2", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde_clust, method = "nearest",
               unit.id = ~ clust, ratio = 2)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 2L)
  expect_matchit_snapshot(m)
})

test_that("unit.id + exact (clusters within strata)", {
  set.seed(12345)
  expect_wrn(
    m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde_clust, method = "nearest",
               unit.id = ~ clust, exact = ~ race),
    "Fewer control unit IDs than treated units"
  )
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

test_that("unit.id + reuse.max=2 + ratio=2 (bounded reuse of clusters)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde_clust, method = "nearest",
               unit.id = ~ clust, reuse.max = 2, ratio = 2)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = FALSE, ratio = 2L,
                      replace = structure(TRUE, reuse.max = 2))
  expect_matchit_snapshot(m)
})

# ===== additional distance matrix specifications =====

#A full n x n matrix is subset to the treated-by-control block internally.
test_that("distance supplied as a full n x n matrix", {
  set.seed(12345)
  X <- scale(lalonde_sub[c("age", "educ", "re74", "re75")])
  m <- matchit(treat ~ age + educ + re74 + re75, data = lalonde_sub,
               distance = as.matrix(dist(X)))
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

#Reaches `nn_matchC_distmat_closest()` with reuse, which `replace = TRUE` cannot.
test_that("distmat + reuse.max=2 + m.order='closest' (distmat_closest bounded)", {
  set.seed(12345)
  d <- scaled_euclidean_dist(treat ~ age + educ + re74 + re75, data = lalonde_sub)
  m <- matchit(treat ~ age + educ + re74 + re75, data = lalonde_sub,
               distance = d, m.order = "closest", ratio = 2, reuse.max = 2)
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = TRUE,
                      expect_subclass = FALSE, ratio = 2L,
                      replace = structure(TRUE, reuse.max = 2))
  expect_matchit_snapshot(m)
})

# ===== multi-constraint combinations =====

test_that("exact + antiexact + caliper + ratio=2 + m.order='closest'", {
  set.seed(12345)
  expect_warning(
    expect_warning(
      m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
                   data = lalonde, method = "nearest",
                   exact = ~ race, antiexact = ~ married,
                   caliper = 0.2, std.caliper = FALSE,
                   ratio = 2, m.order = "closest"),
      "Fewer control units than treated units"
    ),
    "Not all treated units will get 2 matches"
  )
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 2L)
  expect_matchit_snapshot(m)
})

test_that("multi-variable exact + antiexact", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               exact = ~ married + nodegree, antiexact = ~ race)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

test_that("variable ratio + m.order='closest'", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               ratio = 2, min.controls = 1, max.controls = 4,
               m.order = "closest")
  ratio_attr <- structure(2L, min.controls = 1, max.controls = 4)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = ratio_attr)
  expect_matchit_snapshot(m)
})

# ===== invariants the snapshots above rely on =====

#Several tests above pass `replace = TRUE`; this documents why they cannot also be
#used to test `m.order`.
test_that("m.order is ignored when matching with replacement", {
  mm <- lapply(c("data", "closest", "farthest", "random", "largest", "smallest"),
               function(mo) {
                 set.seed(12345)
                 m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
                              data = lalonde, method = "nearest",
                              replace = TRUE, ratio = 2, m.order = mo)
                 m$match.matrix
               })

  for (i in seq_along(mm)[-1L]) {
    expect_identical(mm[[i]], mm[[1L]])
  }
})

#Explains why "baseline: Mahalanobis, m.order='data'" and "full Mahalanobis
#(distance='mahalanobis')" record identical snapshots despite taking different
#branches through `matchit2nearest()`.
test_that("mahvars and distance='mahalanobis' agree on the same covariates", {
  set.seed(12345)
  m1 <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
                data = lalonde, method = "nearest",
                mahvars = ~ age + educ + re74 + re75, m.order = "data")
  set.seed(12345)
  m2 <- matchit(treat ~ age + educ + re74 + re75,
                data = lalonde, method = "nearest",
                distance = "mahalanobis")

  expect_identical(m1$match.matrix, m2$match.matrix)
  expect_identical(m1$weights, m2$weights)
})

# ===== constraints under m.order = "closest" =====

# The `_closest` algorithms re-find and re-rank matches as controls are used up, so
# they apply each constraint by a different mechanism than their counterparts above.
# These pin every constraint against the mahcovs_closest and distmat_closest engines.

test_that("mahvars + m.order='closest' + distance caliper", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               mahvars = ~ age + educ + re74 + re75,
               m.order = "closest", caliper = 0.2, std.caliper = FALSE)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

test_that("mahvars + m.order='closest' + covariate caliper", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               mahvars = ~ age + educ + re74 + re75,
               m.order = "closest", caliper = c(age = 2), std.caliper = FALSE)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

test_that("mahvars + m.order='closest' + antiexact", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               mahvars = ~ age + educ + re74 + re75,
               m.order = "closest", antiexact = ~ married)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

test_that("mahvars + m.order='closest' + discard", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               mahvars = ~ age + educ + re74 + re75,
               m.order = "closest", discard = lalonde$re74 > 15000)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

test_that("mahvars + unit.id + ratio=2", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde_clust, method = "nearest",
               mahvars = ~ age + educ + re74 + re75,
               unit.id = ~ clust, ratio = 2)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 2L)
  expect_matchit_snapshot(m)
})

test_that("mahvars + m.order='closest' + unit.id + ratio=2", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde_clust, method = "nearest",
               mahvars = ~ age + educ + re74 + re75,
               m.order = "closest", unit.id = ~ clust, ratio = 2)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 2L)
  expect_matchit_snapshot(m)
})

# ===== variable ratio with mahvars =====

#Variable ratio matching needs a propensity score, so it is available with
#`mahvars` even though it is not with `distance` supplied as a matrix.

test_that("mahvars + variable ratio", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               mahvars = ~ age + educ + re74 + re75,
               ratio = 2, min.controls = 1, max.controls = 4)
  ratio_attr <- structure(2L, min.controls = 1, max.controls = 4)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = ratio_attr)
  expect_matchit_snapshot(m)
})

test_that("mahvars + variable ratio + m.order='closest'", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               mahvars = ~ age + educ + re74 + re75,
               ratio = 2, min.controls = 1, max.controls = 4,
               m.order = "closest")
  ratio_attr <- structure(2L, min.controls = 1, max.controls = 4)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = ratio_attr)
  expect_matchit_snapshot(m)
})

# ===== distance matrix under constraints =====

test_that("distmat + m.order='closest' + covariate caliper", {
  set.seed(12345)
  d <- scaled_euclidean_dist(treat ~ age + educ + re74 + re75, data = lalonde_sub)
  m <- matchit(treat ~ age + educ + re74 + re75, data = lalonde_sub,
               distance = d, m.order = "closest",
               caliper = c(age = 3), std.caliper = FALSE)
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

test_that("distmat + m.order='closest' + antiexact", {
  set.seed(12345)
  d <- scaled_euclidean_dist(treat ~ age + educ + re74 + re75, data = lalonde_sub)
  m <- matchit(treat ~ age + educ + re74 + re75, data = lalonde_sub,
               distance = d, m.order = "closest", antiexact = ~ married)
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

test_that("distmat + m.order='closest' + exact (stratum loop)", {
  set.seed(12345)
  d <- scaled_euclidean_dist(treat ~ age + educ + re74 + re75, data = lalonde_sub)
  expect_wrn(
    m <- matchit(treat ~ age + educ + re74 + re75, data = lalonde_sub,
               distance = d, m.order = "closest", exact = ~ race),
    "Fewer control units than treated units"
  )
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

test_that("distmat + unit.id", {
  set.seed(12345)
  d <- scaled_euclidean_dist(treat ~ age + educ + re74 + re75, data = lalonde_sub)
  m <- matchit(treat ~ age + educ + re74 + re75, data = lalonde_sub_clust,
               distance = d, unit.id = ~ clust)
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

test_that("distmat + m.order='closest' + unit.id", {
  set.seed(12345)
  d <- scaled_euclidean_dist(treat ~ age + educ + re74 + re75, data = lalonde_sub)
  m <- matchit(treat ~ age + educ + re74 + re75, data = lalonde_sub_clust,
               distance = d, unit.id = ~ clust, m.order = "closest")
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

test_that("distmat + discard", {
  set.seed(12345)
  d <- scaled_euclidean_dist(treat ~ age + educ + re74 + re75, data = lalonde_sub)
  m <- matchit(treat ~ age + educ + re74 + re75, data = lalonde_sub,
               distance = d, discard = lalonde_sub$re74 > 15000)
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})

test_that("distmat + m.order='closest' + discard", {
  set.seed(12345)
  d <- scaled_euclidean_dist(treat ~ age + educ + re74 + re75, data = lalonde_sub)
  m <- matchit(treat ~ age + educ + re74 + re75, data = lalonde_sub,
               distance = d, discard = lalonde_sub$re74 > 15000,
               m.order = "closest")
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_matchit_snapshot(m)
})
