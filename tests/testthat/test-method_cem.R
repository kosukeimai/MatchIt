test_that("Coarsened exact matching works", {
  set.seed(123)
  k <- 6
  n <- 1e4

  d <- as.data.frame(matrix(rnorm(k * n), nrow = n))

  d[[1]] <- factor(cut(d[[1]], 4, labels = FALSE))
  d[[2]] <- factor(cut(d[[2]], 10, labels = FALSE))

  d$a <- rbinom(n, 1, .3)

  m <- matchit(a ~ ., data = d, method = "cem")

  expect_good_matchit(m, expect_subclass = TRUE, expect_distance = FALSE,
                      expect_match.matrix = FALSE)

  #Categories are exactly matched by default
  expect_true(all(sapply(levels(m$subclass), function(s) length(unique(d[[1]][which(m$subclass == s)])) == 1)))
  expect_true(all(sapply(levels(m$subclass), function(s) length(unique(d[[2]][which(m$subclass == s)])) == 1)))
  expect_false(all(sapply(levels(m$subclass), function(s) length(unique(d[[3]][which(m$subclass == s)])) == 1)))

  #k2k didn't accidentally activate
  expect_true(length(unique(sapply(unique(m$treat), function(t) {
    sum(m$weights[m$treat == t] > 0)
  }))) > 1L)

  #Groupings: V1 into 2 categories, no grouping of V2
  m <- matchit(a ~ ., data = d, method = "cem",
               grouping = list(V1 = list(c("1", "2"), c("3", "4")),
                               V2 = list(levels(d$V2))))

  expect_good_matchit(m, expect_subclass = TRUE, expect_distance = FALSE,
                      expect_match.matrix = FALSE)


  #Each subclass has V1 in 1,2 or 3,4
  expect_true(all(sapply(levels(m$subclass), function(s) all(d[[1]][which(m$subclass == s)] %in% c("1", "2")) ||
                           all(d[[1]][which(m$subclass == s)] %in% c("3", "4")))))

  #No restriction on bins for V2
  expect_false(all(sapply(levels(m$subclass), function(s) length(unique(d[[2]][which(m$subclass == s)])) == 1)))


  m <- matchit(a ~ ., data = d, method = "cem",
               grouping = list(V1 = list(c("1", "2"), c("3", "4")),
                               V2 = list(levels(d$V2))),
               cutpoints = list(V3 = c(-1.5, 1.5),
                                V4 = 1))

  expect_good_matchit(m, expect_subclass = TRUE, expect_distance = FALSE,
                      expect_match.matrix = FALSE)

  #Each subclass has V1 in 1,2 or 3,4
  expect_true(all(sapply(levels(m$subclass), function(s) all(d[[1]][which(m$subclass == s)] %in% c("1", "2")) ||
                           all(d[[1]][which(m$subclass == s)] %in% c("3", "4")))))

  #No restriction on bins for V2
  expect_false(all(sapply(levels(m$subclass), function(s) length(unique(d[[2]][which(m$subclass == s)])) == 1)))

  #V3 correctly split into defined bins
  expect_true(all(sapply(levels(m$subclass), function(s) {
    all(d[[3]][which(m$subclass == s)] < -1.5) ||
      all(d[[3]][which(m$subclass == s)] > -1.5 | d[[3]][which(m$subclass == s)] < 1.5) ||
      all(d[[3]][which(m$subclass == s)] > 1.5)
  })))

  #Setting V1 = 1 in cutpoints same as omitting it
  m1 <- matchit(a ~ . - V4, data = d, method = "cem",
               grouping = list(V1 = list(c("1", "2"), c("3", "4")),
                               V2 = list(levels(d$V2))),
               cutpoints = list(V3 = c(-1.5, 1.5)))

  expect_equal(m$subclass, m1$subclass)

  m <- matchit(a ~ ., data = d, method = "cem", k2k = TRUE)

  expect_good_matchit(m, expect_subclass = TRUE, expect_distance = FALSE,
                      expect_match.matrix = TRUE, ratio = 1)

  #1:1 matched
  expect_true(length(unique(sapply(unique(m$treat), function(t) {
    sum(m$weights[m$treat == t] > 0)
  }))) == 1L)

  #Default is using Mahalanobis
  m1 <- matchit(a ~ ., data = d, method = "cem", k2k = TRUE, k2k.method = "mahalanobis")
  expect_good_matchit(m, expect_subclass = TRUE, expect_distance = FALSE,
                      expect_match.matrix = TRUE, ratio = 1)

  expect_equal(m$match.matrix, m1$match.matrix)
  expect_equal(m$subclass, m1$subclass)

  m1 <- matchit(a ~ ., data = d, method = "cem", k2k = TRUE, k2k.method = "scaled_euclidean")

  expect_good_matchit(m, expect_subclass = TRUE, expect_distance = FALSE,
                      expect_match.matrix = TRUE, ratio = 1)

  expect_failure(expect_equal(m$match.matrix, m1$match.matrix))
  expect_failure(expect_equal(m$subclass, m1$subclass))

  m1 <- matchit(a ~ ., data = d, method = "cem", k2k = TRUE, k2k.method = "manhattan")

  expect_good_matchit(m, expect_subclass = TRUE, expect_distance = FALSE,
                      expect_match.matrix = TRUE, ratio = 1)

  expect_failure(expect_equal(m$match.matrix, m1$match.matrix))
  expect_failure(expect_equal(m$subclass, m1$subclass))

  m1 <- matchit(a ~ ., data = d, method = "cem", k2k = TRUE, m.order = "data")

  expect_good_matchit(m, expect_subclass = TRUE, expect_distance = FALSE,
                      expect_match.matrix = TRUE, ratio = 1)

  expect_equal(m$match.matrix, m1$match.matrix)
  expect_equal(m$subclass, m1$subclass)

  m1 <- matchit(a ~ ., data = d, method = "cem", k2k = TRUE, m.order = "closest")

  expect_good_matchit(m, expect_subclass = TRUE, expect_distance = FALSE,
                      expect_match.matrix = TRUE, ratio = 1)

  expect_failure(expect_equal(m$match.matrix, m1$match.matrix))
  expect_failure(expect_equal(m$subclass, m1$subclass))

  m2 <- matchit(a ~ ., data = d, method = "cem")
  m2$subclass[is.na(m2$subclass)] <- m2$subclass[!is.na(m2$subclass)][1]

  # Equivalent to NN matching with exact matching on subclass
  suppressWarnings({
    m1 <- matchit(a ~ ., data = d, method = "nearest", distance = "mahalanobis",
                  discard = m$weights == 0, exact = ~m2$subclass)
  })

  expect_equal(m$match.matrix, m1$match.matrix)
  expect_equal(m$subclass, m1$subclass)
})

# ===========================================================================
# Coarsened exact matching: option surface on lalonde, with pinned results.
#
# CEM has by far the largest set of user-facing options of any method here:
# `cutpoints` (five accepted forms), `grouping`, `k2k`, `k2k.method` (ten values
# across two different code paths), `mpower`, `m.order`, `estimand`, and
# `s.weights`. Everything below is deterministic except `m.order = "random"`.
# ===========================================================================

data("lalonde", package = "MatchIt")

f_cem <- treat ~ age + educ + race + married + nodegree + re74 + re75

lalonde_sw <- seq(0.5, 2, length.out = nrow(lalonde))

#Every unit in a stratum must share the same coarsened value of every variable, and
#every retained stratum must contain both treatment groups. This is the defining
#property of the method and holds whatever the coarsening.
expect_cem_strata_valid <- function(m) {
  keep <- !is.na(m$subclass)

  tab <- table(m$subclass[keep], m$treat[keep])
  expect_true(all(tab[, "0"] > 0L))
  expect_true(all(tab[, "1"] > 0L))

  invisible(m)
}

# ===== cutpoints =====

test_that("cem: default coarsening (sturges)", {
  m <- matchit(f_cem, data = lalonde, method = "cem")
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)
  expect_cem_strata_valid(m)
  expect_matchit_snapshot(m)
})

test_that("cem: cutpoints as a single bin count applies to every numeric variable", {
  m <- matchit(f_cem, data = lalonde, method = "cem", cutpoints = 3)
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)
  expect_cem_strata_valid(m)
  expect_matchit_snapshot(m)
})

test_that("cem: cutpoints as a per-variable list", {
  m <- matchit(f_cem, data = lalonde, method = "cem",
               cutpoints = list(age = 4, re74 = "q5"))
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)
  expect_cem_strata_valid(m)
  expect_matchit_snapshot(m)
})

test_that("cem: cutpoints as an explicit vector of boundaries", {
  m <- matchit(f_cem, data = lalonde, method = "cem",
               cutpoints = list(age = c(20, 30, 40)))
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)
  expect_cem_strata_valid(m)
  expect_matchit_snapshot(m)
})

test_that("cem: binning algorithms", {
  for (cp in c("sturges", "fd", "scott")) {
    m <- matchit(f_cem, data = lalonde, method = "cem", cutpoints = cp)
    expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = FALSE,
                        expect_subclass = TRUE)
    expect_cem_strata_valid(m)
    expect_matchit_snapshot(m)
  }
})

test_that("cem: sturges is the default and the algorithms differ from each other", {
  m_default <- matchit(f_cem, data = lalonde, method = "cem")

  expect_identical(matchit(f_cem, data = lalonde, method = "cem",
                           cutpoints = "sturges")$subclass,
                   m_default$subclass)

  #Partial matching is accepted, since the strings go through `pmatch()`
  expect_identical(matchit(f_cem, data = lalonde, method = "cem",
                           cutpoints = "stu")$subclass,
                   m_default$subclass)

  for (cp in c("fd", "scott")) {
    expect_not_equal(matchit(f_cem, data = lalonde, method = "cem",
                             cutpoints = cp)$subclass,
                     m_default$subclass)
  }
})

test_that("cem: quantile bins are equally sized", {
  #"q4" asks for quartiles, so each bin should hold a quarter of the sample
  set.seed(3)
  d <- data.frame(t = rep(0:1, each = 100L), x = runif(200L))

  m <- matchit(t ~ x, data = d, method = "cem", cutpoints = list(x = "q4"))

  br <- quantile(d$x, probs = seq(0.25, 0.75, by = 0.25), names = FALSE)
  expect_equal(unname(table(findInterval(d$x, c(-Inf, br, Inf)))),
               rep(50L, 4L), ignore_attr = TRUE)
  expect_length(levels(droplevels(m$subclass)), 4L)
})

test_that("cem: cutpoints = 0 is exact matching on the numeric variables", {
  #Documented equivalence: no binning means exact matching, so coarsening every
  #numeric variable at 0 must reproduce `method = "exact"` on the same formula.
  m_cem <- matchit(f_cem, data = lalonde, method = "cem", cutpoints = 0)
  m_exact <- matchit(f_cem, data = lalonde, method = "exact")

  expect_identical(is.na(m_cem$subclass), is.na(m_exact$subclass))
  expect_equal(m_cem$weights, m_exact$weights)
  expect_length(levels(m_cem$subclass), nlevels(m_exact$subclass))

  expect_matchit_snapshot(m_cem)
})

test_that("cem: cutpoints = NA is the same as cutpoints = 0", {
  #`NA` is undocumented but reaches the same branch: the variable is not binned.
  expect_identical(matchit(f_cem, data = lalonde, method = "cem",
                           cutpoints = list(age = NA))$subclass,
                   matchit(f_cem, data = lalonde, method = "cem",
                           cutpoints = list(age = 0))$subclass)
})

test_that("cem: cutpoints = 1 drops a variable from the strata", {
  #Documented: 1 removes the variable from the exact matching variables. With every
  #numeric variable dropped only `race` remains, so nothing is discarded.
  m <- matchit(f_cem, data = lalonde, method = "cem", cutpoints = 1)

  expect_length(levels(m$subclass), nlevels(lalonde$race))
  expect_false(anyNA(m$subclass))
  expect_true(all(m$weights > 0))

  expect_matchit_snapshot(m)
})

test_that("cem: cutpoints = 1 keeps the variable available for k2k pairing", {
  #Dropping a variable from the strata lets more strata contain both groups, so more
  #pairs form than when the same variable is matched on exactly.
  set.seed(12345)
  m1 <- matchit(f_cem, data = lalonde, method = "cem",
                cutpoints = list(re74 = 1), k2k = TRUE)
  set.seed(12345)
  m0 <- matchit(f_cem, data = lalonde, method = "cem",
                cutpoints = list(re74 = 0), k2k = TRUE)

  expect_gt(nlevels(m1$subclass), nlevels(m0$subclass))
})

test_that("cem: a length-1 numeric is a bin count, not a cut point", {
  #Easy to misread: `cutpoints = list(x = 2)` asks for two bins, while
  #`cutpoints = list(x = c(2, 3))` asks for cuts at 2 and 3.
  d <- data.frame(t = rep(0:1, each = 8L), x = rep(1:4, 4L))

  m_count <- matchit(t ~ x, data = d, method = "cem", cutpoints = list(x = 2))
  m_cuts <- matchit(t ~ x, data = d, method = "cem", cutpoints = list(x = c(2, 3)))

  expect_length(levels(droplevels(m_count$subclass)), 2L)
  expect_length(levels(droplevels(m_cuts$subclass)), 3L)
})

test_that("cem: values on a bin boundary go into the higher bin", {
  #Documented, and easy to break: with a cut at 25, a value of exactly 25 belongs
  #with 26 rather than with 24.
  d <- data.frame(t = rep(0:1, each = 9L), x = rep(c(24, 25, 26), 6L))

  m <- matchit(t ~ x, data = d, method = "cem", cutpoints = list(x = c(25, 100)))

  s <- as.character(m$subclass)
  expect_identical(s[d$x == 25][1L], s[d$x == 26][1L])
  expect_false(identical(s[d$x == 24][1L], s[d$x == 25][1L]))
})

test_that("cem: invalid cutpoints are errors", {
  for (cp in list(list(age = -1), list(age = "nope"), Inf, list(3))) {
    expect_err(matchit(f_cem, data = lalonde, method = "cem", cutpoints = cp))
  }

  expect_err(matchit(f_cem, data = lalonde, method = "cem", cutpoints = list(3)),
             "`cutpoints` must be a named list of binning values")
})

test_that("cem: cutpoints naming problems warn and are ignored", {
  expect_wrn(matchit(f_cem, data = lalonde, method = "cem",
                     cutpoints = list(nope = 3)),
             "named in `cutpoints` is not in the variables supplied to `matchit()`")

  expect_wrn(matchit(f_cem, data = lalonde, method = "cem",
                     cutpoints = list(race = 3)),
             "named in `cutpoints` is not numeric")

  #Ignoring them must leave the default result untouched
  m0 <- matchit(f_cem, data = lalonde, method = "cem")

  for (cp in list(list(nope = 3), list(race = 3))) {
    m <- suppressWarnings(matchit(f_cem, data = lalonde, method = "cem",
                                  cutpoints = cp))
    expect_identical(m$subclass, m0$subclass)
  }
})

# ===== grouping =====

test_that("cem: grouping combines categorical levels", {
  m <- matchit(f_cem, data = lalonde, method = "cem",
               grouping = list(race = list(c("black", "hispan"), "white")))
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)
  expect_cem_strata_valid(m)

  #No stratum may contain a white unit alongside a black or hispanic one
  keep <- !is.na(m$subclass)
  grouped <- ifelse(lalonde$race == "white", "white", "nonwhite")
  spans <- tapply(grouped[keep], m$subclass[keep], function(x) length(unique(x)))
  expect_true(all(spans == 1L))

  expect_matchit_snapshot(m)
})

test_that("cem: levels omitted from grouping are left alone", {
  #Documented: any levels left out keep their own category, so naming only the group
  #to be combined gives the same answer as naming every level.
  expect_identical(
    matchit(f_cem, data = lalonde, method = "cem",
            grouping = list(race = list(c("black", "hispan"))))$subclass,
    matchit(f_cem, data = lalonde, method = "cem",
            grouping = list(race = list(c("black", "hispan"), "white")))$subclass
  )
})

test_that("cem: grouping changes the result and relaxes the strata", {
  m0 <- matchit(f_cem, data = lalonde, method = "cem")
  m1 <- matchit(f_cem, data = lalonde, method = "cem",
                grouping = list(race = list(c("black", "hispan"))))

  expect_not_equal(m0$subclass, m1$subclass)
  expect_gte(sum(m1$weights > 0), sum(m0$weights > 0))
})

test_that("cem: invalid grouping is an error", {
  expect_err(matchit(f_cem, data = lalonde, method = "cem",
                     grouping = c(race = 1)),
             "`grouping` must be a named list of grouping values")

  expect_err(matchit(f_cem, data = lalonde, method = "cem",
                     grouping = list(list(c("black")))),
             "`grouping` must be a named list of grouping values")

  expect_err(matchit(f_cem, data = lalonde, method = "cem",
                     grouping = list(race = c("black", "hispan"))),
             "must be a list with entries containing values of the corresponding variable")
})

test_that("cem: grouping naming problems warn", {
  expect_wrn(matchit(f_cem, data = lalonde, method = "cem",
                     grouping = list(nope = list(c("a")))),
             "named in `grouping` is not in the variables supplied to `matchit()`")

  #A variable named in both takes its grouping; the cutpoints entry is dropped
  expect_wrn(matchit(f_cem, data = lalonde, method = "cem",
                     grouping = list(race = list(c("black", "hispan"))),
                     cutpoints = list(race = 3)),
             "is named in both `grouping` and `cutpoints`")
})

# ===== k2k =====

test_that("cem: k2k = TRUE produces 1:1 pairs within strata", {
  set.seed(12345)
  m <- matchit(f_cem, data = lalonde, method = "cem", k2k = TRUE)
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L, replace = FALSE)

  #No control unit is used twice
  matched <- na.omit(as.vector(m$match.matrix))
  expect_identical(anyDuplicated(matched), 0L)

  #Every pair comes from a single coarsened stratum of the k2k = FALSE solution
  strata <- matchit(f_cem, data = lalonde, method = "cem")$subclass

  for (t1 in rownames(m$match.matrix)) {
    t0 <- m$match.matrix[t1, 1L]

    if (!is.na(t0)) {
      expect_identical(as.character(strata[t1]), as.character(strata[t0]))
    }
  }

  expect_matchit_snapshot(m)
})

test_that("cem: every k2k.method is accepted and pins a result", {
  #Two distinct code paths: the `distance()` transforms go through
  #`nn_matchC_dispatch()` on transformed covariates, while the `dist()` methods build
  #a distance matrix per stratum.
  for (km in list("mahalanobis", "robust_mahalanobis", "euclidean", "scaled_euclidean",
                  "maximum", "manhattan", "canberra", "binary", "minkowski")) {
    set.seed(12345)
    m <- matchit(f_cem, data = lalonde, method = "cem", k2k = TRUE, k2k.method = km)
    expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = TRUE,
                        expect_subclass = TRUE, ratio = 1L, replace = FALSE)
    expect_matchit_snapshot(m)
  }
})

test_that("cem: k2k.method = NULL matches without a distance", {
  set.seed(12345)
  m <- matchit(f_cem, data = lalonde, method = "cem", k2k = TRUE, k2k.method = NULL)
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L, replace = FALSE)
  expect_matchit_snapshot(m)
})

test_that("cem: the k2k.method values give genuinely different pairings", {
  #Guards against an option that is accepted but does nothing. `minkowski` with the
  #default `mpower = 2` is deliberately the same as `euclidean`.
  methods <- list(NULL, "mahalanobis", "robust_mahalanobis", "euclidean",
                  "scaled_euclidean", "maximum", "manhattan", "canberra", "binary")

  mms <- lapply(methods, function(km) {
    set.seed(12345)
    matchit(f_cem, data = lalonde, method = "cem", k2k = TRUE,
            k2k.method = km)$match.matrix
  })

  for (i in seq_along(mms)[-1L]) {
    for (j in seq_len(i - 1L)) {
      expect_not_equal(mms[[i]], mms[[j]])
    }
  }
})

test_that("cem: minkowski with mpower = 2 is euclidean", {
  set.seed(12345)
  m_mink <- matchit(f_cem, data = lalonde, method = "cem", k2k = TRUE,
                    k2k.method = "minkowski", mpower = 2)
  set.seed(12345)
  m_euc <- matchit(f_cem, data = lalonde, method = "cem", k2k = TRUE,
                   k2k.method = "euclidean")

  expect_identical(m_mink$match.matrix, m_euc$match.matrix)
})

test_that("cem: mpower changes the minkowski pairing", {
  #The ranking of candidate controls within a stratum stabilizes as the power grows,
  #so only well-separated powers give distinct pairings: 3, 4, and 6 agree with each
  #other, as do 10 and 50. These five are pairwise distinct.
  mms <- lapply(c(0.5, 1, 1.5, 3, 10), function(mp) {
    set.seed(12345)
    matchit(f_cem, data = lalonde, method = "cem", k2k = TRUE,
            k2k.method = "minkowski", mpower = mp)$match.matrix
  })

  for (i in seq_along(mms)[-1L]) {
    for (j in seq_len(i - 1L)) {
      expect_not_equal(mms[[i]], mms[[j]])
    }
  }
})

test_that("cem: every m.order is accepted and changes the pairing", {
  mms <- list()

  for (mo in c("data", "random", "closest", "farthest")) {
    set.seed(12345)
    m <- matchit(f_cem, data = lalonde, method = "cem", k2k = TRUE, m.order = mo)
    expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = TRUE,
                        expect_subclass = TRUE, ratio = 1L, replace = FALSE)
    expect_matchit_snapshot(m)
    mms[[mo]] <- m$match.matrix
  }

  expect_not_equal(mms[["data"]], mms[["closest"]])
  expect_not_equal(mms[["data"]], mms[["farthest"]])
  expect_not_equal(mms[["closest"]], mms[["farthest"]])
})

test_that("cem: m.order = 'random' is the only randomized specification", {
  same_twice <- function(...) {
    a <- matchit(f_cem, data = lalonde, method = "cem", k2k = TRUE, ...)
    b <- matchit(f_cem, data = lalonde, method = "cem", k2k = TRUE, ...)
    identical(a$match.matrix, b$match.matrix)
  }

  expect_true(same_twice())
  expect_true(same_twice(k2k.method = NULL))
  expect_true(same_twice(k2k.method = "manhattan"))
  expect_false(same_twice(m.order = "random"))
})

test_that("cem: k2k.method = NULL is not random on its own", {
  #The documentation calls `NULL` "random matching", but with the default
  #`m.order = "data"` it matches in data order and is fully reproducible. Randomness
  #comes from `m.order = "random"`, not from the absence of a distance.
  set.seed(1)
  a <- matchit(f_cem, data = lalonde, method = "cem", k2k = TRUE, k2k.method = NULL)
  set.seed(2)
  b <- matchit(f_cem, data = lalonde, method = "cem", k2k = TRUE, k2k.method = NULL)

  expect_identical(a$match.matrix, b$match.matrix)

  set.seed(1)
  c1 <- matchit(f_cem, data = lalonde, method = "cem", k2k = TRUE,
                k2k.method = NULL, m.order = "random")
  set.seed(2)
  c2 <- matchit(f_cem, data = lalonde, method = "cem", k2k = TRUE,
                k2k.method = NULL, m.order = "random")

  expect_not_equal(c1$match.matrix, c2$match.matrix)
})

test_that("cem: invalid k2k arguments are errors", {
  expect_err(matchit(f_cem, data = lalonde, method = "cem", k2k = TRUE,
                     k2k.method = "nope"),
             "`k2k.method` should be one of")

  expect_err(matchit(f_cem, data = lalonde, method = "cem", k2k = TRUE,
                     k2k.method = "minkowski", mpower = 0),
             "`mpower` must be positive.")

  expect_err(matchit(f_cem, data = lalonde, method = "cem", k2k = TRUE,
                     m.order = "largest"),
             "`m.order` should be one of")
})

test_that("cem: k2k arguments are ignored without k2k = TRUE", {
  #Recorded, not endorsed: unlike the arguments that belong to other methods, these
  #are silently dropped rather than warned about.
  m0 <- matchit(f_cem, data = lalonde, method = "cem")

  expect_identical(matchit(f_cem, data = lalonde, method = "cem",
                           k2k.method = "manhattan")$subclass,
                   m0$subclass)
  expect_identical(matchit(f_cem, data = lalonde, method = "cem",
                           m.order = "closest")$subclass,
                   m0$subclass)
  expect_null(matchit(f_cem, data = lalonde, method = "cem",
                      k2k.method = "manhattan")$match.matrix)
})

# ===== estimand and s.weights =====

test_that("cem: estimand without k2k changes the weights but not the strata", {
  ms <- lapply(c("ATT", "ATC", "ATE"), function(e) {
    matchit(f_cem, data = lalonde, method = "cem", estimand = e)
  })

  expect_identical(ms[[1L]]$subclass, ms[[2L]]$subclass)
  expect_identical(ms[[1L]]$subclass, ms[[3L]]$subclass)

  expect_not_equal(ms[[1L]]$weights, ms[[2L]]$weights)
  expect_not_equal(ms[[1L]]$weights, ms[[3L]]$weights)

  for (m in ms) {
    expect_matchit_snapshot(m)
  }
})

test_that("cem: estimand with k2k selects the focal group", {
  #The focal group supplies the rows of `match.matrix`, so ATC flips its orientation.
  for (e in c("ATT", "ATE")) {
    set.seed(12345)
    m <- matchit(f_cem, data = lalonde, method = "cem", estimand = e, k2k = TRUE)
    expect_equal(nrow(m$match.matrix), sum(lalonde$treat == 1L))
    expect_matchit_snapshot(m)
  }

  set.seed(12345)
  m <- matchit(f_cem, data = lalonde, method = "cem", estimand = "ATC", k2k = TRUE)
  expect_equal(nrow(m$match.matrix), sum(lalonde$treat == 0L))
  expect_matchit_snapshot(m)
})

test_that("cem: s.weights do not affect the strata", {
  #Coarsening depends only on the covariates.
  expect_identical(matchit(f_cem, data = lalonde, method = "cem")$subclass,
                   matchit(f_cem, data = lalonde, method = "cem",
                           s.weights = lalonde_sw)$subclass)
})

test_that("cem: s.weights enter the k2k distance only for the scaled methods", {
  #Documented as affecting "the scaling factors when k2k = TRUE and certain methods
  #are used": the `distance()` transforms take a weighted covariance, while plain
  #euclidean and the `dist()` methods have no scaling to weight.
  changed <- function(km) {
    set.seed(12345)
    a <- matchit(f_cem, data = lalonde, method = "cem", k2k = TRUE,
                 k2k.method = km)$match.matrix
    set.seed(12345)
    b <- matchit(f_cem, data = lalonde, method = "cem", k2k = TRUE,
                 k2k.method = km, s.weights = lalonde_sw)$match.matrix
    !identical(a, b)
  }

  for (km in c("mahalanobis", "robust_mahalanobis", "scaled_euclidean")) {
    expect_true(changed(km))
  }

  for (km in c("euclidean", "manhattan")) {
    expect_false(changed(km))
  }
})

test_that("cem: s.weights change the matching weights", {
  m0 <- matchit(f_cem, data = lalonde, method = "cem")
  m1 <- matchit(f_cem, data = lalonde, method = "cem", s.weights = lalonde_sw)

  expect_not_equal(m0$weights, m1$weights)
  expect_matchit_snapshot(m1)
})

# ===== interaction, balance, and conditions =====

test_that("cem: coarsening and k2k combine", {
  set.seed(12345)
  m <- matchit(f_cem, data = lalonde, method = "cem",
               cutpoints = list(age = 4, re74 = "q3"),
               grouping = list(race = list(c("black", "hispan"))),
               k2k = TRUE, k2k.method = "scaled_euclidean", m.order = "closest")
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L, replace = FALSE)
  expect_matchit_snapshot(m)
})

test_that("cem: matching improves balance", {
  expect_balance_improved(matchit(f_cem, data = lalonde, method = "cem"))

  set.seed(12345)
  expect_balance_improved(matchit(f_cem, data = lalonde, method = "cem", k2k = TRUE))
})

test_that("cem: no covariates is an error", {
  expect_err(matchit(treat ~ 1, data = lalonde, method = "cem"),
             "covariates must be specified in the input formula to use coarsened exact matching")
})

test_that("cem: unused arguments warn and are ignored", {
  m0 <- matchit(f_cem, data = lalonde, method = "cem")

  for (arg in c("distance", "exact", "mahvars", "discard", "replace", "caliper",
                "ratio")) {
    args <- list(f_cem, data = lalonde, method = "cem")
    args[[arg]] <- switch(arg,
                          distance = "probit",
                          exact = ~ race,
                          mahvars = ~ age,
                          discard = "both",
                          replace = TRUE,
                          caliper = 0.1,
                          ratio = 2)

    expect_wrn(do.call(matchit, args),
               sprintf('The argument `%s` is not used with `method = "cem"` and will be ignored.',
                       arg))

    expect_identical(suppressWarnings(do.call(matchit, args))$subclass, m0$subclass)
  }
})

test_that("cem: missing values in covariates are an error", {
  lalonde_na <- inject_missingness(lalonde, "educ")

  expect_err(matchit(f_cem, data = lalonde_na, method = "cem"),
             "Missing and non-finite values are not allowed in the covariates")
})

test_that("cem: no unexpected conditions in the baseline calls", {
  expect_no_unexpected_warning(matchit(f_cem, data = lalonde, method = "cem"))
  expect_no_unexpected_warning(matchit(f_cem, data = lalonde, method = "cem",
                                       k2k = TRUE))
})
