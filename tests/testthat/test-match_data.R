# Tests for match_data(), its match.data() alias, and get_matches().
#
# These are the functions that hand the matched dataset to the user, so the columns they
# add, the rows they keep, and the arithmetic in the `weights` column are the parts that
# matter. Row content is checked directly rather than snapshotted; the matched sets
# themselves are already pinned by the per-method test files.

data("lalonde", package = "MatchIt")

f_md <- treat ~ age + educ + race + married + nodegree + re74 + re75

lalonde_sw <- seq(0.5, 2, length.out = nrow(lalonde))

# ===== match_data() =====

test_that("match_data: adds the expected columns and attributes", {
  m <- matchit(f_md, data = lalonde, method = "nearest")
  md <- match_data(m)

  expect_s3_class(md, "matchdata")
  expect_s3_class(md, "data.frame")

  expect_identical(setdiff(names(md), names(lalonde)),
                   c("distance", "weights", "subclass"))

  #The attributes record which column is which, for downstream functions
  expect_identical(attr(md, "distance"), "distance")
  expect_identical(attr(md, "weights"), "weights")
  expect_identical(attr(md, "subclass"), "subclass")

  #The added columns carry the values from the matchit object
  i <- match(rownames(md), rownames(lalonde))
  expect_equal(md$distance, m$distance[i], ignore_attr = TRUE)
  expect_equal(md$weights, m$weights[i], ignore_attr = TRUE)
  expect_equal(as.character(md$subclass), as.character(m$subclass[i]))
})

test_that("match_data: only columns the object supports are added", {
  #`nearest` with replacement produces no subclass, and the methods with no distance
  #measure produce no distance column.
  m_rep <- matchit(f_md, data = lalonde, method = "nearest", replace = TRUE)
  expect_null(m_rep$subclass)
  expect_identical(setdiff(names(match_data(m_rep)), names(lalonde)),
                   c("distance", "weights"))

  m_cem <- matchit(f_md, data = lalonde, method = "cem")
  expect_null(m_cem$distance)
  expect_identical(setdiff(names(match_data(m_cem)), names(lalonde)),
                   c("weights", "subclass"))

  m_card <- matchit(f_md, data = lalonde, method = "cardinality")
  expect_identical(setdiff(names(match_data(m_card)), names(lalonde)), "weights")
})

test_that("match_data: drop.unmatched", {
  m <- matchit(f_md, data = lalonde, method = "nearest")

  md <- match_data(m)
  expect_identical(nrow(md), sum(m$weights > 0))
  expect_true(all(md$weights > 0))

  md_all <- match_data(m, drop.unmatched = FALSE)
  expect_identical(nrow(md_all), nrow(lalonde))
  expect_true(any(md_all$weights == 0))
})

test_that("match_data: group selects a treatment group", {
  m <- matchit(f_md, data = lalonde, method = "nearest")

  md <- match_data(m)

  md_t <- match_data(m, group = "treated")
  expect_true(all(md_t$treat == 1L))
  expect_identical(nrow(md_t), sum(md$treat == 1L))

  md_c <- match_data(m, group = "control")
  expect_true(all(md_c$treat == 0L))
  expect_identical(nrow(md_c), sum(md$treat == 0L))

  expect_identical(nrow(md_t) + nrow(md_c), nrow(md))

  expect_err(match_data(m, group = "nope"),
             '`group` should be one of "all", "treated", or "control".')
})

test_that("match_data: the added columns can be renamed", {
  m <- matchit(f_md, data = lalonde, method = "nearest")
  md <- match_data(m, distance = "ps", weights = "w", subclass = "strat")

  expect_true(all(c("ps", "w", "strat") %in% names(md)))
  expect_false(any(c("distance", "weights", "subclass") %in% names(md)))

  #The attributes follow the new names
  expect_identical(attr(md, "distance"), "ps")
  expect_identical(attr(md, "weights"), "w")
  expect_identical(attr(md, "subclass"), "strat")
})

test_that("match_data: a name that already exists in the data is an error", {
  m <- matchit(f_md, data = lalonde, method = "nearest")

  for (nm in c("distance", "weights", "subclass")) {
    clash <- lalonde
    clash[[nm]] <- 1

    expect_err(match_data(m, data = clash),
               sprintf('"%s" is already the name of a variable in the data.', nm))
  }

  #Renaming resolves the clash
  clash <- lalonde
  clash$weights <- 1
  expect_no_error(match_data(m, data = clash, weights = "matching_weights"))
})

test_that("match_data: include.s.weights multiplies the matching weights", {
  m <- matchit(f_md, data = lalonde, method = "nearest", s.weights = lalonde_sw)

  md_with <- match_data(m, include.s.weights = TRUE)
  md_without <- match_data(m, include.s.weights = FALSE)

  i <- match(rownames(md_with), rownames(lalonde))

  expect_equal(md_with$weights, m$weights[i] * lalonde_sw[i], ignore_attr = TRUE)
  expect_equal(md_without$weights, m$weights[i], ignore_attr = TRUE)

  #TRUE is the default
  expect_equal(match_data(m)$weights, md_with$weights)
})

test_that("match_data: include.s.weights is a no-op without sampling weights", {
  m <- matchit(f_md, data = lalonde, method = "nearest")

  expect_equal(match_data(m, include.s.weights = TRUE)$weights,
               match_data(m, include.s.weights = FALSE)$weights)
})

test_that("match_data: a supplied data frame is used", {
  m <- matchit(f_md, data = lalonde, method = "nearest")

  lalonde_extra <- lalonde
  lalonde_extra$extra <- lalonde$age * 3

  md <- match_data(m, data = lalonde_extra)
  expect_true("extra" %in% names(md))
  expect_identical(nrow(md), sum(m$weights > 0))
})

test_that("match_data: a wrongly sized `data` is an error", {
  m <- matchit(f_md, data = lalonde, method = "nearest")

  #The message names both sizes so the mismatch is diagnosable
  expect_err(match_data(m, data = lalonde[1:10, ]),
             "`data` must be the original dataset")
  expect_err(match_data(m, data = lalonde[1:10, ]), "614 units")
  expect_err(match_data(m, data = lalonde[1:10, ]), "has 10 rows")

  #Even when the original dataset is in scope and would have been found
  expect_err(match_data(m, data = rbind(lalonde, lalonde)),
             "`data` must be the original dataset")

  #get_matches() validates through match_data()
  expect_err(get_matches(m, data = lalonde[1:10, ]),
             "`data` must be the original dataset")
})

test_that("match_data: a `data` that is not two-dimensional is an error", {
  m <- matchit(f_md, data = lalonde, method = "nearest")

  expect_err(match_data(m, data = as.list(lalonde)), "must be a data frame")
  expect_err(match_data(m, data = lalonde$age), "must be a data frame")

  #A matrix of the right size is accepted and coerced
  md <- match_data(m, data = as.matrix(lalonde[c("age", "educ")]))
  expect_s3_class(md, "data.frame")
  expect_named(md, c("age", "educ", "distance", "weights", "subclass"))
})

test_that("match_data: a supplied `data` takes precedence over the search", {
  m <- matchit(f_md, data = lalonde, method = "nearest")

  #A frame with the right number of rows is honored even if it holds other columns
  md <- match_data(m, data = lalonde[, 1:2])
  expect_identical(nrow(md), sum(m$weights > 0))
  expect_identical(ncol(md), 5L)

  #Supplying the original dataset is equivalent to letting it be found
  expect_identical(match_data(m, data = lalonde), match_data(m))
})

test_that("match_data: row order in a supplied `data` is trusted, not checked", {
  #Recorded so the assumption is explicit: the rows of `data` are assumed to line up
  #with the original, and a reordered frame of the right size is accepted as-is.
  m <- matchit(f_md, data = lalonde, method = "nearest")

  reversed <- lalonde[rev(seq_len(nrow(lalonde))), ]
  md <- match_data(m, data = reversed)

  expect_identical(nrow(md), sum(m$weights > 0))
  expect_not_equal(md$age, match_data(m)$age)
})

test_that("match_data: match.data() is an alias", {
  m <- matchit(f_md, data = lalonde, method = "nearest")

  expect_identical(match.data(m), match_data(m))
  expect_identical(match.data(m, group = "treated"),
                   match_data(m, group = "treated"))
})

test_that("match_data: works for every matching method", {
  fits <- list(
    nearest = matchit(f_md, data = lalonde, method = "nearest"),
    cem = matchit(f_md, data = lalonde, method = "cem"),
    exact = matchit(treat ~ age + educ + race, data = lalonde, method = "exact"),
    subclass = matchit(f_md, data = lalonde, method = "subclass"),
    `no matching` = matchit(f_md, data = lalonde, method = NULL)
  )

  if (rlang::is_installed("optmatch")) {
    fits$full <- matchit(f_md, data = lalonde, method = "full")
  }

  if (rlang::is_installed("highs")) {
    fits$cardinality <- matchit(f_md, data = lalonde, method = "cardinality")
  }

  for (nm in names(fits)) {
    md <- match_data(fits[[nm]])

    expect_s3_class(md, "matchdata")
    expect_true("weights" %in% names(md))
    expect_true(all(md$weights > 0))
    expect_gt(nrow(md), 0L)
  }
})

test_that("match_data: the weighted group totals match the matchit object", {
  #Whatever subsetting happens, the weights that come out must be the ones that went in.
  for (method in c("nearest", "cem", "subclass")) {
    m <- matchit(f_md, data = lalonde, method = method)
    md <- match_data(m)

    for (g in c(0L, 1L)) {
      expect_equal(sum(md$weights[md$treat == g]),
                   sum(m$weights[m$treat == g]))
    }
  }
})

# ===== get_matches() =====

test_that("get_matches: one row per unit per pair", {
  m <- matchit(f_md, data = lalonde, method = "nearest")
  gm <- get_matches(m)

  expect_s3_class(gm, "getmatches")
  expect_s3_class(gm, "data.frame")

  expect_identical(setdiff(names(gm), names(lalonde)),
                   c("id", "subclass", "weights", "distance"))

  expect_identical(attr(gm, "id"), "id")
  expect_identical(attr(gm, "subclass"), "subclass")
  expect_identical(attr(gm, "weights"), "weights")
  expect_identical(attr(gm, "distance"), "distance")

  #1:1 matching gives two rows per pair, and every subclass holds one of each group
  expect_identical(nrow(gm), 2L * sum(!is.na(m$match.matrix)))
  expect_true(all(table(gm$subclass) == 2L))
  expect_true(all(tapply(gm$treat, gm$subclass, sum) == 1L))
})

test_that("get_matches: reused units appear once per pair", {
  #This is the difference from match_data(): with replacement a control matched to two
  #treated units gets two rows, distinguished by `subclass` and sharing an `id`.
  m <- matchit(f_md, data = lalonde, method = "nearest", replace = TRUE)
  gm <- get_matches(m)

  expect_identical(nrow(gm), 2L * sum(!is.na(m$match.matrix)))
  expect_gt(nrow(gm), length(unique(gm$id)))
  expect_true(all(table(gm$subclass) == 2L))

  #Rows sharing an id are the same unit, so their covariates agree
  dup <- gm$id[duplicated(gm$id)][1L]
  rows <- gm[gm$id == dup, ]
  expect_length(unique(rows$age), 1L)
})

test_that("get_matches: ratio > 1 gives one row per treated unit per match", {
  m <- matchit(f_md, data = lalonde, method = "nearest", ratio = 2)
  gm <- get_matches(m)

  #Each subclass holds one treated unit and its matched controls
  expect_true(all(tapply(gm$treat, gm$subclass, sum) == 1L))
  expect_true(all(table(gm$subclass) <= 3L))
})

test_that("get_matches: requires a match.matrix", {
  skip_if_not_installed("optmatch")

  m <- matchit(f_md, data = lalonde, method = "full")
  expect_null(m$match.matrix)

  expect_err(get_matches(m),
             "component must be present in the <matchit> object")
})

test_that("get_matches: renaming and s.weights behave as in match_data", {
  m <- matchit(f_md, data = lalonde, method = "nearest", s.weights = lalonde_sw)

  gm <- get_matches(m, distance = "ps", weights = "w", subclass = "strat", id = "unit")
  expect_true(all(c("ps", "w", "strat", "unit") %in% names(gm)))

  gm_with <- get_matches(m, include.s.weights = TRUE)
  gm_without <- get_matches(m, include.s.weights = FALSE)
  expect_not_equal(gm_with$weights, gm_without$weights)
})

test_that("get_matches: no unexpected conditions", {
  m <- matchit(f_md, data = lalonde, method = "nearest")

  expect_no_unexpected_warning(get_matches(m))
  expect_no_unexpected_warning(match_data(m))
})
