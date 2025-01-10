test_that("distance vector, mah vars, and distance matrix yield identical results", {
  set.seed(1234)
  n <- 1e3
  p <- runif(n, 0, .4)
  x <- runif(n)
  g <- sample(1:5, n, TRUE)
  a <- rbinom(n, 1, p)
  u <- 1:n; u[a == 0] <- sample(u[a == 0][1:round(sum(a == 0)/5)], sum(a == 0), replace = TRUE)
  dis <- as.logical(rbinom(n, 1, .1))
  d <- data.frame(p, x, a, g, u, dis)
  d$p_ <- d$p

  dd <- euclidean_dist(a ~ p, data = d)

  test_all <- function(..., which = 1:4) {

    M <- list()
    if (any(which == 1)) {
      matchit_try({
        m <- matchit(a ~ p + p_, data = d,
                     distance = d$p,
                     ...)
      }, dont_warn_if = "Fewer control")
      expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                          expect_subclass = !m$info$replace, replace = m$info$replace,
                          ratio = m$info$ratio)
      M <- c(M, list(m))
    }
    if (any(which == 2)) {
      matchit_try({
        m <- matchit(a ~ p + p_, data = d,
                     distance = "euclidean",
                     ...)
      }, dont_warn_if = "Fewer control")
      expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = TRUE,
                          expect_subclass = !m$info$replace, replace = m$info$replace,
                          ratio = m$info$ratio)
      M <- c(M, list(m))
    }
    if (any(which == 3)) {
      matchit_try({
        m <- matchit(a ~ p + p_, data = d,
                     distance = d$p,
                     mahvars = ~p + p_,
                     ...)
      }, dont_warn_if = "Fewer control")
      expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                          expect_subclass = !m$info$replace, replace = m$info$replace,
                          ratio = m$info$ratio)
      M <- c(M, list(m))
    }
    if (any(which == 4)) {
      matchit_try({
        m <- matchit(a ~ p + p_, data = d,
                     distance = dd,
                     ...)
      }, dont_warn_if = "Fewer control")
      expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = TRUE,
                          expect_subclass = !m$info$replace, replace = m$info$replace,
                          ratio = m$info$ratio)
      M <- c(M, list(m))
    }

    all(unlist(lapply(M[-1], function(m) isTRUE(all.equal(M[[1]]$match.matrix,
                                                          m$match.matrix)))))
  }

  expect_true(test_all(m.order = "data"))
  expect_true(test_all(m.order = "closest"))
  expect_true(test_all(m.order = "farthest"))
  expect_true(test_all(m.order = "largest", which = c(1, 3)))
  expect_true(test_all(m.order = "smallest", which = c(1, 3)))

  expect_true(test_all(m.order = "data", ratio = 2))
  expect_true(test_all(m.order = "closest", ratio = 2))

  expect_true(test_all(m.order = "data", ratio = 2, max.controls = 3, which = c(1, 3)))
  expect_true(test_all(m.order = "closest", ratio = 2, max.controls = 3, which = c(1, 3)))
  expect_true(test_all(m.order = "largest", ratio = 2, max.controls = 3, which = c(1, 3)))

  expect_true(test_all(m.order = "data", ratio = 2, replace = TRUE))
  expect_true(test_all(m.order = "closest", ratio = 2, replace = TRUE))

  expect_true(test_all(m.order = "data", ratio = 2, reuse.max = 3))
  expect_true(test_all(m.order = "closest", ratio = 2, reuse.max = 3))

  expect_true(test_all(m.order = "data", ratio = 2, caliper = .001, std.caliper = FALSE, which = c(1, 3)))
  expect_true(test_all(m.order = "closest", ratio = 2, caliper = .001, std.caliper = FALSE, which = c(1, 3)))
  expect_true(test_all(m.order = "largest", ratio = 2, caliper = .001, std.caliper = FALSE, which = c(1, 3)))

  expect_true(test_all(m.order = "data", ratio = 2, caliper = -.001, std.caliper = FALSE, which = c(1, 3)))
  expect_true(test_all(m.order = "closest", ratio = 2, caliper = -.001, std.caliper = FALSE, which = c(1, 3)))
  expect_true(test_all(m.order = "largest", ratio = 2, caliper = -.001, std.caliper = FALSE, which = c(1, 3)))

  expect_true(test_all(m.order = "data", ratio = 2, caliper = c(p = .001), std.caliper = FALSE))
  expect_true(test_all(m.order = "closest", ratio = 2, caliper = c(p = .001), std.caliper = FALSE))

  expect_true(test_all(m.order = "data", ratio = 2, caliper = c(p = -.001), std.caliper = FALSE))
  expect_true(test_all(m.order = "closest", ratio = 2, caliper = c(p = -.001), std.caliper = FALSE))

  expect_true(test_all(m.order = "data", ratio = 2, caliper = c(p = .001), std.caliper = FALSE, reuse.max = 3))
  expect_true(test_all(m.order = "closest", ratio = 2, caliper = c(p = .001), std.caliper = FALSE, reuse.max = 3))

  expect_true(test_all(m.order = "data", ratio = 2, caliper = c(p = -.001), std.caliper = FALSE, reuse.max = 3))
  expect_true(test_all(m.order = "closest", ratio = 2, caliper = c(p = -.001), std.caliper = FALSE, reuse.max = 3))

  expect_true(test_all(m.order = "data", ratio = 2, exact = ~g))
  expect_true(test_all(m.order = "closest", ratio = 2, exact = ~g))

  expect_true(test_all(m.order = "data", ratio = 2, exact = ~g, replace = TRUE))
  expect_true(test_all(m.order = "closest", ratio = 2, exact = ~g, replace = TRUE))

  expect_true(test_all(m.order = "data", ratio = 2, antiexact = ~g))
  expect_true(test_all(m.order = "closest", ratio = 2, antiexact = ~g))

  expect_true(test_all(m.order = "data", ratio = 2, antiexact = ~g, replace = TRUE))
  expect_true(test_all(m.order = "closest", ratio = 2, antiexact = ~g, replace = TRUE))

  expect_true(test_all(m.order = "data", ratio = 2, discard = dis))
  expect_true(test_all(m.order = "closest", ratio = 2, discard = dis))

  expect_true(test_all(m.order = "data", ratio = 2, unit.id = ~u))
  expect_true(test_all(m.order = "closest", ratio = 2, unit.id = ~u))

  expect_true(test_all(m.order = "data", ratio = 2, unit.id = ~u, reuse.max = 3))
  expect_true(test_all(m.order = "closest", ratio = 2, unit.id = ~u, reuse.max = 3))

  expect_true(test_all(m.order = "data", ratio = 2, unit.id = ~u, replace = TRUE))
  expect_true(test_all(m.order = "closest", ratio = 2, unit.id = ~u, replace = TRUE))
})

test_that("calipers work, positive and negative", {
  set.seed(1234)
  n <- 1e3
  p <- runif(n, 0, .4)
  x <- runif(n)
  g <- sample(1:5, n, TRUE)
  a <- rbinom(n, 1, p)
  u <- 1:n; u[a == 0] <- sample(u[a == 0][1:round(sum(a == 0)/5)], sum(a == 0), replace = TRUE)
  dis <- as.logical(rbinom(n, 1, .1))
  d <- data.frame(p, x, a, g, u, dis)

  #Positive calipers
  m <- matchit(a ~ p + x, data = d, distance = d$p, caliper = .001,
               std.caliper = FALSE, replace = TRUE, ratio = 3)

  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = !m$info$replace,
                      ratio = m$info$ratio)

  m <- matchit(a ~ p + x, data = d, distance = d$p, caliper = c(x = .001),
               std.caliper = FALSE, replace = TRUE, ratio = 3)

  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = !m$info$replace,
                      ratio = m$info$ratio)

  m <- matchit(a ~ p + x, data = d, distance = d$p, caliper = c(.02, x = .01),
               std.caliper = FALSE, replace = TRUE, ratio = 3)

  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = !m$info$replace,
                      ratio = m$info$ratio)

  #Negative calipers
  m <- matchit(a ~ p + x, data = d, distance = d$p, caliper = -.001,
               std.caliper = FALSE, replace = TRUE, ratio = 3)

  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = !m$info$replace,
                      ratio = m$info$ratio)

  m <- matchit(a ~ p + x, data = d, distance = d$p, caliper = c(x = -.001),
               std.caliper = FALSE, replace = TRUE, ratio = 3)

  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = !m$info$replace,
                      ratio = m$info$ratio)

  m <- matchit(a ~ p + x, data = d, distance = d$p, caliper = c(-.02, x = -.01),
               std.caliper = FALSE, replace = TRUE, ratio = 3)

  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = !m$info$replace,
                      ratio = m$info$ratio)

  m <- matchit(a ~ p + x, data = d, distance = d$p, caliper = c(-.02, x = .01),
               std.caliper = FALSE, replace = TRUE, ratio = 3)

  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = !m$info$replace,
                      ratio = m$info$ratio)

  m <- matchit(a ~ p + x, data = d, distance = d$p, caliper = c(.02, x = -.01),
               std.caliper = FALSE, replace = TRUE, ratio = 3)

  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = !m$info$replace,
                      ratio = m$info$ratio)
})
