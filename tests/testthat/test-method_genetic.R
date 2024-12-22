test_that("calipers work, positive and negative", {
  set.seed(1234)
  n <- 1e3
  p <- runif(n, 0, .4)
  x <- matrix(runif(n * 4), nrow = n)
  g <- sample(1:5, n, TRUE)
  a <- rbinom(n, 1, p)
  u <- 1:n; u[a == 0] <- sample(u[a == 0][1:round(sum(a == 0)/5)], sum(a == 0), replace = TRUE)
  dis <- as.logical(rbinom(n, 1, .1))
  d <- data.frame(p, x, a, g, u, dis)

  #Positive calipers
  m <- matchit(a ~ X1 + X2 + X3 + X4, data = d, distance = d$p, method = "genetic",
               pop.size = 10, max.generations = 5,
               caliper = .001, std.caliper = FALSE, replace = TRUE, ratio = 3)

  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = !m$info$replace,
                      ratio = m$info$ratio)

  m <- matchit(a ~ X1 + X2 + X3 + X4, data = d, distance = d$p, method = "genetic",
               pop.size = 10, max.generations = 5, caliper = c(X1 = .001),
               std.caliper = FALSE, replace = TRUE, ratio = 3)

  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = !m$info$replace,
                      ratio = m$info$ratio)

  m <- matchit(a ~ X1 + X2 + X3 + X4, data = d, distance = d$p, method = "genetic",
               pop.size = 10, max.generations = 5, caliper = c(.02, X1 = .01),
               std.caliper = FALSE, replace = TRUE, ratio = 3)

  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = !m$info$replace,
                      ratio = m$info$ratio)

  #Negative calipers
  m <- matchit(a ~ X1 + X2 + X3 + X4, data = d, distance = d$p, method = "genetic",
               pop.size = 10, max.generations = 5, caliper = -.1,
               std.caliper = FALSE, replace = TRUE, ratio = 3)

  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = !m$info$replace,
                      ratio = m$info$ratio)

  m <- matchit(a ~ X1 + X2 + X3 + X4, data = d, distance = d$p, method = "genetic",
               pop.size = 10, max.generations = 5, caliper = c(X1 = -.1),
               std.caliper = FALSE, replace = TRUE, ratio = 3)

  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = !m$info$replace,
                      ratio = m$info$ratio)

  m <- matchit(a ~ X1 + X2 + X3 + X4, data = d, distance = d$p, method = "genetic",
               pop.size = 10, max.generations = 5, caliper = c(-.2, X1 = -.1),
               std.caliper = FALSE, replace = TRUE, ratio = 3)

  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = !m$info$replace,
                      ratio = m$info$ratio)

  m <- matchit(a ~ X1 + X2 + X3 + X4, data = d, distance = d$p, method = "genetic",
               pop.size = 10, max.generations = 5, caliper = c(-.02, X1 = .001),
               std.caliper = FALSE, replace = TRUE, ratio = 3)

  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = !m$info$replace,
                      ratio = m$info$ratio)

  m <- matchit(a ~ X1 + X2 + X3 + X4, data = d, distance = d$p, method = "genetic",
               pop.size = 10, max.generations = 5, caliper = c(.002, X1 = -.1),
               std.caliper = FALSE, replace = TRUE, ratio = 3)

  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = !m$info$replace,
                      ratio = m$info$ratio)
})
