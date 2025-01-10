test_that("calipers work, positive and negative", {
  set.seed(1234)
  n <- 1e3
  p <- runif(n, 0, .4)
  x <- matrix(runif(n * 4), nrow = n)
  g <- sample(1:5, n, TRUE)
  a <- rbinom(n, 1, p)
  dis <- as.logical(rbinom(n, 1, .1))
  d <- data.frame(p, x, a, g, dis)

  #Positive calipers
  m <- matchit(a ~ X1 + X2 + X3 + X4, data = d, distance = d$p, method = "full",
               caliper = .001, std.caliper = FALSE)

  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)

  m <- matchit(a ~ X1 + X2 + X3 + X4, data = d, distance = d$p, method = "full",
               caliper = c(X1 = .001),
               std.caliper = FALSE)

  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)

  m <- matchit(a ~ X1 + X2 + X3 + X4, data = d, distance = d$p, method = "full",
               caliper = c(.02, X1 = .01),
               std.caliper = FALSE)

  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)

  #Negative calipers
  m <- matchit(a ~ X1 + X2 + X3 + X4, data = d, distance = d$p, method = "full",
               caliper = -.3,
               std.caliper = FALSE)

  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)

  m <- matchit(a ~ X1 + X2 + X3 + X4, data = d, distance = d$p, method = "full",
               caliper = c(X1 = -.5),
               std.caliper = FALSE)

  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)

  m <- matchit(a ~ X1 + X2 + X3 + X4, data = d, distance = d$p, method = "full",
               caliper = c(-.3, X1 = -.5),
               std.caliper = FALSE)

  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)

  m <- matchit(a ~ X1 + X2 + X3 + X4, data = d, distance = d$p, method = "full",
               caliper = c(-.3, X1 = .001),
               std.caliper = FALSE)

  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)

  m <- matchit(a ~ X1 + X2 + X3 + X4, data = d, distance = d$p, method = "full",
               caliper = c(.002, X1 = -.5),
               std.caliper = FALSE)

  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)
})