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
