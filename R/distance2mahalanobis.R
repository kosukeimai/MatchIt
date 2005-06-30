distance2mahalanobis <- function(formula, data, ...) {
  res <- list()
  X <- model.matrix(formula, data)
  Sigma <- var(X)
  res$dis <- mahalanobis(X, colMeans(X), cov(X))
  return(res)
}

