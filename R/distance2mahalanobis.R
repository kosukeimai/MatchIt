distance2mahalanobis <- function(formula, data, ...) {
  X <- model.matrix(formula, data)
  Sigma <- var(X)
  return(list(model = NULL,
              distance = mahalanobis(X, colMeans(X), cov(X))))
}

