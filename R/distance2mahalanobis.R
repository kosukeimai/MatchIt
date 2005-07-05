distance2mahalanobis <- function(formula, data, ...) {
  X <- model.matrix(formula, data)
  Sigma <- var(X)
  return(list(assign.model = NULL, pscore = mahalanobis(X,
                                     colMeans(X), cov(X))))
}

