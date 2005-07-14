distance2nnet <- function(formula, data, ...) {
  require(nnet)
  res <- nnet(formula, data, ...)
  return(list(model = res, distance = fitted(res)))
}

