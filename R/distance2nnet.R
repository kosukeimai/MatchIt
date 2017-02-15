distance2nnet <- function(formula, data, ...) {
  requireNamesapce(nnet)
  res <- nnet(formula, data, ...)
  return(list(model = res, distance = fitted(res)))
}

