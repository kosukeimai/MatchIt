distance2nnet <- function(formula, data, ...) {
  requireNamesapce(nnet)
  res <- nnet::nnet(formula, data, ...)
  return(list(model = res, distance = fitted(res)))
}

