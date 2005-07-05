distance2nnet <- function(formula, data, ...) {
  require(nnet)
  res <- nnet(formula, data, ...)
  return(list(assign.model = res, pscore = fitted(res)))
}

