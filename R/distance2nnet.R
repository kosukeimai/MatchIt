distance2nnet <- function(formula, data, ...) {
  require(nnet)
  res <- list()
  res$out <- nnet(formula, data, ...)
  res$dis <- predict(res$out)
  return(res)
}

