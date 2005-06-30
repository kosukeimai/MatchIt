distance2rpart <- function(formula, data, ...) {
  require(rpart)
  res <- list()
  res$out <- rpart(formula, data, ...)
  res$dis <- predict(res$out)
  return(res)
}

