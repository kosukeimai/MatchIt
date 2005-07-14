distance2rpart <- function(formula, data, ...) {
  require(rpart)
  res <- rpart(formula, data, ...)
  return(list(model = res, distance = predict(res)))
}

