distance2rpart <- function(formula, data, ...) {
  requireNamespace(rpart)
  res <- rpart(formula, data, ...)
  return(list(model = res, distance = predict(res)))
}

