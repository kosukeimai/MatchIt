distance2rpart <- function(formula, data, ...) {
  require(rpart)
  res <- rpart(formula, data, ...)
  return(list(assign.model = res, pscore = predict(res)))
}

