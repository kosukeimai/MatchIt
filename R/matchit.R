matchit <- function(formula, data, method = "nearest",
                    distance = "logit", distance.option=list(), ...) { 

  ## check inputs
  fn1 <- paste("distance2", distance, sep = "")
  if (!exists(fn1))
    stop(distance, "not supported.")
  fn2 <- paste("matchit2", method, sep = "")
  if (!exists(fn2))
    stop(method, "not supported.")

  ## obtain T and X
  tt <- terms(formula)
  treat <- model.response(model.frame(tt, data))
  X <- model.matrix(tt, data)
  
  ## estimate the distance measure
  distance.option$formula <- formula
  distance.option$data <- data
  dist <- do.call(fn1, distance.option)

  ## matching!
  res <- do.call(fn2, list(treat, dist, ...))
  
  return(res)
}
