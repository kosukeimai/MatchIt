matchit <- function(formula, data, method = "nearest",
                    distance = "logit", distance.options=list(), ...) { 

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
  if (method == "exact") {
    dis <-  NULL
    if (!is.null(distance))
      warning("distance is set to `NULL' when exact matching is used.")
  }
  else {
    distance.options$formula <- formula
    distance.options$data <- data
    out1 <- do.call(fn1, distance.options)
    dis <- out1$dis
  }

  ## matching!
  out2 <- do.call(fn2, list(treat, X, dis, ...))
  
  return(out2)
}
