matchit <- function(formula, data, method = "nearest", distance = "logit",
                    distance.options=list(), verbose = FALSE, ...) { 

  ## check inputs
  fn1 <- paste("distance2", distance, sep = "")
  if (!exists(fn1))
    stop(distance, "not supported.")
  fn2 <- paste("matchit2", method, sep = "")
  if (!exists(fn2))
    stop(method, "not supported.")
  
  ## obtain T and X
  tt <- terms(formula)
  mf <- model.frame(tt, data)
  treat <- model.response(mf)
  X <- model.matrix(tt, data=mf)

  ## estimate the distance measure
  if (method == "exact") {
    dis <-  out1 <- NULL
    if (!is.null(distance))
      warning("distance is set to `NULL' when exact matching is used.")
  }
  else {
    if (verbose)
    cat("Calculating distance measure via", distance, "\n")
    distance.options$formula <- formula
    distance.options$data <- data
    out1 <- do.call(fn1, distance.options)
    dis <- out1$dis
  }

  ## matching!
  if (verbose)
    cat("Matching via", method, "\n")
  if (method == "exact") out2 <- do.call(fn2, list(treat, X, dis, ...)) 
  else out2 <- do.call(fn2, list(treat, data, dis, ...))

  ## putting all the results together
  out2$call <- match.call()
  out2$assign <- out1
  out2$formula <- formula
  out2$treat <- treat
  out2$covariates <- X
  out2$pscore <- dis
  
  return(out2)
}
