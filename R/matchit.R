matchit <- function(formula, data, method = "nearest", distance = "logit",
                    distance.options=list(), discard = "none",
                    reestiamte = FALSE, verbose = FALSE, ...) { 

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
    pscore <- out1 <- discarded <- NULL
    if (!is.null(distance))
      warning("distance is set to `NULL' when exact matching is used.")
  }
  else {
    if (verbose)
      cat("Calculating distance measure via", distance, "\n")
    if (is.null(distance.options$formula))
      distance.options$formula <- formula
    if (is.null(distance.options$data))
      distance.options$data <- data
    out1 <- do.call(fn1, distance.options)
    discarded <- discard(treat, out1$pscore, discard)
    if (reestimate) {
      distance.options$data <- data[!discarded,]
      out1 <- do.call(fn1, distance.options)
    }
    pscore <- out1$pscore
  }

  ## matching!
  if (verbose)
    cat("Matching via", method, "\n")
  out2 <- do.call(fn2, list(treat, X, data, pscore, discarded, ...)) 

  ## putting all the results together
  out2$call <- match.call()
  out2$assign.model <- out1$assign.model
  out2$formula <- formula
  out2$treat <- treat
  out2$X <- X
  out2$pscore <- pscore
  out2$discarded <- discarded
  
  return(out2)
}
