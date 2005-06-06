matchit <- function(formula, method="nearest", data, ...) {

  fn <- paste("matchit2", method, sep = "")
  if (!exists(fn))
    stop(method, "not supported.")
  mf <- match.call(expand.dots = TRUE)  
  res <- do.call(fn, list(formula, method, data, ...))
  
  return(res)
}
