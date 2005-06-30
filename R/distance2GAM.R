distance2GAMlogit <- function(formula, data, ...) {
  require(mgcv)
  res <- list()
  res$out <- gam(formula, data, family=binomial(logit), ...)
  res$dis <- predict(res$out, type="response")
  return(res)
}

distance2GAMprobit <- function(formula, data, ...) {
  require(mgcv)
  res <- list()
  res$out <- gam(formula, data, family=binomial(probit), ...)
  res$dis <- predict(res$out, type="response")
  return(res)
}

distance2GAMcloglog <- function(formula, data, ...) {
  require(mgcv)
  res <- list()
  res$out <- gam(formula, data, family=binomial(cloglog), ...)
  res$dis <- fitted(res$out)
  return(res)
}

distance2GAMlog <- function(formula, data, ...) {
  require(mgcv)
  res <- list()
  res$out <- gam(formula, data, family=binomial(log), ...)
  res$dis <- fitted(res$out)
  return(res)
}

distance2GAMcauchit <- function(formula, data, ...) {
  res <- list()
  res$out <- gam(formula, data, family=binomial(cauchit), ...)
  res$dis <- fitted(res$out)
  return(res)
}

