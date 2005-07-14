distance2GAMlogit <- function(formula, data, ...) {
  require(mgcv)
  res <- gam(formula, data, family=binomial(logit), ...)
  return(list(model = res, distance = fitted(res)))
}

distance2GAMprobit <- function(formula, data, ...) {
  require(mgcv)
  res <- gam(formula, data, family=binomial(probit), ...)
  return(list(model = res, distance = fitted(res)))
}

distance2GAMcloglog <- function(formula, data, ...) {
  require(mgcv)
  res <- gam(formula, data, family=binomial(cloglog), ...)
  return(list(model = res, distance = fitted(res)))
}

distance2GAMlog <- function(formula, data, ...) {
  require(mgcv)
  res <- gam(formula, data, family=binomial(log), ...)
  return(list(model = res, distance = fitted(res)))
}

distance2GAMcauchit <- function(formula, data, ...) {
  require(mgcv)
  res <- gam(formula, data, family=binomial(cauchit), ...)
  return(list(model = res, distance = fitted(res)))
}

