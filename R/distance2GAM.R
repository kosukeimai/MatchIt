distance2GAMlogit <- function(formula, data, ...) {
  require(mgcv)
  res <- gam(formula, data, family=binomial(logit), ...)
  return(list(assign.model = res, pscore = fitted(res)))
}

distance2GAMprobit <- function(formula, data, ...) {
  require(mgcv)
  res <- gam(formula, data, family=binomial(probit), ...)
  return(list(assign.model = res, pscore = fitted(res)))
}

distance2GAMcloglog <- function(formula, data, ...) {
  require(mgcv)
  res <- gam(formula, data, family=binomial(cloglog), ...)
  return(list(assign.model = res, pscore = fitted(res)))
}

distance2GAMlog <- function(formula, data, ...) {
  require(mgcv)
  res <- gam(formula, data, family=binomial(log), ...)
  return(list(assign.model = res, pscore = fitted(res)))
}

distance2GAMcauchit <- function(formula, data, ...) {
  require(mgcv)
  res <- gam(formula, data, family=binomial(cauchit), ...)
  return(list(assign.model = res, pscore = fitted(res)))
}

