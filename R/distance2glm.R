distance2logit <- function(formula, data, ...) {
  res <- glm(formula, data, family=binomial(logit), ...)
  return(list(assign.model = res, pscore = fitted(res)))
}

distance2probit <- function(formula, data, ...) {
  res <- glm(formula, data, family=binomial(probit), ...)
  return(list(assign.model = res, pscore = fitted(res)))
}

distance2cloglog <- function(formula, data, ...) {
  res <- glm(formula, data, family=binomial(cloglog), ...)
  return(list(assign.model = res, pscore = fitted(res)))
}

distance2log <- function(formula, data, ...) {
  res <- glm(formula, data, family=binomial(log), ...)
  return(list(assign.model = res, pscore = fitted(res)))
}

distance2cauchit <- function(formula, data, ...) {
  res <- glm(formula, data, family=binomial(cauchit), ...)
  return(list(assign.model = res, pscore = fitted(res)))
}

