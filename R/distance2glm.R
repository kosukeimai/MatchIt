distance2logit <- function(formula, data, ...) {
  res <- glm(formula, data, family=binomial(logit), ...)
  return(list(model = res, distance = fitted(res)))
}

distance2linear.logit <- function(formula, data, ...) {
  res <- glm(formula, data, family=binomial(logit), ...)
  return(list(model = res, distance = predict(res)))
}

distance2probit <- function(formula, data, ...) {
  res <- glm(formula, data, family=binomial(probit), ...)
  return(list(model = res, distance = fitted(res)))
}

distance2linear.probit <- function(formula, data, ...) {
  res <- glm(formula, data, family=binomial(probit), ...)
  return(list(model = res, distance = predict(res)))
}

distance2cloglog <- function(formula, data, ...) {
  res <- glm(formula, data, family=binomial(cloglog), ...)
  return(list(model = res, distance = fitted(res)))
}

distance2linear.cloglog <- function(formula, data, ...) {
  res <- glm(formula, data, family=binomial(cloglog), ...)
  return(list(model = res, distance = predict(res)))
}

distance2log <- function(formula, data, ...) {
  res <- glm(formula, data, family=binomial(log), ...)
  return(list(model = res, distance = fitted(res)))
}

distance2linear.log <- function(formula, data, ...) {
  res <- glm(formula, data, family=binomial(log), ...)
  return(list(model = res, distance = predict(res)))
}

distance2cauchit <- function(formula, data, ...) {
  res <- glm(formula, data, family=binomial(cauchit), ...)
  return(list(model = res, distance = fitted(res)))
}

distance2linearcauchit <- function(formula, data, ...) {
  res <- glm(formula, data, family=binomial(cauchit), ...)
  return(list(model = res, distance = predict(res)))
}

