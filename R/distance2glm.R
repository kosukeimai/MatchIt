distance2logit <- function(formula, data, ...) {
  res <- list()
  res$out <- glm(formula, data, family=binomial(logit), ...)
  res$dis <- predict(res$out, type="response")
  return(res)
}

distance2probit <- function(formula, data, ...) {
  res <- list()
  res$out <- glm(formula, data, family=binomial(probit), ...)
  res$dis <- predict(res$out, type="response")
  return(res)
}

distance2cloglog <- function(formula, data, ...) {
  res <- list()
  res$out <- glm(formula, data, family=binomial(cloglog), ...)
  res$dis <- predict(res$out, type="response")
  return(res)
}

distance2log <- function(formula, data, ...) {
  res <- list()
  res$out <- glm(formula, data, family=binomial(log), ...)
  res$dis <- predict(res$out, type="response")
  return(res)
}

distance2cauchit <- function(formula, data, ...) {
  res <- list()
  res$out <- glm(formula, data, family=binomial(cauchit), ...)
  res$dis <- predict(res$out, type="response")
  return(res)
}

