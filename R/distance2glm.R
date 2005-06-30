distance2logit <- function(formula, data, ...) {
  res <- list()
  res$out <- glm(formula, data, family=binomial(logit), ...)
  res$dis <- fitted(res$out)
  return(res)
}

distance2probit <- function(formula, data, ...) {
  res <- list()
  res$out <- glm(formula, data, family=binomial(probit), ...)
  res$dis <- fitted(res$out)
  return(res)
}

distance2cloglog <- function(formula, data, ...) {
  res <- list()
  res$out <- glm(formula, data, family=binomial(cloglog), ...)
  res$dis <- fitted(res$out)
  return(res)
}

distance2log <- function(formula, data, ...) {
  res <- list()
  res$out <- glm(formula, data, family=binomial(log), ...)
  res$dis <- fitted(res$out)
  return(res)
}

distance2cauchit <- function(formula, data, ...) {
  res <- list()
  res$out <- glm(formula, data, family=binomial(cauchit), ...)
  res$dis <- fitted(res$out)
  return(res)
}

