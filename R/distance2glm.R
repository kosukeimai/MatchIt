distance2logit <- function(formula, data, ...) 
  return(predict(glm(formula, data, family=binomial(logit), ...)))

distance2probit <- function(formula, data, ...) 
  return(predict(glm(formula, data, family=binomial(probit), ...)))

distance2cloglog <- function(formula, data, ...) 
  return(predict(glm(formula, data, family=binomial(cloglog), ...)))

distance2log <- function(formula, data, ...) 
  return(predict(glm(formula, data, family=binomial(log), ...)))

distance2cauchit <- function(formula, data, ...) 
  return(predict(glm(formula, data, family=binomial(cauchit), ...)))
