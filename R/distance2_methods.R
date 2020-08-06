#distance2glm-----------------
distance2glm <- function(formula, data, link = "logit", ...) {

  if (!is.null(link) && startsWith(as.character(link), "linear")) {
    linear <- TRUE
    link <- sub("linear.", "", as.character(link))
  }
  else linear <- FALSE

  A <- list(...)
  A[!names(A) %in% c(names(formals(glm)), names(formals(glm.control)))] <- NULL

  res <- do.call("glm", c(list(formula = formula, data = data, family = binomial(link = link)), A))

  if (linear) pred <- predict(res, type = "link")
  else pred <- predict(res, type = "response")

  return(list(model = res, distance = pred))
}
# distance2logit <- function(formula, data, ...) {
#   res <- glm(formula, data, family=binomial(logit), ...)
#   return(list(model = res, distance = fitted(res)))
# }
#
# distance2linear.logit <- function(formula, data, ...) {
#   res <- glm(formula, data, family=binomial(logit), ...)
#   return(list(model = res, distance = predict(res)))
# }
#
# distance2probit <- function(formula, data, ...) {
#   res <- glm(formula, data, family=binomial(probit), ...)
#   return(list(model = res, distance = fitted(res)))
# }
#
# distance2linear.probit <- function(formula, data, ...) {
#   res <- glm(formula, data, family=binomial(probit), ...)
#   return(list(model = res, distance = predict(res)))
# }
#
# distance2cloglog <- function(formula, data, ...) {
#   res <- glm(formula, data, family=binomial(cloglog), ...)
#   return(list(model = res, distance = fitted(res)))
# }
#
# distance2linear.cloglog <- function(formula, data, ...) {
#   res <- glm(formula, data, family=binomial(cloglog), ...)
#   return(list(model = res, distance = predict(res)))
# }
#
# distance2log <- function(formula, data, ...) {
#   res <- glm(formula, data, family=binomial(log), ...)
#   return(list(model = res, distance = fitted(res)))
# }
#
# distance2linear.log <- function(formula, data, ...) {
#   res <- glm(formula, data, family=binomial(log), ...)
#   return(list(model = res, distance = predict(res)))
# }
#
# distance2cauchit <- function(formula, data, ...) {
#   res <- glm(formula, data, family=binomial(cauchit), ...)
#   return(list(model = res, distance = fitted(res)))
# }
#
# distance2linearcauchit <- function(formula, data, ...) {
#   res <- glm(formula, data, family=binomial(cauchit), ...)
#   return(list(model = res, distance = predict(res)))
# }

#distance2gam-----------------
distance2gam <- function(formula, data, link = "logit", ...) {
  check.package("mgcv")

  if (!is.null(link) && startsWith(as.character(link), "linear")) {
    linear <- TRUE
    link <- sub("linear.", "", as.character(link))
  }
  else linear <- FALSE

  res <- mgcv::gam(formula, data, family = binomial(link), ...)

  if (linear) pred <- predict(res, type = "link")
  else pred <- predict(res, type = "response")

  return(list(model = res, distance = pred))
}

# distance2GAMlogit <- function(formula, data, ...) {
#   check.package("mgcv")
#   res <- mgcv::gam(formula, data, family=binomial(logit), ...)
#   return(list(model = res, distance = fitted(res)))
# }
#
# distance2GAMprobit <- function(formula, data, ...) {
#   check.package("mgcv")
#   res <- mgcv::gam(formula, data, family=binomial(probit), ...)
#   return(list(model = res, distance = fitted(res)))
# }
#
# distance2GAMcloglog <- function(formula, data, ...) {
#   check.package("mgcv")
#   res <- mgcv::gam(formula, data, family=binomial(cloglog), ...)
#   return(list(model = res, distance = fitted(res)))
# }
#
#
# distance2GAMlog <- function(formula, data, ...) {
#   check.package("mgcv")
#   res <- mgcv::gam(formula, data, family=binomial(log), ...)
#   return(list(model = res, distance = fitted(res)))
#
# }
#
# distance2GAMcauchit <- function(formula, data, ...) {
#   check.package("mgcv")
#   res <- mgcv::gam(formula, data, family=binomial(cauchit), ...)
#   return(list(model = res, distance = fitted(res)))
# }

#distance2rpart-----------------
distance2rpart <- function(formula, data, link = NULL, ...) {
  check.package("rpart")
  A <- list(...)
  A[!names(A) %in% c(names(formals(rpart::rpart)), names(formals(rpart::rpart.control)))] <- NULL
  A$formula <- formula
  A$data <- data
  A$method <- "class"

  res <- do.call(rpart::rpart, A)
  return(list(model = res, distance = predict(res, type = "prob")[,"1"]))
}

#distance2nnet-----------------
distance2nnet <- function(formula, data, link = NULL, ...) {
  check.package("nnet")
  res <- nnet::nnet(formula, data, entropy = TRUE, ...)
  return(list(model = res, distance = drop(fitted(res))))
}
#distance2cbps-----------------
distance2cbps <- function(formula, data, link = NULL, ...) {
  check.package("CBPS")

  if (!is.null(link) && startsWith(as.character(link), "linear")) {
    linear <- TRUE
  }
  else linear <- FALSE

  A <- list(...)

  A[["standardized"]] <- FALSE
  if (is.null(A[["ATT"]])) {
    if (is.null(A[["estimand"]])) A[["ATT"]] <- 1
    else {
      estimand <- toupper(A[["estimand"]])
      estimand <- match_arg(estimand, c("ATT", "ATC", "ATE"))
      A[["ATT"]] <- switch(estimand, "ATT" = 1, "ATC" = 2, 0)
    }
  }
  if (is.null(A[["method"]])) {
    if (is.null(A[["over"]])) A[["method"]] <- "over"
    else {
      A[["method"]] <- if (isFALSE(A[["over"]])) "exact" else "over"
    }
  }
  A[c("estimand", "over")] <- NULL
  res <- do.call(CBPS::CBPS, c(list(formula, data), A), quote = TRUE)

  pred <- fitted(res)
  if (linear) pred <- qlogis(pred)

  return(list(model = res, distance = pred))
}

#distance2bart----------------
distance2bart <- function(formula, data, link = NULL, ...) {
  check.package("dbarts")

  if (!is.null(link) && startsWith(as.character(link), "linear")) {
    linear <- TRUE
  }
  else linear <- FALSE

  A <- list(...)
  A[!names(A) %in% c(names(formals(dbarts::bart2)), names(formals(dbarts::dbartsControl)))] <- NULL
  A$formula <- formula
  A$data <- data

  res <- do.call(dbarts::bart2, A)

  if (linear) pred <- fitted(res, type = "link")
  else pred <- fitted(res, type = "response")

  return(list(model = res, distance = pred))
}

#distance2randomforest-----------------
distance2randomforest <- function(formula, data, link = NULL, ...) {
  check.package("randomForest")
  newdata <- get_all_vars(formula, data)
  treatvar <- as.character(formula[[2]])
  newdata[[treatvar]] <- factor(newdata[[treatvar]], levels = c("0", "1"))
  res <- randomForest::randomForest(formula, data = newdata, ...)
  return(list(model = res, distance = predict(res, type = "prob")[,"1"]))
}