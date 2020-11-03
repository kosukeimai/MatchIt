#distance2glm-----------------
distance2glm <- function(formula, data, link = "logit", ...) {

  if (!is.null(link) && startsWith(as.character(link), "linear")) {
    linear <- TRUE
    link <- sub("linear.", "", as.character(link), fixed = TRUE)
  }
  else linear <- FALSE

  A <- list(...)
  A[!names(A) %in% c(names(formals(glm)), names(formals(glm.control)))] <- NULL

  res <- do.call("glm", c(list(formula = formula, data = data, family = quasibinomial(link = link)), A))

  if (linear) pred <- predict(res, type = "link")
  else pred <- predict(res, type = "response")

  return(list(model = res, distance = pred))
}

#distance2gam-----------------
distance2gam <- function(formula, data, link = "logit", ...) {
  check.package("mgcv")

  if (!is.null(link) && startsWith(as.character(link), "linear")) {
    linear <- TRUE
    link <- sub("linear.", "", as.character(link), fixed = TRUE)
  }
  else linear <- FALSE

  A <- list(...)
  weights <- A$weights
  A$weights <- NULL

  res <- do.call(mgcv::gam, c(list(formula, data, family = quasibinomial(link),
                                   weights = weights), A),
                 quote = TRUE)

  if (linear) pred <- predict(res, type = "link")
  else pred <- predict(res, type = "response")

  return(list(model = res, distance = pred))
}

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

  A <- list(...)
  weights <- A$weights
  A$weights <- NULL

  res <- do.call(nnet::nnet, c(list(formula, data, entropy = TRUE), A), quote = TRUE)
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

  if (!is.null(A[["weights"]])) {
    A[["sample.weights"]] <- A[["weights"]]
    A[["weights"]] <- NULL
  }

  capture.output({ #Keeps from printing message about treatment
    res <- do.call(CBPS::CBPS, c(list(formula, data), A), quote = TRUE)
  })

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

# distance2bart <- function(formula, data, link = NULL, ...) {
#   check.package("BART")
#
#   if (!is.null(link) && startsWith(as.character(link), "linear")) {
#     linear <- TRUE
#     link <- sub("linear.", "", as.character(link), fixed = TRUE)
#   }
#   else linear <- FALSE
#
#   #Keep link probit because default in matchit is logit but probit is much faster with BART
#   link <- "probit"
#
#   # if (is.null(link)) link <- "probit"
#   # else if (!link %in% c("probit", "logit")) {
#   #   stop("'link' must be \"probit\" or \"logit\" with distance = \"bart\".", call. = FALSE)
#   # }
#
#   data <- model.frame(formula, data)
#
#   treat <- binarize(data[[1]])
#   X <- data[-1]
#
#   chars <- vapply(X, is.character, logical(1L))
#   X[chars] <- lapply(X[chars], factor)
#
#   A <- list(...)
#
#   if (!is.null(A[["mc.cores"]]) && A[["mc.cores"]][1] > 1) fun <- BART::mc.gbart
#   else fun <- BART::gbart
#
#   res <- do.call(fun, c(list(X,
#                              y.train = treat,
#                              type = switch(link, "logit" = "lbart", "pbart")),
#                         A[intersect(names(A), setdiff(names(formals(fun)),
#                                                       c("x.train", "y.train", "x.test", "type", "ntype")))]))
#
#   pred <- res$prob.train.mean
#   if (linear) pred <- switch(link, logit = qlogis, probit = qnorm)(pred)
#
#   return(list(model = res, distance = pred))
# }

#distance2randomforest-----------------
distance2randomforest <- function(formula, data, link = NULL, ...) {
  check.package("randomForest")
  newdata <- get_all_vars(formula, data)
  treatvar <- as.character(formula[[2]])
  newdata[[treatvar]] <- factor(newdata[[treatvar]], levels = c("0", "1"))
  res <- randomForest::randomForest(formula, data = newdata, ...)
  return(list(model = res, distance = predict(res, type = "prob")[,"1"]))
}