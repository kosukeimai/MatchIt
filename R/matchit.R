matchit <- function(formula, data = NULL, method = "nearest", distance = "glm",
                    link = "logit", distance.options = list(), estimand = "ATT",
                    exact = NULL, mahvars = NULL, antiexact = NULL, discard = "none",
                    reestimate = FALSE, s.weights = NULL, replace = FALSE, m.order = NULL,
                    caliper = NULL, std.caliper = TRUE, ratio = 1, verbose = FALSE, ...) {

  #Checking input format
  #data input
  mcall <- match.call()

  ## Process method
  if (length(method) == 1 && is.character(method)) {
    method <- tolower(method)
    method <- match_arg(method, c("exact", "cem", "nearest", "optimal", "full", "genetic", "subclass", "cardinality"))
    fn2 <- paste0("matchit2", method)
  }
  else if (is.null(method)) {
    fn2 <- "matchit2null"
  }
  else {
    stop("'method' must be the name of a supported matching method. See ?matchit for allowable options.", call. = FALSE)
  }

  #Process formula and data inputs
  if (!inherits(formula, "formula")) stop("'formula' must be a formula object.", call. = FALSE)

  treat.form <- update(terms(formula, data = data), . ~ 0)
  treat.mf <- model.frame(treat.form, data = data, na.action = "na.pass")
  treat <- model.response(treat.mf)
  if (anyNA(treat)) stop("Missing values are not allowed in the treatment.", call. = FALSE)
  if (length(unique(treat)) != 2) stop("The treatment must be a binary variable.", call. = FALSE)
  treat <- binarize(treat) #make 0/1
  names(treat) <- rownames(treat.mf)

  n.obs <- length(treat)

  #Process inputs
  ignored.inputs <- check.inputs(mcall = mcall, method = method, distance = distance, exact = exact,
                                 mahvars = mahvars, antiexact = antiexact, caliper = caliper, discard = discard,
                                 reestimate = reestimate, s.weights = s.weights, replace = replace,
                                 ratio = ratio, m.order = m.order, estimand = estimand)

  if (length(ignored.inputs) > 0) {
    for (i in ignored.inputs) assign(i, NULL)
  }

  #Process s.weights
  if (!is.null(s.weights)) {
    if (is.character(s.weights)) {
      if (is.null(data) || !is.data.frame(data)) {
        stop("If 's.weights' is specified a string, a data frame containing the named variable must be supplied to 'data'.", call. = FALSE)
      }
      if (!all(s.weights %in% names(data))) {
        stop("The name supplied to 's.weights' must be a variable in 'data'.", call. = FALSE)
      }
      s.weights.form <- reformulate(s.weights)
      s.weights <- model.frame(s.weights.form, data, na.action = "na.pass")
      if (ncol(s.weights) != 1) stop("'s.weights' can only contain one named variable.", call. = FALSE)
      s.weights <- s.weights[[1]]
    }
    else if (inherits(s.weights, "formula")) {
      s.weights.form <- update(s.weights, NULL ~ .)
      s.weights <- model.frame(s.weights.form, data, na.action = "na.pass")
      if (ncol(s.weights) != 1) stop("'s.weights' can only contain one named variable.", call. = FALSE)
      s.weights <- s.weights[[1]]
    }
    else if (!is.numeric(s.weights)) {
      stop("'s.weights' must be supplied as a numeric vector, string, or one-sided formula.", call. = FALSE)
    }

    if (anyNA(s.weights)) stop("Missing values are not allowed in 's.weights'.", call. = FALSE)
    if (length(s.weights) != n.obs) stop("'s.weights' must be the same length as the treatment vector.", call. = FALSE)

    names(s.weights) <- names(treat)

  }

  #Process distance function
  if (!is.null(method) && method %in% c("exact", "cem", "cardinality")) {
    fn1 <- NULL
  }
  else {
    distance <- process.distance(distance, method, treat)

    if (is.numeric(distance)) {
      fn1 <- "distance2user"
    }
    else {
      fn1 <- paste0("distance2", distance)
    }
  }
  is.full.mahalanobis <- identical(fn1, "distance2mahalanobis")

  #Process exact, mahvars, and antiexact
  exactcovs <- mahcovs <- antiexactcovs <- NULL

  if (!is.null(exact)) {
    if (is.character(exact)) {
      if (is.null(data) || !is.data.frame(data)) {
        stop("If 'exact' is specified as strings, a data frame containing the named variables must be supplied to 'data'.", call. = FALSE)
      }
      if (!all(exact %in% names(data))) {
        stop("All names supplied to 'exact' must be variables in 'data'.", call. = FALSE)
      }
      exact <- reformulate(exact)
    }
    else if (inherits(exact, "formula")) {
      exact <- update(exact, NULL ~ .)
    }
    else {
      stop("'exact' must be supplied as a character vector of names or a one-sided formula.", call. = FALSE)
    }
    exactcovs <- model.frame(exact, data, na.action = "na.pass")
    if (anyNA(exactcovs)) stop("Missing values are not allowed in the covariates named in 'exact'.", call. = FALSE)
  }

  if (!is.null(mahvars)) {
    if (is.character(mahvars)) {
      if (is.null(data) || !is.data.frame(data)) {
        stop("If 'mahvars' is specified as strings, a data frame containing the named variables must be supplied to 'data'.", call. = FALSE)
      }
      if (!all(mahvars %in% names(data))) {
        stop("All names supplied to 'mahvars' must be variables in 'data'.", call. = FALSE)
      }
      mahvars <- reformulate(mahvars)
    }
    else if (inherits(mahvars, "formula")) {
      mahvars <- update(mahvars, NULL ~ .)
    }
    else {
      stop("'mahvars' must be supplied as a character vector of names or a one-sided formula.", call. = FALSE)
    }
    mahcovs <- model.frame(mahvars, data, na.action = "na.pass")
    if (anyNA(mahcovs)) stop("Missing values are not allowed in the covariates named in 'mahvars'.", call. = FALSE)
  }

  if (!is.null(antiexact)) {
    if (is.character(antiexact)) {
      if (is.null(data) || !is.data.frame(data)) {
        stop("If 'antiexact' is specified as strings, a data frame containing the named variables must be supplied to 'data'.", call. = FALSE)
      }
      if (!all(antiexact %in% names(data))) {
        stop("All names supplied to 'antiexact' must be variables in 'data'.", call. = FALSE)
      }
      antiexact <- reformulate(antiexact)
    }
    else if (inherits(antiexact, "formula")) {
      antiexact <- update(antiexact, NULL ~ .)
    }
    else {
      stop("'antiexact' must be supplied as a character vector of names or a one-sided formula.", call. = FALSE)
    }
    antiexactcovs <- model.frame(antiexact, data, na.action = "na.pass")
    if (anyNA(antiexactcovs)) stop("Missing values are not allowed in the covariates named in 'antiexact'.", call. = FALSE)
  }

  #Estimate distance, discard from common support, optionally re-estimate distance
  if (is.null(fn1) || fn1 == "distance2mahalanobis") {
    #No distance measure
    dist.model <- distance <- link <- NULL
  }
  else if (fn1 == "distance2user") {
    dist.model <- link <- NULL
  }
  else {
    if (!is.null(s.weights)) attr(s.weights, "in_ps") <- !distance %in% c("bart", "randomforest")

    #Estimate distance
    if (is.null(distance.options$formula)) distance.options$formula <- formula
    if (is.null(distance.options$data)) distance.options$data <- data
    if (is.null(distance.options$verbose)) distance.options$verbose <- verbose
    if (is.null(distance.options$estimand)) distance.options$estimand <- estimand
    if (is.null(distance.options$weights) && !fn1 %in% c("distance2bart", "distance2randomforest")) {
      distance.options$weights <- s.weights
    }

    if (!is.null(attr(distance, "link"))) distance.options$link <- attr(distance, "link")
    else distance.options$link <- link

    dist.out <- do.call(fn1, distance.options, quote = TRUE)

    dist.model <- dist.out$model
    distance <- dist.out$distance

    #Remove smoothing terms from gam formula
    if (inherits(dist.model, "gam")) {
      formula <- mgcv::interpret.gam(formula)$fake.formula
    }
  }

  #Process covs
  covs.formula <- delete.response(terms(formula, data = data))
  covs <- model.frame(covs.formula, data = data, na.action = "na.pass")
  if (anyNA(covs)) stop("Missing values are not allowed in the covariates.", call. = FALSE)
  for (i in which(vapply(covs, is.character, logical(1L)))) covs[[i]] <- factor(covs[[i]])

  #Process discard
  if (is.null(fn1) || fn1 %in% c("distance2mahalanobis", "distance2user")) {
    discarded <- discard(treat, distance, discard)
  }
  else {
    discarded <- discard(treat, dist.out$distance, discard)

    #Optionally reestimate
    if (reestimate && any(discarded)) {
      for (i in seq_along(distance.options)) {
        if (length(distance.options[[i]]) == n.obs) {
          distance.options[[i]] <- distance.options[[i]][!discarded]
        }
        else if (length(dim(distance.options[[i]])) == 2 && nrow(distance.options[[i]]) == n.obs) {
          distance.options[[i]] <- distance.options[[i]][!discarded,,drop = FALSE]
        }
      }
      dist.out <- do.call(fn1, distance.options, quote = TRUE)
      dist.model <- dist.out$model
      distance[!discarded] <- dist.out$distance
    }
  }

  #Process caliper
  calcovs <- NULL

  if (!is.null(caliper)) {
    caliper <- process.caliper(caliper, method, data, covs, mahcovs, distance, discarded, std.caliper)

    if (!is.null(attr(caliper, "cal.formula"))) {
      calcovs <- model.frame(attr(caliper, "cal.formula"), data, na.action = "na.pass")
      if (anyNA(calcovs)) stop("Missing values are not allowed in the covariates named in 'caliper'.", call. = FALSE)
      attr(caliper, "cal.formula") <- NULL
    }
  }

  #Matching!
  match.out <- do.call(fn2, list(treat = treat, covs = covs, data = data, distance = distance,
                                 discarded = discarded, exact = exact, mahvars = mahvars,
                                 replace = replace, m.order = m.order, caliper = caliper,
                                 s.weights = s.weights, ratio = ratio, is.full.mahalanobis = is.full.mahalanobis,
                                 formula = formula, estimand = estimand, verbose = verbose,
                                 antiexact = antiexact, ...),
                       quote = TRUE)

  info <- create_info(method, fn1, link, discard, replace, ratio, mcall,
                      mahalanobis = is.full.mahalanobis || !is.null(mahvars),
                      subclass = match.out$subclass,
                      antiexact = colnames(antiexactcovs),
                      distance_is_matrix = !is.null(distance) && is.matrix(distance))

  #Create X.list for X output, removing duplicate variables
  X.list <- list(covs, exactcovs, mahcovs, calcovs, antiexactcovs)
  all.covs <- lapply(X.list, names)
  for (i in seq_along(X.list)[-1]) if (!is.null(X.list[[i]])) X.list[[i]][names(X.list[[i]]) %in% unlist(all.covs[1:(i-1)])] <- NULL
  X.list[lengths(X.list) == 0] <- NULL

  ## putting all the results together
  match.out$model <- dist.model
  match.out$X <- do.call("cbind", X.list)
  match.out$call <- mcall
  match.out$info <- info
  match.out$estimand <- estimand
  match.out$formula <- formula
  match.out$treat <- treat
  match.out$distance <- if (!is.null(distance) && !is.matrix(distance)) setNames(distance, names(treat))
  match.out$discarded <- discarded
  match.out$s.weights <- s.weights
  match.out$exact <- exact
  match.out$mahvars <- mahvars
  match.out$caliper <- caliper
  match.out$nn <- nn(treat, match.out$weights, discarded, s.weights)

  return(match.out)
}
