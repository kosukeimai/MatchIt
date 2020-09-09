#' @export
matchit <- function(formula, data = NULL, method = "nearest", distance = "glm",
                    link = "logit", distance.options = list(), estimand = "ATT",
                    exact = NULL, mahvars = NULL, discard = "none",
                    reestimate = FALSE, replace = FALSE, m.order = NULL,
                    caliper = NULL, std.caliper = TRUE, ratio = 1, verbose = FALSE, ...) {

  #Checking input format
  #data input
  mcall <- match.call()

  #Process formula and data inputs
  if (!inherits(formula, "formula")) stop("'formula' must be a formula object.", call. = FALSE)
  f.env <- environment(formula)

  treat.formula <- update(formula, . ~ 0)
  treat.mf <- model.frame(treat.formula, data = data, na.action = "na.pass")
  treat <- model.response(treat.mf)
  if (anyNA(treat)) stop("Missing values are not allowed in the treatment.", call. = FALSE)
  treat <- binarize(treat) #make 0/1
  names(treat) <- rownames(treat.mf)

  n.obs <- length(treat)

  ## Process method
  if (length(method) == 1 && is.atomic(method)) {
    method <- tolower(as.character(method))
    fn2 <- paste0("matchit2", method)
    if (!exists(fn2)) stop(paste0("\"", method, "\" is not a supported method."), call. = FALSE)
  }
  else {
    stop("method must be the name of a supported matching method. See ?matchit for allowable options.", call. = FALSE)
  }

  #Process distance and discard
  check.inputs(method, distance = distance, mcall = mcall, exact = exact,
               mahvars = mahvars, caliper = caliper, discard = discard,
               reestimate = reestimate, replace = replace, ratio = ratio,
               m.order = m.order, estimand = estimand)

  exactcovs <- mahcovs <- calcovs <- NULL

  if (method %in% c("exact", "cem")) {
    fn1 <- NULL
  }
  else {
    distance <- process.distance(distance, method)

    if (is.numeric(distance)) {
      if (length(distance) != length(treat)) stop("'distance' must be the same length as the dataset if specified as a numeric vector.", call. = FALSE)
      if (anyNA(distance)) stop("Missing values are not allowed in 'distance'.", call. = FALSE)
      fn1 <- "distance2user"
    }
    else {
      fn1 <- paste0("distance2", distance)
    }

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
        stop("'exact' must be supplied as a character vector of names or a one-sided formula.", call. = FALSE)
      }
      mahcovs <- model.frame(mahvars, data, na.action = "na.pass")
    }
  }

  is.full.mahalanobis <- identical(fn1, "distance2mahalanobis")

  #Estimate distance, discard form common support, optionally re-estimate distance

  if (is.null(fn1)) {
    #Only for cem and exact
    dist.model <- distance <- link <- NULL
    discarded <- rep(FALSE, length(treat))
  }
  else {
    if (fn1 == "distance2user") {
      dist.model <- link <- NULL
      discarded <- discard(treat, distance, discard)
    }
    else if (fn1 == "distance2mahalanobis") {
      distance <- link <- dist.model <- NULL
      discarded <- discard(treat, distance, discard)
    }
    else {
      #Estimate distance
      if (is.null(distance.options$formula)) distance.options$formula <- formula
      if (is.null(distance.options$data)) distance.options$data <- data
      if (is.null(distance.options$verbose)) distance.options$verbose <- verbose
      if (is.null(distance.options$estimand)) distance.options$estimand <- estimand
      if (!is.null(attr(distance, "link"))) distance.options$link <- attr(distance, "link")
      else distance.options$link <- link

      dist.out <- do.call(fn1, distance.options, quote = TRUE)

      dist.model <- dist.out$model
      distance <- dist.out$distance

      #Discard
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

      #Remove smoothing terms from gam formula
      if (inherits(dist.model, "gam")) {
        formula <- mgcv::interpret.gam(formula)$fake.formula
      }
    }
  }

  covs.formula <- delete.response(terms(formula))
  covs <- model.frame(covs.formula, data = data, na.action = "na.pass")
  if (anyNA(covs)) stop("Missing values are not allowed in the covariates.", call. = FALSE)
  for (i in which(vapply(covs, is.character, logical(1L)))) covs[[i]] <- factor(covs[[i]])

  caliper <- process.caliper(caliper, method, data, covs, mahcovs, distance, discarded, std.caliper)

  if (!is.null(attr(caliper, "cal.formula"))) {
    calcovs <- model.frame(attr(caliper, "cal.formula"), data, na.action = "na.pass")
    attr(caliper, "cal.formula") <- NULL
  }

  ## matching!
  match.out <- do.call(fn2, list(treat = treat, covs = covs, data = data, distance = distance,
                                 discarded = discarded, exact = exact, mahvars = mahvars,
                                 replace = replace, m.order = m.order, caliper = caliper,
                                 ratio = ratio, is.full.mahalanobis = is.full.mahalanobis,
                                 formula = formula, estimand = estimand, verbose = verbose, ...),
                       quote = TRUE)

  ## Sample size
  nn <- matrix(0, ncol=2, nrow=4, dimnames = list(c("All","Matched","Unmatched","Discarded"),
                                                  c("Control", "Treated")))
  nn["All",] <- c(sum(treat==0), sum(treat==1))
  nn["Matched",] <- c(sum(treat==0 & match.out$weights > 0), sum(treat==1 & match.out$weights > 0))
  nn["Unmatched",] <- c(sum(treat==0 & match.out$weights==0 & !discarded), sum(treat==1 & match.out$weights==0 & !discarded))
  nn["Discarded",] <- c(sum(treat==0 & discarded), sum(treat==1 & discarded))

  info <- list(method = method,
               distance = if (!is.null(fn1)) sub("distance2", "", fn1, fixed = TRUE) else NULL,
               link = if (!is.null(link)) link else NULL,
               discard = discard,
               replace = if (method %in% c("nearest", "genetic")) replace else NULL,
               ratio = if (method %in% c("nearest", "optimal", "genetic")) ratio else NULL,
               mahalanobis = is.full.mahalanobis || !is.null(mahvars),
               subclass = if (method == "subclass") length(unique(match.out$subclass[!is.na(match.out$subclass)])) else NULL)

  #Create X.list for X output, removing duplicate variables
  X.list <- list(covs, exactcovs, mahcovs, calcovs)
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
  match.out$distance <- distance
  match.out$discarded <- discarded
  match.out$exact <- exact
  match.out$mahvars <- mahvars
  match.out$caliper <- caliper
  match.out$nn <- nn

  return(match.out)
}
