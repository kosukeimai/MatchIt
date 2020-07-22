#' @export
matchit <- function(formula, data = NULL, method = "nearest", distance = "glm",
                    link = "logit", distance.options=list(), exact = NULL, mahvars = NULL,
                    discard = "none", replace = FALSE, m.order = NULL, caliper = NULL,
                    ratio = 1, reestimate = FALSE, verbose = FALSE, ...) {

  #Checking input format
  #data input
  mcall <- match.call()

  #Process formula and data inputs
  if (!inherits(formula, "formula")) stop("formula must be a formula object.", call. = FALSE)
  f.env <- environment(formula)

  data0 <- get_all_vars(formula, data = data)
  for (i in names(data0)[vapply(data0, is.character, logical(1L))]) data0[[i]] <- factor(data0[[i]])

  treat.formula <- update(formula, . ~ 0)
  treat.mf <- model.frame(treat.formula, data = data0, na.action = "na.pass")
  treat <- model.response(treat.mf)
  if (anyNA(treat)) stop("Missing values are not allowed in the treatment.", call. = FALSE)
  names(treat) <- rownames(data0)

  covs.formula <- delete.response(terms(formula))
  covs <- get_all_vars(covs.formula, data = data0)
  if (anyNA(covs)) stop("Missing values are not allowed in the covariates.", call. = FALSE)

  n.obs <- length(treat)

  ## Process method
  if (length(method) == 1 && is.atomic(method)) {
    method <- as.character(method)
    fn2 <- paste0("matchit2", method)
    if (!exists(fn2)) stop(paste0("\"", method, "\" is not a supported method."), call. = FALSE)
  }
  else {
    stop("method must be the name of a supported matching method. See ?matchit for allowable options.", call. = FALSE)
  }

  #Process distance and discard
  check.inputs(method, distance, mcall, exact, mahvars, caliper, discard, replace, ratio, m.order)

  exactcovs <- mahcovs <- NULL

  if (method %in% c("exact", "cem")) {
    distance <- fn1 <- NULL
  }
  else {
    distance <- process.distance(distance, method)

    if (is.numeric(distance)) {
      if (length(distance) != length(treat)) stop("distance must be the same length as the dataset.", call. = FALSE)
      if (anyNA(distance)) stop("Missing values are not allowed in distance.", call. = FALSE)
      fn1 <- "distance2user"
    }
    else {
      fn1 <- paste0("distance2", distance)
    }

    if (!is.null(exact)) {
      if (is.character(exact)) {
        if (is.null(data) || !is.data.frame(data)) {
          stop("If exact is specified as strings, a data frame containing the named variables must be supplied to data.", call. = FALSE)
        }
        if (!all(exact %in% names(data))) {
          stop("All names supplied to exact must be variables in data.", call. = FALSE)
        }
        exact <- reformulate(exact)
      }
      else if (inherits(exact, "formula")) {
        exact <- update(exact, NULL ~ .)
      }
      else {
        stop("exact must be supplied as a character vector of names or a one-sided formula.", call. = FALSE)
      }
      exactcovs <- get_all_vars(exact, data)
    }

    if (!is.null(mahvars)) {
      if (is.character(mahvars)) {
        if (is.null(data) || !is.data.frame(data)) {
          stop("If mahvars is specified as strings, a data frame containing the named variables must be supplied to data.", call. = FALSE)
        }
        if (!all(mahvars %in% names(data))) {
          stop("All names supplied to mahvars must be variables in data.", call. = FALSE)
        }
        mahvars <- reformulate(mahvars)
      }
      else if (inherits(mahvars, "formula")) {
        mahvars <- update(mahvars, NULL ~ .)
      }
      else {
        stop("exact must be supplied as a character vector of names or a one-sided formula.", call. = FALSE)
      }
      mahcovs <- get_all_vars(mahvars, data)
    }
  }

  is.full.mahalanobis <- identical(fn1, "distance2mahalanobis")

  #Estimate distance, discard form common support, optionally re-estimate distance

  if (is.null(fn1)) {
    #Only for cem and exact
    dist.model <- NULL
    discarded <- rep(FALSE, length(treat))
  }
  else if (fn1 == "distance2user") {
    dist.model <- NULL
    discarded <- discard(treat, distance, discard, covs)
  }
  else if (fn1 == "distance2mahalanobis") {
    distance <- NULL
    dist.model <- NULL
    discarded <- discard(treat, distance, discard, covs)
  }
  else {
    #Estimate distance
    if (is.null(distance.options$formula)) distance.options$formula <- formula
    if (is.null(distance.options$data)) distance.options$data <- data
    if (is.null(distance.options$verbose)) distance.options$verbose <- verbose
    if (!is.null(attr(distance, "link"))) distance.options$link <- attr(distance, "link")
    else distance.options$link <- link

    dist.out <- do.call(fn1, distance.options, quote = TRUE)

    dist.model <- dist.out$model
    distance <- dist.out$distance

    #Discard
    discarded <- discard(treat, dist.out$distance, discard, covs)

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


  # ## obtain T and X
  # tryerror <- try(model.frame(formula), TRUE)
  # if (is.character(distance) &&
  #       distance %in% c("GAMlogit", "GAMprobit", "GAMcloglog", "GAMlog", "GAMcauchit")) {
  #   requireNamespace("mgcv")
  #   tt <- terms(mgcv::interpret.gam(formula)$fake.formula)
  # } else {
  #   tt <- terms(formula)
  # }
  # attr(tt, "intercept") <- 0
  # mf <- model.frame(tt, data)
  # treat <- model.response(mf)
  # X <- model.matrix(tt, data=mf)

  ## matching!
  match.out <- do.call(fn2, list(treat = treat, covs = covs, data = data, distance = distance,
                                 discarded = discarded, exact = exact, mahvars = mahvars,
                                 replace = replace, m.order = m.order, caliper = caliper,
                                 ratio = ratio, is.full.mahalanobis = is.full.mahalanobis,
                                 formula = formula, ...), quote = TRUE)

  ## basic summary
  nn <- matrix(0, ncol=2, nrow=4)
  nn[1,] <- c(sum(treat==0), sum(treat==1))
  nn[2,] <- c(sum(treat==0 & match.out$weights > 0), sum(treat==1 & match.out$weights > 0))
  nn[3,] <- c(sum(treat==0 & match.out$weights==0 & !discarded), sum(treat==1 & match.out$weights==0 & !discarded))
  nn[4,] <- c(sum(treat==0 & discarded), sum(treat==1 & discarded))
  dimnames(nn) <- list(c("All","Matched","Unmatched","Discarded"),
                       c("Control", "Treated"))

  #Create X.list for X output, removing duplicate variables
  X.list <- list(covs, exactcovs, mahcovs)
  all.covs <- lapply(X.list, names)
  for (i in seq_along(X.list)[-1]) if (!is.null(X.list[[i]])) X.list[[i]][names(X.list[[i]]) %in% unlist(all.covs[1:(i-1)])] <- NULL
  X.list[lengths(X.list) == 0] <- NULL

  ## putting all the results together
  match.out$model <- dist.model
  match.out$X <- do.call("data.frame", X.list)
  match.out$call <- mcall
  match.out$formula <- formula
  match.out$treat <- treat
  match.out$distance <- distance
  match.out$discarded <- discarded
  match.out$exact <- exact
  match.out$mahvars <- mahvars
  match.out$nn <- nn

  return(match.out)
}
