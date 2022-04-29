#Functions to compute distance matrices
eucdist_internal <- function(X, treat = NULL) {

  if (is.null(dim(X))) X <- as.matrix(X)

  if (is.null(treat)) {
    if (ncol(X) == 1) {
      d <- abs(outer(drop(X), drop(X), "-"))
    }
    else {
      d <- dist_to_matrixC(dist(X))
    }
    dimnames(d) <- list(rownames(X), rownames(X))
  }
  else {
    treat_l <- as.logical(treat)
    X1 <- X[treat_l,, drop = FALSE]
    X0 <- X[!treat_l,, drop = FALSE]

    if (ncol(X) == 1) {
      d <- abs(outer(drop(X1), drop(X0), "-"))
    }
    else {
      d <- dist_to_matrixC(dist(X))[treat == 1, treat == 0, drop = FALSE]
    }
    dimnames(d) <- list(rownames(X1), rownames(X0))
  }

  d
}

mahalanobis_dist <- function(formula = NULL, data = NULL, treat = NULL, s.weights = NULL, var = NULL, ...) {
  X <- transform_covariates(formula, data, method = "mahalanobis",
                            s.weights = s.weights, var = var,
                            treat = treat, ...)
  eucdist_internal(X, attr(X, "treat"))
}

scaled_euclidean_dist <- function(formula = NULL, data = NULL, treat = NULL, s.weights = NULL, ...) {
  X <- transform_covariates(formula, data, method = "scaled_euclidean",
                            s.weights = s.weights, treat = treat, ...)
  eucdist_internal(X, attr(X, "treat"))
}

robust_mahalanobis_dist <- function(formula = NULL, data = NULL, treat = NULL, s.weights = NULL, ...) {
  X <- transform_covariates(formula, data = data, method = "robust_mahalanobis",
                            s.weights = s.weights, treat = treat, ...)
  eucdist_internal(X, attr(X, "treat"))
}

euclidean_dist <- function(formula = NULL, data = NULL, treat = NULL, ...) {
  X <- transform_covariates(formula, data = data, method = "euclidean",
                            treat = treat, ...)
  eucdist_internal(X, attr(X, "treat"))
}

#Transforms covariates so that Euclidean distance computed on transforms covariates is equivalent to
#requested distance. When discarded is not NULL, statistics relevant to transformation are computed
#using retained units, but full covariate matrix is returned.
transform_covariates <- function(formula = NULL, data = NULL, method = "mahalanobis", s.weights = NULL, var = NULL, treat = NULL,
                                 discarded = NULL) {
  X <- get.covs.matrix.for.dist(formula, data)

  X <- check_X(X)
  treat <- check_treat(treat, X)

  method <- match_arg(method, matchit_distances())

  if (is.null(discarded)) discarded <- rep(FALSE, nrow(X))

  if (method == "mahalanobis") {
    # X <- sweep(X, 2, colMeans(X))

    if (is.null(var)) {
      X <- scale(X)
      #NOTE: optmatch and Rubin (1980) use pooled within-group covariance matrix
      if (!is.null(treat)) {
        var <- pooled_cov(X[!discarded,, drop = FALSE], treat[!discarded], s.weights[!discarded])
      }
      else if (is.null(s.weights)) {
        var <- cov(X[!discarded,, drop = FALSE])
      }
      else {
        var <- cov.wt(X[!discarded,, drop = FALSE], s.weights[!discarded])$cov
      }
    }
    else if (!is.cov_like(var)) {
      stop("bad var")
    }

    inv_var <- try(solve(var), silent = TRUE)
    if (inherits(inv_var, "try-error")) {
      inv_var <- generalized_inverse(var)
    }

    X <- mahalanobize(X, inv_var)
  }
  else if (method == "robust_mahalanobis") {
    #Rosenbaum (2010, ch8)
    X_r <- matrix(0, nrow = sum(!discarded), ncol = ncol(X),
                  dimnames = list(rownames(X)[!discarded], colnames(X)))
    for (i in seq_len(ncol(X_r))) X_r[,i] <- rank(X[,i])

    if (is.null(s.weights)) var_r <- cov(X_r[!discarded,, drop = FALSE])
    else var_r <- cov.wt(X_r[!discarded,, drop = FALSE], s.weights[!discarded])$cov

    multiplier <- sd(seq_len(sum(!discarded)))/sqrt(diag(var_r))
    var_r <- var_r * outer(multiplier, multiplier, "*")

    inv_var <- try(solve(var_r), silent = TRUE)
    if (inherits(inv_var, "try-error")) {
      inv_var <- generalized_inverse(var_r)
    }

    if (any(discarded)) {
      X_r <- array(0, dim = dim(X), dimnames = dimnames(X))
      for (i in seq_len(ncol(X_r))) X_r[,i] <- rank(X[,i])
    }

    X <- mahalanobize(X_r, inv_var)
  }
  else if (method == "euclidean") {
    #Do nothing
  }
  else if (method == "scaled_euclidean") {
    if (!is.null(treat)) {
      sds <- pooled_sd(X[!discarded,, drop = FALSE], treat[!discarded], s.weights[!discarded])
    }
    else {
      sds <- sqrt(apply(X[!discarded,, drop = FALSE], 2, wvar, w = s.weights))
    }

    for (i in seq_len(ncol(X))) {
      X[,i] <- X[,i]/sds[i]
    }
  }

  attr(X, "treat") <- treat
  X
}

#Get covariates (RHS) vars from formula; factor variable contrasts divided by sqrt(2)
#to ensure same result as when non-factor binary variable supplied (see optmatch:::contr.match_on)
get.covs.matrix.for.dist <- function(formula = NULL, data = NULL) {

  if (is.null(formula)) {
    fnames <- colnames(data)
    fnames[!startsWith(fnames, "`")] <- paste0("`", fnames[!startsWith(fnames, "`")], "`")
    data <- as.data.frame(data)
    formula <- terms(reformulate(fnames), data = data)
  }
  else {
    data <- as.data.frame(data)
    formula <- terms(formula, data = data)
  }

  mf <- model.frame(formula, data, na.action = na.pass)

  chars.in.mf <- vapply(mf, is.character, logical(1L))
  mf[chars.in.mf] <- lapply(mf[chars.in.mf], factor)

  formula <- update(formula, NULL ~ . + 1)

  X <- model.matrix(formula, data = mf,
                    contrasts.arg = lapply(Filter(is.factor, mf),
                                           function(x) contrasts(x, contrasts = FALSE)/sqrt(2)))
  assign <- attr(X, "assign")[-1]
  X <- X[,-1, drop=FALSE]
  attr(X, "assign") <- assign

  attr(X, "treat") <-  model.response(mf)

  return(X)
}
check_X <- function(X) {
  if (isTRUE(attr(X, "checked"))) return(X)

  if (is.data.frame(X)) X <- as.matrix(X)
  else if (is.numeric(X) && is.null(dim(X))) {
    X <- matrix(X, nrow = length(X),
                dimnames = list(names(X), NULL))
  }

  if (!is.numeric(X) || length(dim(X)) != 2) {
    stop("bad X")
  }
  attr(X, "checked") <- TRUE
  X
}
check_treat <- function(treat = NULL, X) {

  if (is.null(treat)) {
    if (is.null(attr(X, "treat"))) return(NULL)
    treat <- attr(X, "treat")
  }
  if (isTRUE(attr(treat, "checked"))) return(treat)

  if (!is.atomic(treat) || !is.null(dim(treat)) || length(treat) != nrow(X)) {
    stop("bad treat")
  }
  treat <- binarize(treat) #make 0/1
  attr(treat, "checked") <- TRUE
  treat
}
is.cov_like <- function(var, X) {
  is.numeric(var) &&
    length(dim(var)) == 2 &&
    (missing(X) || all(dim(var) == ncol(X))) &&
    isSymmetric(var) &&
    all(diag(var) >= 0)
}
matchit_distances <- function() {
  c("mahalanobis", "robust_mahalanobis", "euclidean", "scaled_euclidean")
}
mahalanobize <- function(X, inv_var) {
  ## Mahalanobize covariates by computing cholesky decomp,
  ## allowing for NPD cov matrix by pivoting
  ch <- suppressWarnings(chol(inv_var, pivot = TRUE))
  p <- order(attr(ch, "pivot"))
  # r <- seq_len(attr(ch, "rank"))

  tcrossprod(X, ch[, p, drop = FALSE])
}