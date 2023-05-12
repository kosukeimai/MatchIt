#' Compute a Distance Matrix
#' @name mahalanobis_dist
#'
#' @description
#' The functions compute a distance matrix, either for a single dataset (i.e.,
#' the distances between all pairs of units) or for two groups defined by a
#' splitting variable (i.e., the distances between all units in one group and
#' all units in the other). These distance matrices include the Mahalanobis
#' distance, Euclidean distance, scaled Euclidean distance, and robust
#' (rank-based) Mahalanobis distance. These functions can be used as inputs to
#' the `distance` argument to [matchit()] and are used to compute the
#' corresponding distance matrices within `matchit()` when named.
#'
#' @aliases euclidean_dist scaled_euclidean_dist mahalanobis_dist
#' robust_mahalanobis_dist
#'
#' @param formula a formula with the treatment (i.e., splitting variable) on
#' the left side and the covariates used to compute the distance matrix on the
#' right side. If there is no left-hand-side variable, the distances will be
#' computed between all pairs of units. If `NULL`, all the variables in
#' `data` will be used as covariates.
#' @param data a data frame containing the variables named in `formula`.
#' If `formula` is `NULL`, all variables in `data` will be used
#' as covariates.
#' @param s.weights when `var = NULL`, an optional vector of sampling
#' weights used to compute the variances used in the Mahalanobis, scaled
#' Euclidean, and robust Mahalanobis distances.
#' @param var for `mahalanobis_dist()`, a covariance matrix used to scale
#' the covariates. For `scaled_euclidean_dist()`, either a covariance
#' matrix (from which only the diagonal elements will be used) or a vector of
#' variances used to scale the covariates. If `NULL`, these values will be
#' calculated using formulas described in Details.
#' @param discarded a `logical` vector denoting which units are to be
#' discarded or not. This is used only when `var = NULL`. The scaling
#' factors will be computed only using the non-discarded units, but the
#' distance matrix will be computed for all units (discarded and
#' non-discarded).
#' @param \dots ignored. Included to make cycling through these functions
#' easier without having to change the arguments supplied.
#'
#' @return A numeric distance matrix. When `formula` has a left-hand-side
#' (treatment) variable, the matrix will have one row for each treated unit and
#' one column for each control unit. Otherwise, the matrix will have one row
#' and one column for each unit.
#'
#' @details
#' The **Euclidean distance** (computed using `euclidean_dist()`) is
#' the raw distance between units, computed as \deqn{d_{ij} = \sqrt{(x_i -
#' x_j)(x_i - x_j)'}} where \eqn{x_i} and \eqn{x_j} are vectors of covariates
#' for units \eqn{i} and \eqn{j}, respectively. The Euclidean distance is
#' sensitive to the scales of the variables and their redundancy (i.e.,
#' correlation). It should probably not be used for matching unless all of the
#' variables have been previously scaled appropriately or are already on the
#' same scale. It forms the basis of the other distance measures.
#'
#' The **scaled Euclidean distance** (computed using
#' `scaled_euclidean_dist()`) is the Euclidean distance computed on the
#' scaled covariates. Typically the covariates are scaled by dividing by their
#' standard deviations, but any scaling factor can be supplied using the
#' `var` argument. This leads to a distance measure computed as
#' \deqn{d_{ij} = \sqrt{(x_i - x_j)S_d^{-1}(x_i - x_j)'}} where \eqn{S_d} is a
#' diagonal matrix with the squared scaling factors on the diagonal. Although
#' this measure is not sensitive to the scales of the variables (because they
#' are all placed on the same scale), it is still sensitive to redundancy among
#' the variables. For example, if 5 variables measure approximately the same
#' construct (i.e., are highly correlated) and 1 variable measures another
#' construct, the first construct will have 5 times as much influence on the
#' distance between units as the second construct. The Mahalanobis distance
#' attempts to address this issue.
#'
#' The **Mahalanobis distance** (computed using `mahalanobis_dist()`)
#' is computed as \deqn{d_{ij} = \sqrt{(x_i - x_j)S^{-1}(x_i - x_j)'}} where
#' \eqn{S} is a scaling matrix, typically the covariance matrix of the
#' covariates. It is essentially equivalent to the Euclidean distance computed
#' on the scaled principal components of the covariates. This is the most
#' popular distance matrix for matching because it is not sensitive to the
#' scale of the covariates and accounts for redundancy between them. The
#' scaling matrix can also be supplied using the `var` argument.
#'
#' The Mahalanobis distance can be sensitive to outliers and long-tailed or
#' otherwise non-normally distributed covariates and may not perform well with
#' categorical variables due to prioritizing rare categories over common ones.
#' One solution is the rank-based **robust Mahalanobis distance**
#' (computed using `robust_mahalanobis_dist()`), which is computed by
#' first replacing the covariates with their ranks (using average ranks for
#' ties) and rescaling each ranked covariate by a constant scaling factor
#' before computing the usual Mahalanobis distance on the rescaled ranks.
#'
#' The Mahalanobis distance and its robust variant are computed internally by
#' transforming the covariates in such a way that the Euclidean distance
#' computed on the scaled covariates is equal to the requested distance. For
#' the Mahalanobis distance, this involves replacing the covariates vector
#' \eqn{x_i} with \eqn{x_iS^{-.5}}, where \eqn{S^{-.5}} is the Cholesky
#' decomposition of the (generalized) inverse of the covariance matrix \eqn{S}.
#'
#' When a left-hand-side splitting variable is present in `formula` and
#' `var = NULL` (i.e., so that the scaling matrix is computed internally),
#' the covariance matrix used is the "pooled" covariance matrix, which
#' essentially is a weighted average of the covariance matrices computed
#' separately within each level of the splitting variable to capture
#' within-group variation and reduce sensitivity to covariate imbalance. This
#' is also true of the scaling factors used in the scaled Euclidean distance.
#'
#'
#' @author Noah Greifer
#' @seealso [`distance`], [matchit()], [dist()] (which is used
#' internally to compute Euclidean distances)
#'
#' \pkgfun{optmatch}{match_on}, which provides similar functionality but with fewer
#' options and a focus on efficient storage of the output.
#'
#' @references
#'
#' Rosenbaum, P. R. (2010). *Design of observational studies*.
#' Springer.
#'
#' Rosenbaum, P. R., & Rubin, D. B. (1985). Constructing a Control Group Using
#' Multivariate Matched Sampling Methods That Incorporate the Propensity Score.
#' *The American Statistician*, 39(1), 33–38. \doi{10.2307/2683903}
#'
#' Rubin, D. B. (1980). Bias Reduction Using Mahalanobis-Metric Matching.
#' *Biometrics*, 36(2), 293–298. \doi{10.2307/2529981}
#'
#' @examples
#'
#' data("lalonde")
#'
#' # Computing the scaled Euclidean distance between all units:
#' d <- scaled_euclidean_dist(~ age + educ + race + married,
#'                            data = lalonde)
#'
#' # Another interface using the data argument:
#' dat <- subset(lalonde, select = c(age, educ, race, married))
#' d <- scaled_euclidean_dist(data = dat)
#'
#' # Computing the Mahalanobis distance between treated and
#' # control units:
#' d <- mahalanobis_dist(treat ~ age + educ + race + married,
#'                       data = lalonde)
#'
#' # Supplying a covariance matrix or vector of variances (note:
#' # a bit more complicated with factor variables)
#' dat <- subset(lalonde, select = c(age, educ, married, re74))
#' vars <- sapply(dat, var)
#'
#' d <- scaled_euclidean_dist(data = dat, var = vars)
#'
#' # Same result:
#' d <- scaled_euclidean_dist(data = dat, var = diag(vars))
#'
#' # Discard units:
#' discard <- sample(c(TRUE, FALSE), nrow(lalonde),
#'                   replace = TRUE, prob = c(.2, .8))
#'
#' d <- mahalanobis_dist(treat ~ age + educ + race + married,
#'                       data = lalonde, discarded = discard)
#' dim(d) #all units present in distance matrix
#' table(lalonde$treat)
#'

#Functions to compute distance matrices
#' @export
mahalanobis_dist <- function(formula = NULL,
                             data = NULL,
                             s.weights = NULL,
                             var = NULL,
                             discarded = NULL,
                             ...) {
  X <- transform_covariates(formula, data, method = "mahalanobis",
                            s.weights = s.weights, var = var,
                            discarded = discarded)
  eucdist_internal(X, attr(X, "treat"))
}

#' @export
#' @rdname mahalanobis_dist
scaled_euclidean_dist <- function(formula = NULL,
                                  data = NULL,
                                  s.weights = NULL,
                                  var = NULL,
                                  discarded = NULL,
                                  ...) {
  X <- transform_covariates(formula, data, method = "scaled_euclidean",
                            s.weights = s.weights, var = var,
                            discarded = discarded)
  eucdist_internal(X, attr(X, "treat"))
}

#' @export
#' @rdname mahalanobis_dist
robust_mahalanobis_dist <- function(formula = NULL,
                                    data = NULL,
                                    s.weights = NULL,
                                    discarded = NULL,
                                    ...) {
  X <- transform_covariates(formula, data = data, method = "robust_mahalanobis",
                            s.weights = s.weights, discarded = discarded)
  eucdist_internal(X, attr(X, "treat"))
}

#' @export
#' @rdname mahalanobis_dist
euclidean_dist <- function(formula = NULL,
                           data = NULL,
                           ...) {
  X <- transform_covariates(formula, data = data, method = "euclidean")
  eucdist_internal(X, attr(X, "treat"))
}

#Transforms covariates so that Euclidean distance computed on transforms covariates is equivalent to
#requested distance. When discarded is not NULL, statistics relevant to transformation are computed
#using retained units, but full covariate matrix is returned.
transform_covariates <- function(formula = NULL, data = NULL, method = "mahalanobis",
                                 s.weights = NULL, var = NULL, treat = NULL,
                                 discarded = NULL) {
  X <- get.covs.matrix.for.dist(formula, data)

  X <- .check_X(X)
  treat <- check_treat(treat, X)

  #If allvariables have no variance, use Euclidean to avoid errors
  #If some have no variance, removes those to avoid messing up distances
  no_variance <- which(apply(X, 2, function(x) abs(max(x) - min(x)) < sqrt(.Machine$double.eps)))
  if (length(no_variance) == ncol(X)) {
    method <- "euclidean"
    X <- X[, 1, drop = FALSE]
  }
  else if (length(no_variance) > 0) {
    X <- X[, -no_variance, drop = FALSE]
  }

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
      .err("if `var` is not `NULL`, it must be a covariance matrix with as many entries as supplied variables")
    }

    inv_var <- NULL
    d <- det(var)
    if (d > 1e-8) {
      inv_var <- try(solve(var), silent = TRUE)
    }

    if (d <= 1e-8 || inherits(inv_var, "try-error")) {
      inv_var <- generalized_inverse(var)
    }

    X <- mahalanobize(X, inv_var)
  }
  else if (method == "robust_mahalanobis") {
    #Rosenbaum (2010, ch8)
    X_r <- matrix(0, nrow = sum(!discarded), ncol = ncol(X),
                  dimnames = list(rownames(X)[!discarded], colnames(X)))
    for (i in seq_len(ncol(X_r))) X_r[,i] <- rank(X[!discarded, i])

    if (is.null(s.weights)) var_r <- cov(X_r)
    else var_r <- cov.wt(X_r, s.weights[!discarded])$cov

    multiplier <- sd(seq_len(sum(!discarded)))/sqrt(diag(var_r))
    var_r <- var_r * outer(multiplier, multiplier, "*")

    inv_var <- NULL
    d <- det(var_r)
    if (d > 1e-8) {
      inv_var <- try(solve(var_r), silent = TRUE)
    }

    if (d <= 1e-8 || inherits(inv_var, "try-error")) {
      inv_var <- generalized_inverse(var_r)
    }

    if (any(discarded)) {
      X_r <- array(0, dim = dim(X), dimnames = dimnames(X))
      for (i in seq_len(ncol(X_r))) X_r[!discarded,i] <- rank(X[!discarded,i])
    }

    X <- mahalanobize(X_r, inv_var)
  }
  else if (method == "euclidean") {
    #Do nothing
  }
  else if (method == "scaled_euclidean") {
    if (is.null(var)) {
      if (!is.null(treat)) {
        sds <- pooled_sd(X[!discarded,, drop = FALSE], treat[!discarded], s.weights[!discarded])
      }
      else {
        sds <- sqrt(apply(X[!discarded,, drop = FALSE], 2, wvar, w = s.weights))
      }
    }
    else if (is.cov_like(var, X)) {
      sds <- sqrt(diag(var))
    }
    else if (is.numeric(var) && is.cov_like(diag(var), X)) {
      sds <- sqrt(var)
    }
    else {
      .err("if `var` is not `NULL`, it must be a covariance matrix or a vector of variances with as many entries as supplied variables")
    }

    for (i in seq_len(ncol(X))) {
      X[,i] <- X[,i]/sds[i]
    }
  }

  attr(X, "treat") <- treat
  X
}

#Internal function for fast(ish) Euclidean distance
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
    if (ncol(X) == 1) {
      d <- abs(outer(X[treat_l,], X[!treat_l,], "-"))
    }
    else {
      d <- dist_to_matrixC(dist(X))[treat_l, !treat_l, drop = FALSE]
    }
    dimnames(d) <- list(rownames(X)[treat_l], rownames(X)[!treat_l])
  }

  d
}

#Get covariates (RHS) vars from formula; factor variable contrasts divided by sqrt(2)
#to ensure same result as when non-factor binary variable supplied (see optmatch:::contr.match_on)
get.covs.matrix.for.dist <- function(formula = NULL, data = NULL) {

  if (is.null(formula)) {
    if (is.null(colnames(data))) colnames(data) <- paste0("X", seq_len(ncol(data)))
    fnames <- colnames(data)
    fnames[!startsWith(fnames, "`")] <- add_quotes(fnames[!startsWith(fnames, "`")], "`")
    data <- as.data.frame(data)
    formula <- reformulate(fnames)
  }
  else {
    data <- as.data.frame(data)
  }

  formula <- terms(formula, data = data)

  if (rlang::is_formula(formula, lhs = FALSE)) {
    formula <- update(formula, ~ . + 1)
  }
  else {
    formula <- update(formula, . ~ . + 1)
  }

  mf <- model.frame(formula, data, na.action = na.pass)

  chars.in.mf <- vapply(mf, is.character, logical(1L))
  mf[chars.in.mf] <- lapply(mf[chars.in.mf], factor)

  X <- model.matrix(formula, data = mf,
                    contrasts.arg = lapply(Filter(is.factor, mf),
                                           function(x) contrasts(x, contrasts = FALSE)/sqrt(2)))

  if (ncol(X) > 1) {
    assign <- attr(X, "assign")[-1]
    X <- X[, -1, drop = FALSE]
  }
  attr(X, "assign") <- assign

  attr(X, "treat") <-  model.response(mf)

  X
}

.check_X <- function(X) {
  if (isTRUE(attr(X, "checked"))) return(X)

  treat <- attr(X, "treat")

  if (is.data.frame(X)) X <- as.matrix(X)
  else if (is.numeric(X) && is.null(dim(X))) {
    X <- matrix(X, nrow = length(X),
                dimnames = list(names(X), NULL))
  }

  if (anyNA(X)) .err("missing values are not allowed in the covariates")
  if (any(!is.finite(X))) .err("Non-finite values are not allowed in the covariates")

  if (!is.numeric(X) || length(dim(X)) != 2) {
    stop("bad X")
  }
  attr(X, "checked") <- TRUE
  attr(X, "treat") <- treat
  X
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
