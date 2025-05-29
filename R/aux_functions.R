#Function to ensure no subclass is devoid of both treated and control units by "scooting" units
#from other subclasses.
subclass_scoot <- function(sub, treat, x, min.n = 1L) {
  #Reassigns subclasses so there are no empty subclasses
  #for each treatment group.
  subtab <- table(treat, sub)

  if (all(subtab >= min.n)) {
    return(sub)
  }

  nsub <- ncol(subtab)

  if (any(rowSums(subtab) < nsub * min.n)) {
    .err(sprintf("not enough units to fit %s treated and control %s in each subclass",
                 min.n, ngettext(min.n, "unit", "units")))
  }

  subclass_scootC(as.integer(sub), as.integer(treat),
                  as.numeric(x), as.integer(min.n))
}

#Create info component of matchit object
create_info <- function(method, fn1, link, discard, replace, ratio,
                        mahalanobis, transform, subclass, antiexact,
                        distance_is_matrix) {
  list(method = method,
       distance = if (is_not_null(fn1)) sub("distance2", "", fn1, fixed = TRUE) else NULL,
       link = if (is_not_null(link)) link else NULL,
       discard = discard,
       replace = if (is_not_null(method) && method %in% c("nearest", "genetic")) replace else NULL,
       ratio = if (is_not_null(method) && method %in% c("nearest", "optimal", "genetic")) ratio else NULL,
       max.controls = if (is_not_null(method) && method %in% c("nearest", "optimal")) attr(ratio, "max.controls") else NULL,
       mahalanobis = mahalanobis,
       transform = transform,
       subclass = if (is_not_null(method) && method == "subclass") length(unique(subclass[!is.na(subclass)])) else NULL,
       antiexact = antiexact,
       distance_is_matrix = distance_is_matrix)
}

#Function to turn a method name into a phrase describing the method
info_to_method <- function(info) {

  out.list <- setNames(vector("list", 3L), c("kto1", "type", "replace"))

  out.list[["kto1"]] <- {
    if (is_null(info$ratio)) NULL
    else sprintf("%s%s:1",
                 if (is_not_null(info$max.controls)) "variable ratio " else "",
                 round(info$ratio, 2L))
  }

  out.list[["type"]] <- {
    if (is_null(info$method)) "none (no matching)"
    else switch(info$method,
                "exact" = "exact matching",
                "cem" = "coarsened exact matching",
                "nearest" = "nearest neighbor matching",
                "optimal" = "optimal pair matching",
                "full" = "optimal full matching",
                "quick" = "generalized full matching",
                "genetic" = "genetic matching",
                "subclass" = sprintf("subclassification (%s subclasses)", info$subclass),
                "cardinality" = "cardinality matching",
                if (is_null(attr(info$method, "method"))) "an unspecified matching method"
                else attr(info$method, "method"))
  }

  out.list[["replace"]] <- {
    if (is_null(info$replace) || !info$method %in% c("nearest", "genetic")) NULL
    else if (info$replace) "with replacement"
    else "without replacement"
  }

  firstup(do.call("paste", unname(out.list)))
}

info_to_distance <- function(info) {
  distance <- info$distance
  link <- info$link
  if (is_not_null(link) && startsWith(as.character(link), "linear")) {
    linear <- TRUE
    link <- sub("linear.", "", as.character(link))
  }
  else {
    linear <- FALSE
  }

  .dist <- switch(distance,
                 "glm" = switch(link,
                                "logit" = "logistic regression",
                                "probit" = "probit regression",
                                sprintf("GLM with a %s link", link)),
                 "gam" = sprintf("GAM with a %s link", link),
                 "gbm" = "GBM",
                 "elasticnet" = sprintf("an elastic net with a %s link", link),
                 "lasso" = switch(link,
                                  "logit" = "lasso logistic regression",
                                  sprintf("lasso regression with a %s link", link)),
                 "ridge" = switch(link,
                                  "logit" = "ridge logistic regression",
                                  sprintf("ridge regression with a %s link", link)),
                 "rpart" = "CART",
                 "nnet" = "a neural network",
                 "cbps" = "CBPS",
                 "bart" = "BART",
                 "randomforest" = "a random forest")

  if (linear) {
    .dist <- sprintf("%s and linearized", .dist)
  }

  .dist
}

#Make interaction vector out of matrix of covs; similar to interaction()
exactify <- function(X, nam = NULL, sep = "|", include_vars = FALSE, justify = "right") {
  if (is_null(nam)) {
    nam <- rownames(X)
  }

  if (is.matrix(X)) {
    X <- as.data.frame.matrix(X)
  }
  else if (!is.list(X)) {
    stop("X must be a matrix, data frame, or list.")
  }

  X <- X[lengths(X) > 0L]

  if (is_null(X)) {
    return(NULL)
  }

  for (i in seq_along(X)) {
    unique_x <- {
      if (is.factor(X[[i]])) levels(X[[i]])
      else sort(unique(X[[i]]))
    }

    lev <- {
      if (include_vars) sprintf("%s = %s",
                                names(X)[i],
                                add_quotes(unique_x, chk::vld_character_or_factor(X[[i]])))
      else if (is_null(justify)) unique_x
      else format(unique_x, justify = justify)
    }

    X[[i]] <- factor(X[[i]], levels = unique_x, labels = lev)
  }

  out <- interaction2(X, sep = sep, lex.order = if (include_vars) TRUE else NULL)

  if (is_null(nam)) {
    return(out)
  }

  setNames(out, nam)
}

#Get covariates (RHS) vars from formula
get_covs_matrix <- function(formula = NULL, data = NULL) {

  if (is_null(formula)) {
    fnames <- colnames(data)
    fnames[!startsWith(fnames, "`")] <- add_quotes(fnames[!startsWith(fnames, "`")], "`")
    formula <- reformulate(fnames)
  }
  else {
    formula <- update(terms(formula, data = data), NULL ~ . + 1)
  }

  mf <- model.frame(terms(formula, data = data), data,
                    na.action = na.pass)

  chars.in.mf <- vapply(mf, is.character, logical(1L))
  for (i in which(chars.in.mf)) {
    mf[[i]] <- as.factor(mf[[i]])
  }

  mf <- droplevels(mf)

  X <- model.matrix(formula, data = mf,
                    contrasts.arg = lapply(Filter(is.factor, mf),
                                           contrasts, contrasts = FALSE))

  .assign <- attr(X, "assign")[-1L]
  X <- X[, -1L, drop = FALSE]

  attr(X, "assign") <- .assign

  X
}

#Extracts and names the "assign" attribute from get_covs_matrix()
get_assign <- function(mat) {
  if (is_null(attr(mat, "assign"))) {
    return(NULL)
  }

  setNames(attr(mat, "assign"), colnames(mat))
}

#Convert match.matrix (mm) using numerical indices to using char rownames
nummm2charmm <- function(nummm, treat) {
  #Assumes nummm has rownames
  charmm <- array(NA_character_, dim = dim(nummm),
                  dimnames = dimnames(nummm))
  charmm[] <- names(treat)[nummm]
  charmm
}

charmm2nummm <- function(charmm, treat) {
  nummm <- array(NA_integer_, dim = dim(charmm),
                 dimnames = dimnames(charmm))
  n_index <- setNames(seq_along(treat), names(treat))
  nummm[] <- n_index[charmm]
  nummm
}

#Get subclass from match.matrix. Only to be used if replace = FALSE. See subclass2mmC.cpp for reverse.
mm2subclass <- function(mm, treat, focal = NULL) {
  if (!is.integer(mm)) {
    mm <- charmm2nummm(mm, treat)
  }

  mm2subclassC(mm, treat, focal)
}

#Pooled within-group (weighted) covariance by group-mean centering covariates. Used
#in Mahalanobis distance
pooled_cov <- function(X, t, w = NULL) {
  unique_t <- unique(t)
  if (is_null(dim(X))) X <- matrix(X, nrow = length(X))

  if (is_null(w)) {
    n <- nrow(X)
    for (i in unique_t) {
      in_t <- which(t == i)
      for (j in seq_len(ncol(X))) {
        X[in_t, j] <- X[in_t, j] - mean(X[in_t, j])
      }
    }

    return(cov(X) * (n - 1) / (n - length(unique_t)))
  }

  for (i in unique_t) {
    in_t <- which(t == i)
    for (j in seq_len(ncol(X))) {
      X[in_t, j] <- X[in_t, j] - wm(X[in_t, j], w[in_t])
    }
  }

  cov.wt(X, w)$cov
}

pooled_sd <- function(X, t, w = NULL, bin.var = NULL, contribution = "proportional") {
  contribution <- match_arg(contribution, c("proportional", "equal"))
  unique_t <- unique(t)
  if (is_null(dim(X))) X <- matrix(X, nrow = length(X))
  n <- nrow(X)

  if (is_null(bin.var)) {
    bin.var <- apply(X, 2L, function(x) all(x == 0 | x == 1))
  }

  if (contribution == "equal") {
    vars <- do.call("rbind", lapply(unique_t, function(i) {
      in_t <- which(t == i)
      vapply(seq_len(ncol(X)), function(j) {
        wvar(X[in_t, j], w = w[in_t], bin.var = bin.var[j])
      }, numeric(1L))
    }))

    pooled_var <- colMeans(vars)
  }
  else {
    pooled_var <- vapply(seq_len(ncol(X)), function(j) {
      x <- X[, j]
      b <- bin.var[j]

      if (b) {
        v <- {
          if (is_null(w)) vapply(unique_t, function(i) {
            in_i <- which(t == i)
            sxi <- sum(x[in_i])
            ni <- length(in_i)
            sxi * (1 - sxi / ni) / n
          }, numeric(1L))
          else vapply(unique_t, function(i) {
            in_i <- which(t == i)
            sxi <- sum(x[in_i] * w[in_i])
            ni <- sum(w[in_i])
            sxi * (1 - sxi / ni) / sum(w)
          }, numeric(1L))
        }

        return(sum(v))
      }

      if (is_null(w)) {
        for (i in unique_t) {
          in_i <- which(t == i)
          x[in_i] <- x[in_i] - wm(x[in_i])
        }

        return(sum(x^2) / (n - length(unique_t)))
      }

      for (i in unique_t) {
        in_i <- which(t == i)
        x[in_i] <- x[in_i] - wm(x[in_i], w[in_i])
      }
      w_ <- .make_sum_to_1(w)

      sum(w_ * x^2) / (1 - sum(w_^2))

    }, numeric(1L))
  }

  setNames(sqrt(pooled_var), colnames(X))
}

#Effective sample size
ESS <- function(w) {
  sum(w)^2 / sum(w^2)
}

#Compute sample sizes
nn <- function(treat, weights, discarded = NULL, s.weights = NULL) {

  if (is_null(discarded)) {
    discarded <- rep_with(FALSE, treat)
  }

  if (is_null(s.weights)) {
    s.weights <- rep_with(1, treat)
  }

  weights <- weights * s.weights

  n <- matrix(0, ncol = 2L, nrow = 6L,
              dimnames = list(c("All (ESS)", "All", "Matched (ESS)",
                                "Matched", "Unmatched", "Discarded"),
                              c("Control", "Treated")))

  t1 <- treat == 1

  #                       Control                                        Treated
  n["All (ESS)", ] <-     c(ESS(s.weights[!t1]),                  ESS(s.weights[t1]))
  n["All", ] <-           c(sum(!t1),                             sum(t1))
  n["Matched (ESS)", ] <- c(ESS(weights[!t1]),                    ESS(weights[t1]))
  n["Matched", ] <-       c(sum(!t1 & weights > 0),               sum(t1 & weights > 0))
  n["Unmatched", ] <-     c(sum(!t1 & weights == 0 & !discarded), sum(t1 & weights == 0 & !discarded))
  n["Discarded", ] <-     c(sum(!t1 & discarded),                 sum(t1 & discarded))

  n
}

#Compute subclass sample sizes
qn <- function(treat, subclass, discarded = NULL) {

  treat <- factor(treat, levels = 0:1, labels = c("Control", "Treated"))

  if (is_null(discarded)) {
    discarded <- rep.int(FALSE, length(treat))
  }

  qn <- table(treat[!discarded], subclass[!discarded])

  if (any(is.na(subclass) & !discarded)) {
    qn <- cbind(qn, table(treat[is.na(subclass) & !discarded]))
    colnames(qn)[ncol(qn)] <- "Unmatched"
  }

  if (any(discarded)) {
    qn <- cbind(qn, table(treat[discarded]))
    colnames(qn)[ncol(qn)] <- "Discarded"
  }

  qn <- rbind(qn, colSums(qn))
  rownames(qn)[nrow(qn)] <- "Total"

  qn <- cbind(qn, rowSums(qn))
  colnames(qn)[ncol(qn)] <- "All"

  qn
}

#Function to capture and print errors and warnings better
matchit_try <- function(expr, from = NULL, dont_warn_if = NULL) {
  tryCatch({
    withCallingHandlers({
      expr
    },
    warning = function(w) {
      if (is_null(dont_warn_if) || !any(vapply(dont_warn_if, grepl, logical(1L), conditionMessage(w), fixed = TRUE))) {
        if (is_null(from)) .wrn(conditionMessage(w), tidy = FALSE)
        else .wrn(sprintf("(from %s) %s", from, conditionMessage(w)), tidy = FALSE)
      }
      invokeRestart("muffleWarning")
    })},
    error = function(e) {
      if (is_null(from)) .err(conditionMessage(e), tidy = FALSE)
      else .err(sprintf("(from %s) %s", from, conditionMessage(e)), tidy = FALSE)
    })
}
