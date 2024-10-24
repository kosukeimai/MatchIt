#Function to ensure no subclass is devoid of both treated and control units by "scooting" units
#from other subclasses.
subclass_scoot <- function(sub, treat, x, min.n = 1) {
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
  info <- list(method = method,
               distance = if (is_null(fn1)) NULL else sub("distance2", "", fn1, fixed = TRUE),
               link = if (is_null(link)) NULL else link,
               discard = discard,
               replace = if (is_not_null(method) && method %in% c("nearest", "genetic")) replace else NULL,
               ratio = if (is_not_null(method) && method %in% c("nearest", "optimal", "genetic")) ratio else NULL,
               max.controls = if (is_not_null(method) && method %in% c("nearest", "optimal")) attr(ratio, "max.controls") else NULL,
               mahalanobis = mahalanobis,
               transform = transform,
               subclass = if (is_not_null(method) && method == "subclass") length(unique(subclass[!is.na(subclass)])) else NULL,
               antiexact = antiexact,
               distance_is_matrix = distance_is_matrix)
  info
}

#Function to turn a method name into a phrase describing the method
info.to.method <- function(info) {

  out.list <- setNames(vector("list", 3), c("kto1", "type", "replace"))

  out.list[["kto1"]] <- {
    if (is_not_null(info$ratio)) paste0(if (is_not_null(info$max.controls)) "variable ratio ", round(info$ratio, 2), ":1")
    else NULL
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

  firstup(do.call("paste", c(unname(out.list), list(sep = " "))))
}

info.to.distance <- function(info) {
  distance <- info$distance
  link <- info$link
  if (is_not_null(link) && startsWith(as.character(link), "linear")) {
    linear <- TRUE
    link <- sub("linear.", "", as.character(link))
  }
  else {
    linear <- FALSE
  }

  if (distance == "glm") {
    if (link == "logit") dist <- "logistic regression"
    else if (link == "probit") dist <- "probit regression"
    else dist <- sprintf("GLM with a %s link", link)
  }
  else if (distance == "gam") {
    dist <- sprintf("GAM with a %s link", link)
  }
  else if (distance == "gbm") {
    dist <- "GBM"
  }
  else if (distance == "elasticnet") {
    dist <- sprintf("an elastic net with a %s link", link)
  }
  else if (distance == "lasso") {
    if (link == "logit") dist <- "lasso logistic regression"
    else dist <- sprintf("lasso regression with a %s link", link)
  }
  else if (distance == "ridge") {
    dist <- sprintf("ridge regression with a %s link", link)
  }
  else if (distance == "rpart") {
    dist <- "CART"
  }
  else if (distance == "nnet") {
    dist <- "a neural network"
  }
  else if (distance == "cbps") {
    dist <- "CBPS"
  }
  else if (distance == "bart") {
    dist <- "BART"
  }
  else if (distance == "randomforest") {
    dist <- "a random forest"
  }

  if (linear) dist <- paste(dist, "and linearized")

  dist
}

#Make interaction vector out of matrix of covs; similar to interaction()
exactify <- function(X, nam = NULL, sep = "|", include_vars = FALSE, justify = "right") {
  if (is_null(nam)) {
    nam <- rownames(X)
  }

  if (is.matrix(X)) {
    X <- setNames(lapply(seq_len(ncol(X)), function(i) X[,i]), colnames(X))
  }

  if (!is.list(X)) {
    stop("X must be a matrix, data frame, or list.")
  }

  for (i in seq_along(X)) {
    unique_x <- {
      if (is.factor(X[[i]])) levels(X[[i]])
      else sort(unique(X[[i]]))
    }

    lev <- {
      if (include_vars) {
        sprintf("%s = %s",
                names(X)[i],
                add_quotes(unique_x, is.character(X[[i]]) || is.factor(X[[i]])))
      }
      else if (is_null(justify)) unique_x
      else format(unique_x, justify = justify)
    }

    X[[i]] <- factor(X[[i]], levels = unique_x, labels = lev)
  }

  all_levels <- do.call("paste", c(rev(expand.grid(rev(lapply(X, levels)),
                                                   KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)),
                                   sep = sep))

  out <- do.call("paste", c(X, sep = sep))

  out <- factor(out, levels = all_levels[all_levels %in% out])

  if (is_not_null(nam)) names(out) <- nam

  out
}

#Get covariates (RHS) vars from formula
get.covs.matrix <- function(formula = NULL, data = NULL) {

  if (is_null(formula)) {
    fnames <- colnames(data)
    fnames[!startsWith(fnames, "`")] <- paste0("`", fnames[!startsWith(fnames, "`")], "`")
    formula <- reformulate(fnames)
  }
  else {
    formula <- update(terms(formula, data = data), NULL ~ . + 1)
  }

  mf <- model.frame(terms(formula, data = data), data,
                    na.action = na.pass)

  chars.in.mf <- vapply(mf, is.character, logical(1L))
  mf[chars.in.mf] <- lapply(mf[chars.in.mf], factor)

  mf <- droplevels(mf)

  X <- model.matrix(formula, data = mf,
                    contrasts.arg = lapply(Filter(is.factor, mf),
                                           contrasts, contrasts = FALSE))
  assign <- attr(X, "assign")[-1]
  X <- X[,-1, drop = FALSE]
  attr(X, "assign") <- assign

  X
}

#Extracts and names the "assign" attribute from get.covs.matrix()
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
    return(cov(X)*(n-1)/(n-length(unique_t)))
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
  if (is_null(bin.var)) bin.var <- apply(X, 2, function(x) all(x == 0 | x == 1))

  if (contribution == "equal") {
    vars <- matrix(0, nrow = length(unique_t), ncol = ncol(X))
    for (i in seq_along(unique_t)) {
      in_t <- which(t == unique_t[i])
      vars[i,] <- vapply(seq_len(ncol(X)), function(j) {
        x <- X[,j]
        b <- bin.var[j]
        wvar(x[in_t], w = w[in_t], bin.var = b)
      }, numeric(1L))
    }
    pooled_var <- colMeans(vars)
  }
  else {
    pooled_var <- vapply(seq_len(ncol(X)), function(j) {
      x <- X[,j]
      b <- bin.var[j]

      if (b) {
        if (is_null(w)) {
          v <- vapply(unique_t, function(i) {
            sxi <- sum(x[t == i])
            ni <- sum(t == i)
            sxi * (1 - sxi/ni) / n
          }, numeric(1L))
          return(sum(v))
        }
        else {
          v <- vapply(unique_t, function(i) {
            sxi <- sum(x[t == i] * w[t == i])
            ni <- sum(w[t==i])
            sxi * (1 - sxi/ni) / sum(w)
          }, numeric(1L))
          return(sum(v))
        }
      }
      else {
        if (is_null(w)) {
          for (i in unique_t) {
            x[t==i] <- x[t==i] - wm(x[t==i])
          }
          return(sum(x^2)/(n - length(unique_t)))
        }
        else {
          for (i in unique_t) {
            x[t==i] <- x[t==i] - wm(x[t==i], w[t==i])
          }
          w_ <- .make_sum_to_1(w)
          return(sum(w_ * x^2)/(1 - sum(w_^2)))
        }
      }
    }, numeric(1L))
  }

  setNames(sqrt(pooled_var), colnames(X))
}

#Effective sample size
ESS <- function(w) {
  sum(w)^2/sum(w^2)
}

#Compute sample sizes
nn <- function(treat, weights, discarded = NULL, s.weights = NULL) {

  if (is_null(discarded)) discarded <- rep(FALSE, length(treat))
  if (is_null(s.weights)) s.weights <- rep(1, length(treat))
  weights <- weights * s.weights
  n <- matrix(0, ncol=2, nrow=6, dimnames = list(c("All (ESS)", "All", "Matched (ESS)",
                                                   "Matched", "Unmatched","Discarded"),
                                                 c("Control", "Treated")))

  #                      Control                                    Treated
  n["All (ESS)",] <-     c(ESS(s.weights[treat==0]),                ESS(s.weights[treat==1]))
  n["All",] <-           c(sum(treat==0),                           sum(treat==1))
  n["Matched (ESS)",] <- c(ESS(weights[treat==0]),                  ESS(weights[treat==1]))
  n["Matched",] <-       c(sum(treat==0 & weights > 0),             sum(treat==1 & weights > 0))
  n["Unmatched",] <-     c(sum(treat==0 & weights==0 & !discarded), sum(treat==1 & weights==0 & !discarded))
  n["Discarded",] <-     c(sum(treat==0 & discarded),               sum(treat==1 & discarded))

  n
}

#Compute subclass sample sizes
qn <- function(treat, subclass, discarded = NULL) {

  treat <- factor(treat, levels = 0:1, labels = c("Control", "Treated"))
  if (is_null(discarded)) discarded <- rep(FALSE, length(treat))
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
      if (is_null(dont_warn_if) || !grepl(dont_warn_if, conditionMessage(w), fixed = TRUE)) {
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