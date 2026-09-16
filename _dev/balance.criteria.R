#Criterion
#Criteria to use in methods that involve a tuning parameter. Also called "stop.method" in twang.
#Separate functions for different treatment types.

#init functions
init_smd <- function(covs, treat, s.weights = NULL, estimand = "ATT", ...) {

  bin.var <- setNames(logical(ncol(covs)), colnames(covs))
  std <- setNames(numeric(ncol(covs)), colnames(covs))

  if (is.null(s.weights)) s.weights <- rep(1, NROW(covs))

  for (i in seq_len(ncol(covs))) {
    xx <- covs[,i]

    bin.var[i] <- all(xx == 0 | xx == 1)

    std[i] <- switch(estimand,
                  "ATT" = sqrt(wvar(xx[treat==1], bin.var[i], s.weights[treat==1])),
                  "ATC" = sqrt(wvar(xx[treat==0], bin.var[i], s.weights[treat==0])),
                  "ATE" = sqrt(.5*(wvar(xx[treat==1], bin.var[i], s.weights[treat==1]) +
                                     wvar(xx[treat==0], bin.var[i], s.weights[treat==0]))))

    if (std[i] < sqrt(.Machine$double.eps)) std[i] <- sqrt(wvar(xx, bin.var[i], s.weights)) #Avoid divide by zero
  }


  out <- list(treat = treat,
              covs = covs,
              bin.var = bin.var,
              s.weights = s.weights,
              std = std)
  class(out) <- "init_smd"
  out
}
init_ks <- function(covs, treat, s.weights = NULL, ...) {
  bin.var <- vapply(seq_len(ncol(covs)), function(i) {
    xx <- covs[,i]
    all(xx == 0 | xx == 1)
  }, logical(1L))

  if (is.null(s.weights)) s.weights <- rep(1, NROW(covs))

  bin.covs <- covs[,bin.var, drop = FALSE]
  cont.covs <- covs[,!bin.var, drop = FALSE]
  cont.covs_ord <- matrix(0L, nrow = NROW(cont.covs), ncol = NCOL(cont.covs),
                          dimnames = dimnames(cont.covs))
  for (i in seq_len(NCOL(cont.covs))) {
    cont.covs_ord[,i] <- order(cont.covs[,i])
  }

  out <- list(treat = treat,
              bin.covs = bin.covs,
              cont.covs = cont.covs,
              cont.covs_ord = cont.covs_ord,
              s.weights = s.weights)
  class(out) <- "init_ks"
  out
}
init_smd.ks <- function(covs, treat, s.weights = NULL, estimand = "ATT", ...) {
  bin.var <- vapply(seq_len(ncol(covs)), function(i) {
    xx <- covs[,i]
    all(xx == 0 | xx == 1)
  }, logical(1L))

  if (is.null(s.weights)) s.weights <- rep(1, NROW(covs))

  bin.covs <- covs[,bin.var, drop = FALSE]
  cont.covs <- covs[,!bin.var, drop = FALSE]
  cont.covs_ord <- matrix(0L, nrow = NROW(cont.covs), ncol = NCOL(cont.covs),
                          dimnames = dimnames(cont.covs))
  for (i in seq_len(NCOL(cont.covs))) {
    cont.covs_ord[,i] <- order(cont.covs[,i])
  }

  bin.std <- setNames(numeric(ncol(bin.covs)), colnames(bin.covs))
  for (i in seq_len(ncol(bin.covs))) {
    xx <- bin.covs[,i]
    bin.std[i] <- switch(estimand,
                     "ATT" = sqrt(wvar(xx[treat==1], TRUE, s.weights[treat==1])),
                     "ATC" = sqrt(wvar(xx[treat==0], TRUE, s.weights[treat==0])),
                     "ATE" = sqrt(.5*(wvar(xx[treat==1], TRUE, s.weights[treat==1]) +
                                        wvar(xx[treat==0], TRUE, s.weights[treat==0]))))

    if (bin.std[i] < sqrt(.Machine$double.eps)) bin.std[i] <- sqrt(wvar(xx, TRUE, s.weights)) #Avoid divide by zero
  }

  cont.std <- setNames(numeric(ncol(cont.covs)), colnames(cont.covs))
  for (i in seq_len(ncol(cont.covs))) {
    xx <- cont.covs[,i]
    cont.std[i] <- switch(estimand,
                         "ATT" = sqrt(wvar(xx[treat==1], FALSE, s.weights[treat==1])),
                         "ATC" = sqrt(wvar(xx[treat==0], FALSE, s.weights[treat==0])),
                         "ATE" = sqrt(.5*(wvar(xx[treat==1], FALSE, s.weights[treat==1]) +
                                            wvar(xx[treat==0], FALSE, s.weights[treat==0]))))

    if (cont.std[i] < sqrt(.Machine$double.eps)) cont.std[i] <- sqrt(wvar(xx, FALSE, s.weights)) #Avoid divide by zero
  }

  out <- list(treat = treat,
              bin.covs = bin.covs,
              cont.covs = cont.covs,
              cont.covs_ord = cont.covs_ord,
              s.weights = s.weights,
              bin.std = bin.std,
              cont.std = cont.std)
  class(out) <- "init_smd.ks"
  out
}
init_energy.dist <- function(covs, treat, s.weights = NULL, estimand = "ATT", ...) {
  if (is.null(s.weights)) s.weights <- rep(1, NROW(covs))

  covs <- transform_covariates(data = covs, treat = treat,
                               s.weights = s.weights,
                               method = "scaled_euclidean")

  d <- eucdist_internal(covs)

  n <- length(treat)

  J0 <- (treat == 0)
  J1 <- (treat == 1)

  M <- (2 * tcrossprod(J1, J0) - tcrossprod(J1, J1) - tcrossprod(J0, J0)) * d

  #Edist =  w %*% (M10 + M00 + M11) %*% t(w)

  out <- list(M = M,
              s.weights = s.weights,
              treat = treat)
  class(out) <- "init_energy.dist"
  out
}

#Statistics
smd.binary <- function(init, weights = NULL) {
  weights <- weights * init$s.weights

  vapply(seq_len(ncol(init$covs)), function(i) {
    xx <- init$covs[,i]
    m0 <- wm(xx[init$treat==0], weights[init$treat==0], na.rm=TRUE)
    m1 <- wm(xx[init$treat==1], weights[init$treat==1], na.rm=TRUE)
    abs(m1 - m0)/init$std[i]
  }, numeric(1L))
}
ks.binary <- function(init, weights = NULL) {
  weights <- weights * init$s.weights

  for (i in 0:1) weights[init$treat==i] <- weights[init$treat==i]/sum(weights[init$treat==i])

  bin.ks <- vapply(seq_len(NCOL(init$bin.covs)), function(i) {
    x <- init$bin.covs[,i]
    abs(wm(x[init$treat == 1], weights[init$treat == 1]) -
          wm(x[init$treat == 0], weights[init$treat == 0]))
  }, numeric(1L))

  cont.ks <- vapply(seq_len(NCOL(init$cont.covs)), function(i) {
    x <- init$cont.covs[,i]
    ord <- init$cont.covs_ord[,i]

    x_ord <- x[ord]
    w_ord <- weights[ord]
    t_ord <- init$treat[ord]

    w_ord_ <- w_ord
    w_ord_[t_ord==0] <- -w_ord_[t_ord==0]
    ediff <- abs(cumsum(w_ord_))[c(diff(x_ord) != 0, TRUE)]
    max(ediff)
  }, numeric(1L))

  c(bin.ks, cont.ks)
}
smd.ks.binary <- function(init, weights = NULL) {

  ks <- ks.binary(init, weights)

  weights <- weights * init$s.weights

  bin.smds <- vapply(seq_len(ncol(init$bin.covs)), function(i) {
    xx <- init$bin.covs[,i]
    m0 <- wm(xx[init$treat==0], weights[init$treat==0], na.rm=TRUE)
    m1 <- wm(xx[init$treat==1], weights[init$treat==1], na.rm=TRUE)
    abs(m1 - m0)/init$bin.std[i]
  }, numeric(1L))

  cont.smds <- vapply(seq_len(ncol(init$cont.covs)), function(i) {
    xx <- init$bin.covs[,i]
    m0 <- wm(xx[init$treat==0], weights[init$treat==0], na.rm=TRUE)
    m1 <- wm(xx[init$treat==1], weights[init$treat==1], na.rm=TRUE)
    abs(m1 - m0)/init$cont.std[i]
  }, numeric(1L))

  c(cont.smds, bin.smds, ks)
}
energy.dist.binary <- function(init, weights = NULL) {

  if (is.null(weights)) weights <- rep(1, nrow(init[["M2"]]))

  weights <- weights * init[["s.weights"]]

  for (t in 0:1) weights[init[["treat"]] == t] <- weights[init[["treat"]] == t]/sum(weights[init[["treat"]] == t])

  return(drop(t(weights) %*% init[["M"]] %*% weights))
}

initialize_balance <- function(criterion, covs, treat, s.weights = NULL, ...) {
  mc <- match.call()
  init_fun <- switch(criterion,
                     "smd.mean" =, "smd.max" =, "smd.rms" = "init_smd",
                     "ks.mean" =, "ks.max" =, "ks.rms" = "init_ks",
                     "smd.ks.mean" =, "smd.ks.max" =, "smd.ks.rms" = "init_smd.ks",
                     "pvals" = "init_pvals",
                     "energy.dist" = "init_energy.dist")

  mc[[1]] <- call(init_fun)[[1]]

  init <- eval.parent(mc)
  attr(init, "criterion") <- criterion
  attr(init, "lexical") <- FALSE#endsWith(criterion, ".max")
  init
}
compute_balance <- function(init, weights) {
  balance <- switch(attr(init, "criterion"),
      smd.mean =  mean(smd.binary(init, weights)),
      # smd.max = sort(abs(smd.binary(init, weights)), decreasing = TRUE),
      smd.max = max(smd.binary(init, weights)),
      smd.rms = sqrt(mean(smd.binary(init, weights)^2)),
      ks.mean = mean(ks.binary(init, weights)),
      # ks.max = sort(ks.binary(init, weights), decreasing = TRUE),
      ks.max = max(ks.binary(init, weights)),
      ks.rms = sqrt(mean(ks.binary(init, weights)^2)),
      smd.ks.mean = mean(smd.ks.binary(init, weights)),
      # ks.max = sort(ks.binary(init, weights), decreasing = TRUE),
      smd.ks.max = max(smd.ks.binary(init, weights)),
      smd.ks.rms = sqrt(mean(smd.ks.binary(init, weights)^2)),
      energy.dist = energy.dist.binary(init, weights)
    )

  balance
}

bal_criterion.to.phrase <- function(criterion) {

  phrase <- switch(criterion,
                   "es.mean" = "average absolute standardized mean difference",
                   "es.max" = "maximum absolute standardized mean difference",
                   "es.rms" = "root-mean-square absolute standardized mean difference",
                   "ks.mean" = "average Kolmogorov-Smirnov statistic",
                   "ks.max" = "maximum Kolmogorov-Smirnov statistic",
                   "ks.rms" = "root-mean-square Kolmogorov-Smirnov statistic",
                   # "mahalanobis" = "sample Mahalanobis distance",
                   "energy.dist" = "energy distance",
                   # "kernel.dist" = "kernel distance",
                   # "L1.med" = "L1 median",
                   # "r2" = "post-weighting treatment R-squared",
                   NA_character_
  )
  if (anyNA(phrase)) stop(paste0("\"", criterion, "\" is not an allowed criterion."))
  return(phrase)
}