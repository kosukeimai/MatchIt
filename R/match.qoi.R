## Functions to calculate summary stats
qoi <- function(xx, tt, ww = NULL, s.weights, subclass = NULL, mm = NULL, s.d.denom = "treated", standardize = FALSE,
                compute.pair.dist = TRUE) {

  un <- is.null(ww)
  bin.var <- all(xx == 0 | xx == 1)

  xsum <- rep(NA_real_, 7)
  if (standardize)
    names(xsum) <- c("Means Treated","Means Control", "Std. Mean Diff.",
                     "Var. Ratio", "eCDF Mean", "eCDF Max", "Std. Pair Dist.")
  else
    names(xsum) <- c("Means Treated","Means Control", "Mean Diff.",
                     "Var. Ratio", "eQQ Mean", "eQQ Max", "Pair Dist.")

  if (un) ww <- s.weights
  else ww <- ww * s.weights

  too.small <- sum(ww[tt==1] != 0) < 2 || sum(ww[tt==0] != 0) < 2

  xsum["Means Treated"] <- wm(xx[tt==1], ww[tt==1], na.rm=TRUE)
  xsum["Means Control"] <- wm(xx[tt==0], ww[tt==0], na.rm=TRUE)

  mdiff <- xsum["Means Treated"] - xsum["Means Control"]

  if (standardize && abs(mdiff) > sqrt(.Machine$double.eps)) {
    if (!too.small) {
      if (is.numeric(s.d.denom)) {
        std <- s.d.denom
      }
      else {
        s.d.denom <- match_arg(s.d.denom, c("treated", "control", "pooled"))
        std <- switch(s.d.denom,
                      "treated" = sqrt(wvar(xx[tt==1], bin.var, s.weights[tt==1])),
                      "control" = sqrt(wvar(xx[tt==0], bin.var, s.weights[tt==0])),
                      "pooled" = pooled_sd(xx, tt, w = s.weights, bin.var = bin.var, contribution = "equal"))

        #Avoid divide by zero
        if (std < sqrt(.Machine$double.eps)) {
          std <- pooled_sd(xx, tt, w = s.weights, bin.var = bin.var, contribution = "equal")
        }
      }

      xsum[3] <- mdiff/std
      if (!un && compute.pair.dist) xsum[7] <- pair.dist(xx, tt, subclass, mm, std)
    }
  }
  else {
    xsum[3] <- mdiff
    if (!un && compute.pair.dist) xsum[7] <- pair.dist(xx, tt, subclass, mm)
  }

  if (bin.var) {
    xsum[5:6] <- abs(mdiff)
  }
  else if (!too.small) {
    xsum["Var. Ratio"] <- wvar(xx[tt==1], bin.var, ww[tt==1]) / wvar(xx[tt==0], bin.var, ww[tt==0])

    qqmat <- qqsum(xx, tt, ww, standardize = standardize)
    xsum[5:6] <- qqmat[c("meandiff", "maxdiff")]
  }

  xsum
}

qoi.subclass <- function(xx, tt, s.weights, subclass, s.d.denom = "treated", standardize = FALSE, which.subclass = NULL) {
  #Within-subclass balance statistics
  bin.var <- all(xx == 0 | xx == 1)
  in.sub <- !is.na(subclass) & subclass == which.subclass

  xsum <- matrix(NA_real_, nrow = 1, ncol = 6)
  rownames(xsum) <- "Subclass"
  if (standardize)
    colnames(xsum) <- c("Means Treated","Means Control", "Std. Mean Diff.",
                        "Var. Ratio", "eCDF Mean", "eCDF Max")
  else
    colnames(xsum) <- c("Means Treated","Means Control", "Mean Diff",
                        "Var. Ratio", "eQQ Mean", "eQQ Max")

  too.small <- sum(in.sub & tt==1) < 2 || sum(in.sub & tt==0) < 2

  xsum["Subclass","Means Treated"] <- wm(xx[in.sub & tt==1], s.weights[in.sub & tt==1], na.rm=TRUE)
  xsum["Subclass","Means Control"] <- wm(xx[in.sub & tt==0], s.weights[in.sub & tt==0], na.rm=TRUE)

  mdiff <- xsum["Subclass","Means Treated"] - xsum["Subclass","Means Control"]

  if (standardize && abs(mdiff) > 1e-8) {
    if (!too.small) {
      if (is.numeric(s.d.denom)) {
        std <- s.d.denom
      }
      else {
        #SD from full sample, not within subclass
        s.d.denom <- match_arg(s.d.denom, c("treated", "control", "pooled"))
        std <- switch(s.d.denom,
                      "treated" = sqrt(wvar(xx[tt==1], bin.var, s.weights[tt==1])),
                      "control" = sqrt(wvar(xx[tt==0], bin.var, s.weights[tt==0])),
                      "pooled" = pooled_sd(xx, tt, w = s.weights, bin.var = bin.var, contribution = "equal"))
      }

      xsum["Subclass", 3] <- mdiff/std
    }
  }
  else {
    xsum["Subclass", 3] <- mdiff
  }

  if (bin.var) {
    xsum["Subclass", 5:6] <- abs(mdiff)
  }
  else if (!too.small) {
    xsum["Subclass", "Var. Ratio"] <- wvar(xx[in.sub & tt==1], bin.var, s.weights[in.sub & tt==1]) / wvar(xx[in.sub & tt==0], bin.var, s.weights[in.sub & tt==0])

    qqall <- qqsum(xx[in.sub], tt[in.sub], standardize = standardize)
    xsum["Subclass", 5:6] <- qqall[c("meandiff", "maxdiff")]
  }

  xsum
}

#Compute within-pair/subclass distances
pair.dist <- function(xx, tt, subclass = NULL, mm = NULL, std = NULL, fast = TRUE) {

  if (!is.null(mm)) {
    names(xx) <- names(tt)
    xx_t <- xx[rownames(mm)]
    xx_c <- matrix(0, nrow = nrow(mm), ncol = ncol(mm))
    xx_c[] <- xx[mm]

    mpdiff <- mean(abs(xx_t - xx_c), na.rm = TRUE)
  }
  else if (!is.null(subclass)) {
    if (!fast) {
      dists <- unlist(lapply(levels(subclass), function(s) {
        t1 <- which(!is.na(subclass) & subclass == s & tt == 1)
        t0 <- which(!is.na(subclass) & subclass == s & tt == 0)
        if (length(t1) == 1 || length(t0) == 1) {
          xx[t1] - xx[t0]
        }
        else {
          outer(xx[t1], xx[t0], "-")
        }
      }))
      mpdiff <- mean(abs(dists))
    }
    else {
      mpdiff <- pairdistsubC(as.numeric(xx), as.integer(tt),
                             as.integer(subclass), nlevels(subclass))
    }
  }
  else return(NA_real_)

  if (!is.null(std) && abs(mpdiff) > 1e-8) {
      mpdiff <- mpdiff/std
  }

  mpdiff
}

## Function for QQ summary stats
qqsum <- function(x, t, w = NULL, standardize = FALSE) {
  #x = variable, t = treat, w = weights

  n.obs <- length(x)

  if (is.null(w)) w <- rep(1, n.obs)

  if (all(x == 0 | x == 1)) {
    #For binary variables, just difference in means
    ediff <- abs(wm(x[t == t[1]], w[t == t[1]]) - wm(x[t != t[1]], w[t != t[1]]))
    return(c(meandiff = ediff, meddiff = ediff, maxdiff = ediff))
  }
  else {
    for (i in unique(t, nmax = 2)) w[t==i] <- w[t==i]/sum(w[t==i])

    ord <- order(x)
    x_ord <- x[ord]
    w_ord <- w[ord]
    t_ord <- t[ord]

    if (standardize) {
      #Difference between ecdf of x for each group
      w_ord_ <- w_ord
      w_ord_[t_ord==t_ord[1]] <- -w_ord_[t_ord==t_ord[1]]
      ediff <- abs(cumsum(w_ord_))[c(diff(x_ord) != 0, TRUE)]
    }
    else {
      #Horizontal distance of ecdf between groups
      #Need to interpolate larger group to be same size as smaller group

      u <- unique(x_ord)

      wn1 <- sum(w[t == t_ord[1]] > 0)
      wn0 <- sum(w[t != t_ord[1]] > 0)

      w1 <- w_ord[t_ord == t_ord[1]]
      w0 <- w_ord[t_ord != t_ord[1]]

      x1 <- x_ord[t_ord == t_ord[1]][w1 > 0]
      x0 <- x_ord[t_ord != t_ord[1]][w0 > 0]

      if (wn1 < wn0) {
        if (length(u) <= 5) {
          x0probs <- vapply(u, function(u_) wm(x0 == u_, w0[w0 > 0]), numeric(1L))
          x0cumprobs <- c(0, cumsum(x0probs)[-length(u)], 1)
          x0 <- u[findInterval(cumsum(w1[w1 > 0]), x0cumprobs, rightmost.closed = TRUE)]
        }
        else {
          x0 <- approx(cumsum(w0[w0 > 0]), y = x0, xout = cumsum(w1[w1 > 0]), rule = 2,
                       method = "constant", ties = "ordered")$y
        }
      }
      else {
        if (length(u) <= 5) {
          x1probs <- vapply(u, function(u_) wm(x1 == u_, w1[w1 > 0]), numeric(1L))
          x1cumprobs <- c(0, cumsum(x1probs)[-length(u)], 1)
          x1 <- u[findInterval(cumsum(w0[w0 > 0]), x1cumprobs, rightmost.closed = TRUE)]
        }
        else {
          x1 <- approx(cumsum(w1[w1 > 0]), y = x1, xout = cumsum(w0[w0 > 0]), rule = 2,
                       method = "constant", ties = "ordered")$y
        }
      }
      ediff <- abs(x1 - x0)
    }
    return(c(meandiff = mean(ediff), maxdiff = max(ediff)))
  }
}