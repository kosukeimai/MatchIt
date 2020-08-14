## Functions to calculate summary stats
qoi <- function(xx, tt, ww = NULL, subclass = NULL, mm = NULL, s.d.denom = "treated", standardize = FALSE) {

  full.only <- is.null(ww) && is.null(mm) && is.null(subclass)

  xsum <- matrix(NA_real_, nrow = 1, ncol = 7)
  rownames(xsum) <- "Full"
  if (standardize)
    colnames(xsum) <- c("Means Treated","Means Control", "Std. Mean Diff.",
                        "eCDF Med", "eCDF Mean", "eCDF Max", "Std. Pair Diff.")
  else
    colnames(xsum) <- c("Means Treated","Means Control", "Mean Diff",
                        "eQQ Med", "eQQ Mean", "eQQ Max", "Pair Diff")

  too.small <- sum(tt==1) < 2 || sum(tt==0) < 2

  xsum["Full","Means Treated"] <- mean(xx[tt==1], na.rm=TRUE)
  xsum["Full","Means Control"] <- mean(xx[tt==0], na.rm=TRUE)

  mdiff <- xsum["Full","Means Treated"] - xsum["Full","Means Control"]

  if (standardize && abs(mdiff) > 1e-8) {
    if (!too.small) {
      if (is.numeric(s.d.denom)) {
        std <- s.d.denom
      }
      else {
        s.d.denom <- match_arg(s.d.denom, c("treated", "control", "pooled"))
        std <- switch(s.d.denom,
                      "treated" = sd_(xx[tt==1]),
                      "control" = sd_(xx[tt==0]),
                      "pooled" = sqrt(.5*(sd_(xx[tt==1])^2 + sd_(xx[tt==0])^2)))
      }
      xsum["Full", 3] <- mdiff/std
    }
  }
  else {
    xsum["Full", 3] <- mdiff
  }

  if (!too.small) {
    qqall <- qqsum(xx, tt, standardize = standardize)
    xsum["Full", 4:6] <- qqall[c("meddiff", "meandiff", "maxdiff")]
  }

  if (!full.only) {
    xsum.matched <- matrix(NA_real_, nrow = 1, ncol = 7)
    rownames(xsum.matched) <- "Matched"
    colnames(xsum.matched) <- colnames(xsum)

    too.small <- sum(ww[tt==1] > 0) < 2 || sum(ww[tt==0] > 0) < 2

    xsum.matched["Matched","Means Treated"] <- weighted.mean(xx[tt==1], ww[tt==1], na.rm=TRUE)
    xsum.matched["Matched","Means Control"] <- weighted.mean(xx[tt==0], ww[tt==0], na.rm=TRUE)

    if (!(sum(tt==1) < 2 || sum(tt==0) < 2)) {

      mdiff <- xsum.matched["Matched","Means Treated"] - xsum.matched["Matched","Means Control"]

      if (standardize && abs(mdiff) > 1e-8) {
        if (!too.small) {
          xsum.matched["Matched",3] <- mdiff/std
          xsum.matched["Matched", 7] <- pair.dist(xx, tt, subclass, mm, std)
        }
      }
      else {
        xsum.matched["Matched", 3] <- mdiff
        xsum.matched["Matched", 7] <- pair.dist(xx, tt, subclass, mm)

      }

      if (!too.small) {
        qqmat <- qqsum(xx, tt, ww, standardize = standardize)
        xsum.matched["Matched", 4:6] <- qqmat[c("meddiff", "meandiff", "maxdiff")]
      }
    }
    xsum <- rbind(xsum, xsum.matched)
  }
  xsum
}

qoi.subclass <- function(xx, tt, subclass, s.d.denom = "treated", standardize = FALSE, which.subclass = NULL) {
  #Within-subclass balance statistics

  in.sub <- !is.na(subclass) & subclass == which.subclass

  xsum <- matrix(NA_real_, nrow = 1, ncol = 6)
  rownames(xsum) <- "Subclass"
  if (standardize)
    colnames(xsum) <- c("Means Treated","Means Control", "Std. Mean Diff.",
                        "eCDF Med", "eCDF Mean", "eCDF Max")
  else
    colnames(xsum) <- c("Means Treated","Means Control", "Mean Diff",
                        "eQQ Med", "eQQ Mean", "eQQ Max")

  too.small <- sum(in.sub & tt==1) < 2 || sum(in.sub & tt==0) < 2

  xsum["Subclass","Means Treated"] <- mean(xx[in.sub & tt==1], na.rm=TRUE)
  xsum["Subclass","Means Control"] <- mean(xx[in.sub & tt==0], na.rm=TRUE)

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
                      "treated" = sd_(xx[tt==1]),
                      "control" = sd_(xx[tt==0]),
                      "pooled" = sqrt(.5*(sd_(xx[tt==1])^2 + sd_(xx[tt==0])^2)))
      }

      xsum["Subclass", 3] <- mdiff/std
    }
  }
  else {
    xsum["Subclass", 3] <- mdiff
  }

  if (!too.small) {
    qqall <- qqsum(xx[in.sub], tt[in.sub], standardize = standardize)
    xsum["Subclass", 4:6] <- qqall[c("meddiff", "meandiff", "maxdiff")]
  }

  xsum
}

#Compute within-pair/subclass distances
pair.dist <- function(xx, tt, subclass = NULL, mm = NULL, std = NULL) {

  if (!is.null(subclass)) {
    unique.sub <- unique(subclass[!is.na(subclass)])
    dists <- unlist(lapply(unique.sub, function(s) {
      t1 <- which(!is.na(subclass) & subclass == s & tt == 1)
      t0 <- which(!is.na(subclass) & subclass == s & tt == 0)
      if (length(t1) == 1 || length(t0) == 1) {
        xx[t1] - xx[t0]
      }
      else {
        outer(xx[t1], xx[t0], "-")
      }
    }))
  }
  else if (!is.null(mm)) {
    names(xx) <- names(tt)
    dists <- unlist(lapply(rownames(mm)[!is.na(mm[,1])], function(t1) {
      t0 <- mm[t1, !is.na(mm[t1,])]
      if (length(t1) == 1 || length(t0) == 1) {
        xx[t1] - xx[t0]
      }
      else {
        outer(xx[t1], xx[t0], "-")
      }
    }))
  }
  else return(NA_real_)

  mpdiff <- mean(abs(dists))

  if (!is.null(std) && abs(mpdiff) > 1e-8) {
      mpdiff <- mpdiff/std
  }

  mpdiff
}