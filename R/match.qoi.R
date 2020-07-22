## Function to calculate summary stats
qoi <- function(xx, tt, ww = NULL, s.d.denom = "treated", standardize = FALSE) {

  full.only <- is.null(ww)

  xsum <- matrix(NA_real_, nrow = 1, ncol = 6)
  rownames(xsum) <- "Full"
  if (standardize)
    colnames(xsum) <- c("Means Treated","Means Control", "Std. Mean Diff.",
                        "eCDF Med", "eCDF Mean", "eCDF Max")
  else
    colnames(xsum) <- c("Means Treated","Means Control", "Mean Diff",
                        "eQQ Med", "eQQ Mean", "eQQ Max")
  x1 <- xx[tt==1]
  x0 <- xx[tt==0]

  too.small <- length(x1) < 2 || length(x0) < 2

  xsum["Full","Means Treated"] <- mean(x1, na.rm=TRUE)
  xsum["Full","Means Control"] <- mean(x0, na.rm=TRUE)

  mdiff <- xsum["Full","Means Treated"] - xsum["Full","Means Control"]

  if (standardize && abs(mdiff) > 1e-8) {
    if (!too.small) {
      std <- switch(s.d.denom,
                    "treated" = sd(x0, na.rm = TRUE),
                    "control" = sd(x0, na.rm = TRUE),
                    "pooled" = sqrt(.5*(var(x1, na.rm = TRUE) + var(x0, na.rm = TRUE))),
                    stop("s.d.denom must be one of \"treated\", \"control\", or \"pooled\".", call. = FALSE))

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
    xsum.matched <- matrix(NA_real_, nrow = 1, ncol = 6)
    rownames(xsum.matched) <- "Matched"
    colnames(xsum.matched) <- colnames(xsum)

    too.small <- sum(ww[tt==1] > 0) < 2 || sum(ww[tt==0] > 0) < 2

    xsum.matched["Matched","Means Treated"] <- weighted.mean(x1, ww[tt==1], na.rm=TRUE)
    xsum.matched["Matched","Means Control"] <- weighted.mean(x0, ww[tt==0], na.rm=TRUE)

    if (!(length(x1) < 2 || (length(x0) < 2))) {

      mdiff <- xsum.matched["Matched","Means Treated"] - xsum.matched["Matched","Means Control"]

      if (standardize && abs(mdiff) > 1e-8) {
        if (!too.small) {
          std <- switch(s.d.denom,
                        "treated" = sd(x0, na.rm = TRUE),
                        "control" = sd(x0, na.rm = TRUE),
                        "pooled" = sqrt(.5*(var(x1, na.rm = TRUE) + var(x0, na.rm = TRUE))),
                        stop("s.d.denom must be one of \"treated\", \"control\", or \"pooled\".", call. = FALSE))

          xsum.matched["Matched",3] <- mdiff/std
        }
      }
      else {
        xsum.matched["Matched",3] <- mdiff
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