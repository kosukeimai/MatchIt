matchit2genetic <- function(treat, X, data, distance, discarded,
                            ratio = 1, verbose = FALSE, ...) {
  if (!("rgenoud" %in% .packages(all = TRUE)))
    install.packages("rgenoud")
  require(rgenoud)

  if (!("Matching" %in% .packages(all = TRUE)))
    install.packages("Matching")
  require(Matching)

  if (verbose)
    cat("Genetic matching... \n")
  
  tt <- treat[!discarded]
  n <- length(tt)
  n1 <- length(tt[tt==1])
  xx <- X[!discarded,]
  dd <- distance[!discarded]
  tind <- (1:n)[tt==1]
  cind <- (1:n)[tt==0]
  labels <- names(tt)
  tlabels <- names(tt[tt==1])
  clabels <- names(tt[tt==0])
  out <- GenMatch(tt, cbind(dd, xx), M = ratio, ...)$matches
  ## ratio matching does not seem to work with GenMatch
  mm <- matrix(0, nrow = n1, ncol = max(table(out[,1])), dimnames =
               list(tlabels, 1:max(table(out[,1]))))
  for (i in 1:n1) {
    tmp <- labels[c(out[out[,1]==tind[i],2:(ratio+1)])]
    if (length(tmp) < ncol(mm))
      tmp <- c(tmp, rep(NA, ncol(mm)-length(tmp)))
    mm[i,] <- tmp
  }

  if (any(discarded)) {
    tdisc <- discarded[treat==1]
    tmp <- matrix(NA, nrow = sum(tdisc), ncol = ncol(mm), dimnames =
                  list(names(treat[treat == 1 & discarded]),
                       1:ncol(mm)))
    mm <- as.matrix(rbind(mm, tmp)[names(treat[treat==1]),])
  }
  
  res <- list(match.matrix = mm, weights = weights.matrix(mm, treat,
    discarded))
  class(res) <- "matchit"
  return(res)
}
