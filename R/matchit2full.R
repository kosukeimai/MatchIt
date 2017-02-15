matchit2full <- function(treat, X, data, distance, discarded, is.full.mahalanobis,
                         verbose=FALSE, ...) { 
  #if (!("optmatch" %in% .packages(all = TRUE)))
  #  install.packages("optmatch")
  requireNamespace(optmatch)
  
  if(verbose)
    cat("Full matching... \n")
  
  ## full matching for undiscarded units
  ttt <- treat[!discarded]
  ddd <- distance[!discarded]
  n0 <- length(ttt[ttt==0])
  n1 <- length(ttt[ttt==1])
  d1 <- ddd[ttt==1]
  d0 <- ddd[ttt==0]
  d <- matrix(0, ncol=n0, nrow=n1)
  rownames(d) <- names(ttt[ttt==1])
  colnames(d) <- names(ttt[ttt==0])
  for (i in 1:n1) 
    d[i,] <- abs(d1[i]-d0)
  full <- fullmatch(d, ...)
  psclass <- full[pmatch(names(ttt), names(full))]
  psclass <- as.numeric(as.factor(psclass))
  names(psclass) <- names(ttt)

  ## add psclass = NA for discarded units
  if (any(discarded)) {
    tmp <- rep(NA, sum(discarded))
    names(tmp) <- names(treat[discarded])
    psclass <- c(psclass, tmp)[names(treat)]
  }
  
  ## calculate weights and return the results
  res <- list(subclass = psclass, weights = weights.subclass(psclass, treat))
  class(res) <- c("matchit.full", "matchit")
  return(res)
}
