matchit2full <- function(treat, X, data, pscore, discarded, ...) {
  require(optmatch)

  ## full matching for undiscarded units
  ttt <- treat[!discarded]
  n0 <- length(ttt[ttt==0])
  n1 <- length(ttt[ttt==1])
  d1 <- pscore[ttt==1]
  d0 <- pscore[ttt==0]
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
  tmp <- rep(NA, sum(discarded))
  names(tmp) <- names(treat[discarded])
  psclass <- c(psclass, tmp)[names(treat)]

  ## calculate weights and return the results
  psweights <- weights.subclass(psclass, treat)$psweights
  res <- list(psclass = psclass, q.cut = NULL, psweights=psweights)
  class(res) <- c("matchit.subclass", "matchit")
  return(res)
}
