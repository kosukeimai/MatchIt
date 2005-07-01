matchit2full <- function(treat, X, dist, ...) {
  require(optmatch)
  n0 <- length(treat[treat==0])
  n1 <- length(treat[treat==1])
  d1 <- dist[treat==1]
  d0 <- dist[treat==0]
  d <- matrix(0, ncol=n0, nrow=n1)
  rownames(d) <- names(treat[treat==1])
  colnames(d) <- names(treat[treat==0])
  for (i in 1:n1) 
    d[i,] <- abs(d1[i]-d0)
  full <- fullmatch(d, ...)
  psclass <- full[pmatch(names(treat), names(full))]
  psclass <- as.numeric(as.factor(psclass))
  names(psclass) <- names(treat)

  psweights <- weights.subclass(psclass, treat)$psweights

  res <- list(psclass = psclass, q.cut = NULL, psweights=psweights)
  class(res) <- c("matchit.subclass", "matchit")
  return(res)
}
