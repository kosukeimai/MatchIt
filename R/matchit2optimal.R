matchit2optimal <- function(treat, X, data, dist, ratio = 1, ...) {
  require(optmatch)
  n0 <- length(treat[treat==0])
  n1 <- length(treat[treat==1])
  d1 <- dist[treat==1]
  d0 <- dist[treat==0]
  d <- matrix(0, ncol=n0, nrow=n1)
  tlabels <- rownames(d) <- names(treat[treat==1])
  clabels <- colnames(d) <- names(treat[treat==0])
  for (i in 1:n1) 
    d[i,] <- abs(d1[i]-d0)
  full <- fullmatch(d, min.controls = ratio,
                    max.controls = ratio, omit.fraction =
                    (n0-ratio*n1)/n0, ...)
  psclass <- full[pmatch(names(treat), names(full))]
  psclass <- as.numeric(as.factor(psclass))
  names(psclass) <- names(treat)

  mm <- matrix(0, nrow = n1, ncol = ratio, dimnames = list(tlabels, 1:ratio))
  for (i in 1:n1)
    mm[i,] <- names(which(psclass[tlabels[i]] == psclass[-pmatch(tlabels[i],
                                   names(psclass))]))

  psweights <- weights.matrix(match.matrix=as.data.frame(mm), treat=treat)$psweights
  
  res <- list(match.matrix = mm, psclass=psclass, psweights=psweights)

  class(res) <- "matchit"
  return(res)
}
