matchit2full <- function(treat, dist, ...) {
  require(optmatch)
  n0 <- length(treat[treat==0])
  n1 <- length(treat[treat==1])
  d1 <- dist[treat==1]
  d0 <- dist[treat==0]
  d <- matrix(NA, ncol=n0, nrow=n1)
  rownames(d) <- names(treat[treat==1])
  colnames(d) <- names(treat[treat==0])
  for (i in 1:n1)
    d[i,] <- abs(d1[i]-d0)
  return(fullmatch(d, ...))
}
