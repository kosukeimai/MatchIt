summary.matchit.exact <- function(obj, verbose=F, ...) {
  XX <- obj$X
  treat <- obj$treat
  qbins <- max(obj$subclass,na.rm=TRUE)
  if(!verbose){
    q.table <- as.data.frame(matrix(0,qbins,3))
    names(q.table) <- c("Treated","Control","Total")
    for(i in 1:qbins){
      qi <- obj$subclass==i
      q.table[i,] <- c(sum(treat[qi]==1), sum(treat[qi]==0), length(treat[qi]))
    }
  } else {
    kk <- ncol(XX)
    q.table <- as.data.frame(matrix(0,qbins,kk+3))
    names(q.table) <- c("Treated","Control","Total",dimnames(XX)[[2]])
    for(i in 1:qbins){
      qi <- obj$subclass==i
      q.table[i,] <- c(sum(treat[qi]==1), sum(treat[qi]==0), length(treat[qi]),as.numeric(XX[qi,,drop=F][1,]))
    }
  }
  obj$q.table <- q.table
  obj$verbose <- verbose
  class(obj) <- "summary.matchit.exact"
  obj
}
