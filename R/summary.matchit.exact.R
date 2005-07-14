summary.matchit.exact <- function(object, covariates = FALSE, ...) {
  XX <- object$X
  treat <- object$treat
  qbins <- max(object$subclass,na.rm=TRUE)
  if(!covariates){
    q.table <- as.data.frame(matrix(0,qbins,3))
    names(q.table) <- c("Treated","Control","Total")
    for(i in 1:qbins){
      qi <- object$subclass==i
      q.table[i,] <- c(sum(treat[qi]==1, na.rm=T), sum(treat[qi]==0, na.rm=T), length(treat[qi & !is.na(qi)]))
    }
  } else {
    kk <- ncol(XX)
    q.table <- as.data.frame(matrix(0,qbins,kk+3))
    names(q.table) <- c("Treated","Control","Total",dimnames(XX)[[2]])
    for(i in 1:qbins){
      qi <- object$subclass==i
      q.table[i,] <- c(sum(treat[qi]==1, na.rm=T), sum(treat[qi]==0, na.rm=T), length(treat[qi & !is.na(qi)]),as.numeric(XX[qi,,drop=F][1,]))
    }
  }

  ## output
  object$q.table <- q.table
  class(object) <- c("summary.matchit.exact", "summary.matchit")
  return(object)
}
