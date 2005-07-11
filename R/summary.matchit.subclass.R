summary.matchit.subclass <- function(obj, verbose=F, ...) {
  XX <- obj$X
  treat <- obj$treat
  if(is.null(obj$weights)){
    weights <- rep(1,length(treat))
    names(weights) <- names(treat)
  } else {
    weights <- obj$weights
  }
  nam <- dimnames(obj$X)[[2]]
  kk <- ncol(XX)
  ## Summary Stats
  aa <- apply(XX,2,qoi,tt=treat,ww=weights)
  sum.all <- as.data.frame(matrix(0,kk,7))
  sum.matched <- as.data.frame(matrix(0,kk,7))
  row.names(sum.all) <- row.names(sum.matched) <- nam
  names(sum.all) <- names(sum.matched) <- names(aa[[1]])
  sum.all.int <- sum.matched.int <- NULL
  for(i in 1:kk){
    sum.all[i,] <- aa[[i]][1,]
    sum.matched[i,] <- aa[[i]][2,]
    if(verbose){
      for(j in i:kk){
        x2 <- XX[,i]*as.matrix(XX[,j])
        jqoi <- qoi(x2,tt=treat,ww=weights)
        sum.all.int <- rbind(sum.all.int,jqoi[1,])
        sum.matched.int <- rbind(sum.matched.int,jqoi[2,])
        row.names(sum.all.int)[nrow(sum.all.int)] <-
          row.names(sum.matched.int)[nrow(sum.matched.int)] <-
            paste(nam[i],nam[j],sep="x")
      }
    }
  }
  xn <- aa[[1]]$xn
  sum.all <- rbind(sum.all,sum.all.int)
  sum.matched <- rbind(sum.matched,sum.matched.int)
  ## By Subclass
  qbins <- max(obj$subclass,na.rm=TRUE)
  if(verbose){
    q.table <- array(0,dim=c(kk+sum(1:kk),7,qbins))
    ii <- 0
    nn <- NULL
  } else {
    q.table <- array(0,dim=c(kk,7,qbins))
  }
  aa <- apply(XX,2,qoi.by.sub,tt=treat,ww=weights,
                  qq=obj$subclass)
  for(i in 1:kk){
    if(!verbose){
      q.table[i,,] <- as.matrix(aa[[i]]$q.table)
      nn <- names(aa)
        } else {
          ii <- ii + 1 
          q.table[ii,,] <- as.matrix(aa[[i]]$q.table)
          nn <- c(nn,names(aa)[i])
          for(j in i:kk){
            ii <- ii + 1 
            x2 <- covariates[,i]*as.matrix(covariates[,j])
            q.table[ii,,] <- as.matrix(qoi.by.sub(x2,tt=treat,ww=weights,qq=object$psclass)$q.table)
            nn <- c(nn,paste(names(covariates)[i],names(covariates)[j],sep="x"))
          }
        }
  }   
  qn <- aa[[1]]$qn
  dimnames(q.table) <- list(nn,row.names(aa[[i]]$q.table),paste("Subclass",1:qbins))
  obj$verbose <- verbose
  obj$sum.all <- sum.all
  obj$sum.matched <- sum.matched
  obj$q.table <- q.table
  obj$qn <- qn
  class(obj) <- "summary.matchit.subclass"
  obj
}
