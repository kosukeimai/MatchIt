summary.matchit.subclass <- function(object, interactions = FALSE,
                                     addlvariables=NULL, standardize = FALSE,
                                     ...) {

  XX <- cbind(distance=object$distance,object$X)
  if (!is.null(addlvariables)) XX <- cbind(XX, addlvariables)

  treat <- object$treat
  weights <- object$weights
  nam <- dimnames(XX)[[2]]
  kk <- ncol(XX)

  ## Summary Stats
  aa <- apply(XX,2,qoi,tt=treat,ww=as.numeric(weights!=0),standardize=standardize)
  sum.all <- as.data.frame(matrix(0,kk,6))
  sum.matched <- as.data.frame(matrix(0,kk,6))
  row.names(sum.all) <- row.names(sum.matched) <- nam
  names(sum.all) <- names(sum.matched) <- names(aa[[1]])
  sum.all.int <- sum.matched.int <- NULL
  for(i in 1:kk){
    sum.all[i,] <- aa[[i]][1,]
    sum.matched[i,] <- aa[[i]][2,]
    if(interactions){
      for(j in i:kk){
        x2 <- XX[,i]*as.matrix(XX[,j])
        jqoi <- qoi(x2,tt=treat,ww=as.numeric(weights!=0),standardize=standardize)
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
  qbins <- max(object$subclass,na.rm=TRUE)
  if(interactions){
    q.table <- array(0,dim=c(kk+sum(1:kk),7,qbins))
    ii <- 0
    nn <- NULL
  } else {
    q.table <- array(0,dim=c(kk,7,qbins))
  }
  aa <- apply(XX,2,qoi.by.sub,tt=treat,ww=weights,
                  qq=object$subclass,standardize=standardize)
  for(i in 1:kk){
    if(!interactions){
      q.table[i,,] <- as.matrix(aa[[i]]$q.table)
      nn <- names(aa)
        } else {
          ii <- ii + 1 
          q.table[ii,,] <- as.matrix(aa[[i]]$q.table)
          nn <- c(nn,names(aa)[i])
          for(j in i:kk){
            ii <- ii + 1 
            x2 <- XX[,i]*as.matrix(XX[,j])
            q.table[ii,,] <- as.matrix(qoi.by.sub(x2,tt=treat,ww=weights,qq=object$subclass,standardize=standardize)$q.table)
            nn <- c(nn,paste(nam[i],nam[j],sep="x"))
          }
        }
  }   
  qn <- aa[[1]]$qn
  dimnames(q.table) <- list(nn,row.names(aa[[i]]$q.table),paste("Subclass",1:qbins))
  
  ## Aggregate Subclass 
  if(is.null(object$call$sub.by)){
    object$call$sub.by <- "treat"
  }
  if(object$call$sub.by=="treat") {
    wsub <- qn[1,]/sum(qn[1,])
  } else if(sub.by=="control") {
    wsub <- qn[2,]/sum(qn[2,])
  } else if(sub.by=="all") {
    wsub <- qn[3,]/sum(qn[3,])
  }
  sum.subclass <- sum.all
  for(i in 1:kk){
    for(j in 1:7){
      if(j==3) {
        sum.subclass[i,j] <- sqrt(sum((wsub^2)*(q.table[i,j,]^2)))
      } else {
        sum.subclass[i,j] <- sum(wsub*q.table[i,j,])
      }
    }
  }

  ## Imbalance Reduction
  stat0 <- abs(cbind(sum.all[,2]-sum.all[,1],
                     sum.all[,4:6]))
  stat1 <- abs(cbind(sum.subclass[,2]-sum.subclass[,1],
                     sum.subclass[,4:6]))
  reduction <- as.data.frame(100*(stat0-stat1)/stat0)
  if(sum(stat0==0 & stat1==0, na.rm=T)>0){
    reduction[stat0==0 & stat1==0] <- 0
  }
  if(sum(stat0==0 & stat1>0,na.rm=T)>0){
    reduction[stat0==0 & stat1>0] <- -Inf
  }
  if (standardize)
    names(reduction) <- c("Std. Mean Diff.", "eCDF Med","eCDF Mean",
                          "eCDF Max")
  else
    names(reduction) <- c("Mean Diff.", "eQQ Med","eQQ Mean",
                          "eQQ Max")
  ## output
  res <- list(call=object$call, sum.all = sum.all, sum.matched = sum.matched,
              sum.subclass = sum.subclass, reduction = reduction,
              qn = qn, q.table = q.table)
  class(res) <- c("summary.matchit.subclass", "summary.matchit")
  return(res)
}
