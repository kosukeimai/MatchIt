summary.matchit.full <- function(object, interactions = FALSE, addlvariables = NULL, numdraws = 5000, ...) {

  XX <- cbind(distance=object$distance,object$X)
  if (!is.null(addlvariables)) XX <- cbind(XX, addlvariables)

  treat <- object$treat
  weights <- object$weights
  nam <- dimnames(XX)[[2]]
  kk <- ncol(XX)

  ## Get samples of T and C units to send to qqplot
  t.plot <- sample(names(treat)[treat==1], numdraws/2, replace=TRUE, prob=weights[treat==1])
  c.plot <- sample(names(treat)[treat==0], numdraws/2, replace=TRUE, prob=weights[treat==0])

  ## Summary Stats
  aa <- apply(XX,2,qoi,tt=treat,ww=weights, t.plot=t.plot, c.plot=c.plot)
  sum.all <- as.data.frame(matrix(0,kk,7))
  sum.matched <- as.data.frame(matrix(0,kk,7))
  row.names(sum.all) <- row.names(sum.matched) <- nam
  names(sum.all) <- names(sum.matched) <- names(aa[[1]])
  sum.all.int <- sum.matched.int <- NULL
  for(i in 1:kk){
    sum.all[i,] <- aa[[i]][1,]
    sum.matched[i,] <- aa[[i]][2,]
    if(interactions){
      for(j in i:kk){
        x2 <- XX[,i]*as.matrix(XX[,j])
        jqoi <- qoi(x2,tt=treat,ww=weights, t.plot=t.plot, c.plot=c.plot)
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

  ## Imbalance Reduction
  stat0 <- abs(cbind(sum.all[,2]-sum.all[,1],
                     sum.all[,4:7]))
  stat1 <- abs(cbind(sum.matched[,2]-sum.matched[,1],
                     sum.matched[,4:7]))
  reduction <- as.data.frame(100*(stat0-stat1)/stat0)
  if(sum(stat0==0 & stat1==0)>0){
    reduction[stat0==0 & stat1==0] <- 0
  }
  if(sum(stat0==0 & stat1>0)>0){
    reduction[stat0==0 & stat1>0] <- -Inf
  }
  names(reduction) <- c("Mean","QQ Med","QQ Mean", "QQ Max", "Std. Bias")

  ## Sample sizes
  if(sum(object$weights==0)!=0) { 
	nn <- rbind(table(object$treat), table(object$weights!=0,object$treat)[2:1,]) 
  }
  if(sum(object$weights==0)==0) { 
	nn <- rbind(table(object$treat), table(object$weights!=0, object$treat), c(sum(object$weights[object$treat==0]==0,
	sum(object$weights[object$treat==1]==0)))) 
  }
  
  dimnames(nn) <- list(c("Full", "Matched", "Discarded"), c("Control", "Treated"))

  ## output
  res <- list(call=object$call, nn = nn, sum.all = sum.all, sum.matched = sum.matched,
              reduction = reduction)
  class(res) <- c("summary.matchit.full", "summary.matchit")
  return(res)
}
