summary.matchit.full <- function(object, interactions = FALSE,
                                 addlvariables = NULL, numdraws =
                                 5000, standardize = FALSE,
                                 ...) {

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
  aa <- apply(XX,2,qoi,tt=treat,ww=weights, t.plot=t.plot,
              c.plot=c.plot, standardize=standardize)
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
	names(x2) <- names(XX[,1])
        jqoi <- qoi(x2,tt=treat,ww=weights, t.plot=t.plot,
                    c.plot=c.plot, standardize=standardize)
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
                     sum.all[,4:6]))
  stat1 <- abs(cbind(sum.matched[,2]-sum.matched[,1],
                     sum.matched[,4:6]))
  reduction <- as.data.frame(100*(stat0-stat1)/stat0)
  if(sum(stat0==0 & stat1==0, na.rm=T)>0){
    reduction[stat0==0 & stat1==0] <- 0
  }
  if(sum(stat0==0 & stat1>0,na.rm=T)>0){
    reduction[stat0==0 & stat1>0] <- -Inf
  }
  if (standardize)
    names(reduction) <- c("Std. Mean Diff.", "eCDF Med","eCDF Mean", "eCDF Max")
  else
    names(reduction) <- c("Mean Diff.", "eQQ Med","eQQ Mean", "eQQ Max")

  ## Sample sizes
  nn <- matrix(0, ncol=2, nrow=4)
  nn[1,] <- c(sum(object$treat==0), sum(object$treat==1))
  nn[2,] <- c(sum(object$treat==0 & object$weights>0), sum(object$treat==1 & object$weights>0))
  nn[3,] <- c(sum(object$treat==0 & object$weights==0 & object$discarded==0), sum(object$treat==1 & object$weights==0 & object$discarded==0))
  nn[4,] <- c(sum(object$treat==0 & object$weights==0 & object$discarded==1), sum(object$treat==1 & object$weights==0 & object$discarded==1))
   
   dimnames(nn) <- list(c("All","Matched","Unmatched","Discarded"),
                        c("Control","Treated"))

  ## output
  res <- list(call=object$call, nn = nn, sum.all = sum.all, sum.matched = sum.matched,
              reduction = reduction)
  class(res) <- c("summary.matchit.full", "summary.matchit")
  return(res)
}
