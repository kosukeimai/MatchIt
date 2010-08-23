summary.matchit <- function(object, interactions = FALSE,
                            addlvariables = NULL, standardize = FALSE,
                            ...) {

  X <- object$X
  ## Fix X matrix so that it doesn't have any factors  
  varnames <- colnames(X)
  for(var in varnames) {
        if(is.factor(X[,var])) {
                tempX <- X[,!colnames(X)%in%c(var)]   
                form<-formula(substitute(~dummy-1,list(dummy=as.name(var))))
                X <- cbind(tempX, model.matrix(form, X))
        }
  }

  ## No distance output for pure Mahalanobis
  if("matchit.mahalanobis"%in%class(object)){
    XX <- X 
  } else{
    XX <- cbind(distance=object$distance,X)
  }
  if (!is.null(addlvariables)) XX <- cbind(XX, addlvariables)

  treat <- object$treat
  weights <- object$weights
  nam <- dimnames(XX)[[2]]
  dupnam <- duplicated(nam)
  if(sum(dupnam)>0){
    nam[dupnam] <- paste(nam[dupnam],".1",sep="")
  }
  kk <- ncol(XX)

  ## Summary Stats
  aa <- apply(XX,2,qoi,tt=treat,ww=weights,standardize=standardize,std=T)
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
        jqoi <- qoi(x2,tt=treat,ww=weights,standardize=standardize,std=T)
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
                     sum.all[,5:7]))
  stat1 <- abs(cbind(sum.matched[,2]-sum.matched[,1],
                     sum.matched[,5:7]))
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
    names(reduction) <- c("Mean Diff.", "eQQ Med","eQQ Mean", "eQQ Max")
    
  ## output
  res <- list(call=object$call, nn = object$nn, sum.all = sum.all,
              sum.matched = sum.matched, reduction = reduction)
  class(res) <- "summary.matchit"
  return(res)
}
