neyman.est<-function(Y, object){

  weighted.var<-function(x, w)
    sum(w*(x-weighted.mean(x, w))^2)/(sum(w)-1)

  treat<-object$treat
  W<-object$weights
  indx<-object$subclass
  if(is.null(indx)){
    m <- 1
    indx <- rep(1, length(Y))
  } else {
    indx[is.na(indx)] <- 0
    m<-max(indx)
    W <- as.numeric(W>0)
  }
  ate.est<-ate.var<-Ntrt<-Ncont<-matrix(NA, nrow=1, ncol=m)
  for(i in 1:m){
    Ntrt[i]<-length(Y[treat==1 & indx==i])
    Ncont[i]<-length(Y[treat==0 & indx==i])
    ate.est[i]<-weighted.mean(Y[treat==1 & indx==i],
                              W[treat==1 & indx==i]) -
                                weighted.mean(Y[treat==0 & indx==i],
                                              W[treat==0 & indx==i])
    ate.var[i]<-weighted.var(Y[treat==1 & indx==i],
                             W[treat==1 & indx==i])/sum(W[treat==1 & indx == i]) +
                               weighted.var(Y[treat==0 & indx==i],
                                            W[treat==0 & indx==i])/sum(W[treat==0 & indx == i])
  }
  if(m==1) {
    Ntrt=length(Y[treat==1 & W!=0])
    Ncont=length(Y[treat==0 & W!=0])
  } else {
    colnames(ate.est)<-c(1:m)
    for(i in 1:m)
      colnames(ate.est)[i]<-paste("Subclass", i)
  }
    list(ate.est=ate.est, ate.var=ate.var, Ntrt=Ntrt, Ncont=Ncont, bootstrap=NULL)
}
