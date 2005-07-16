neyman.matchit.exactfull <-function(Y, object){

  weighted.var<-function(x, w)
    sum(w*(x-weighted.mean(x, w))^2)/(sum(w)-1)
  
  W <- object$weights
  treat <- object$treat
  
  ate <- weighted.mean(Y[treat==1], W[treat==1]) -
    weighted.mean(Y[treat==0], W[treat==0])
  se <- sqrt(weighted.var(Y[treat==1], W[treat==1])/sum(W[treat==1])
             + weighted.var(Y[treat==0], W[treat==0])/sum(W[treat==0]))
  numtrt <- sum(W[treat==1])
  numcont <- sum(W[treat==0])
  
  ressum<-list(ate=ate, se=se, Ntrt=numtrt, Ncont=numcont)
  class(ressum) <- "neyman"
  ressum
}
