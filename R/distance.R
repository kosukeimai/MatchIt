distance <- function(formula, model="logit", data,
                     discard=0,reestimate=FALSE,counter=TRUE,...){
  if(counter){cat("Calculating propensity score...")}
  mf<-match.call()
  if(!is.vector(discard) | !identical(round(length(discard)),1)){
    warning("discard=", mf$discard,
            " is invalid; used discard=0 instead", call.=FALSE)
    discard <- 0}
  if(!(identical(reestimate,TRUE) | identical(reestimate,FALSE))){
    warning("reestimate=",
            mf$reestimate," is invalid; used reestimate=FALSE instead",
            call.=FALSE)
    mf$reestimate=FALSE}
  mf$counter <- mf$discard<-mf$model<-mf$reestimate<-NULL
  dotsub <- eval(mf$subset, sys.frame(sys.parent()))
  a <- 0 
  while(a<2){
    if(a==1)
      if(is.null(dotsub))
        mf$subset<-in.sample
      else{
        dotsub[dotsub] <- in.sample
        mf$subset <- dotsub
      }
    if(model=="logit" || model=="probit"){
      mf[[1]]<-as.name("glm")
      if(model=="logit"){mf$family<-as.name("binomial")} else
      {mf$family <- binomial(link="probit")}  #is there a better way to define the link function here?
      res<-eval(as.call(mf), sys.frame(sys.parent()))
      if(a==0){
        pscore<-fitted(res)
        treat<-res$y
        covariates <- res$model[,-1,drop=F]
      }
      else{
        pscore[!in.sample] <- NA
        pscore[in.sample] <- fitted(res)
      }
      assign.model<-summary(res)
    }
    else if(model=="nnet"){
      require(nnet)
      mf[[1]] <- as.name("nnet.formula")
      mf$size <- 3  #is this the right default setting?
      res <- eval(as.call(mf),sys.frame(sys.parent()))
      if(a==0){
        pscore <- as.vector(fitted(res))
        nnames <- dimnames(fitted(res))[[1]]
        names(pscore) <- nnames
        if(is.null(mf$data)){
          treat <- model.frame(formula(res))[,1]
          covariates <- model.frame(delete.response(terms(formula(res))),drop=F)
        }else{
          treat <- model.frame(formula(res),data)[,1]
          covariates <- model.frame(delete.response(terms(formula(res))),data,drop=F)
        }
        if(!is.null(dotsub)){
          treat <- treat[dotsub]
          covariates <- covariates[dotsub,,drop=F]
        }
        names(treat) <- nnames
        assign.model <- summary(res)
      }
      else{
        pscore[!in.sample] <- NA
        pscore[in.sample] <- as.vector(fitted(res))
        assign.model <- res
        dan <- mf
      }
    }    
    else if(model=="GAM"){
      require(mgcv)
      mf[[1]] <- as.name("gam")
      mf$family<-as.name("binomial")
      res <- eval(as.call(mf),sys.frame(sys.parent()))
      if(a==0){
        treat <- res$y
        names(treat) <- row.names(data)  #assuming we have a data frame
        pscore <- fitted(res)
        names(pscore) <- names(treat)
        covariates <- res$model[,2:ncol(res$model),drop=F]
      }
      else{
        pscore[!in.sample] <- NA
        pscore[in.sample] <- fitted(res)
      }
      assign.model <- summary(res)
    }
    else if(model=="cart"){
      require(rpart)
      mf[[1]] <- as.name("rpart")
      res <- eval(as.call(mf),sys.frame(sys.parent()))
      if(a==0){
        pscore <- predict(res)
        treat <- res$y
        if(!is.null(dotsub)){
          covariates <- data[dotsub,attr(delete.response(terms(formula)),"term.labels"),drop=F]
        } else {covariates <- data[,attr(delete.response(terms(formula)),"term.labels"),drop=F]}
      }
      else{
        pscore[!in.sample] <- NA
        pscore[in.sample] <- predict(res)
      }
      assign.model <- print(res)
    }
    else {stop("Must specify valid model to estimate propensity score",call.=FALSE)}
    if(a==0){
      maxp0 <- max(pscore[treat==0])
      minp0 <- min(pscore[treat==0])
      maxp1 <- max(pscore[treat==1])
      minp1 <- min(pscore[treat==1])}
    if(discard==0){
      in.sample <- !is.na(treat) #just to make all TRUE
      a <- 2}
    else if(discard==1){
      if(a==0){in.sample <- (pscore>=max(minp1,minp0) & pscore<=min(maxp0,maxp1))}
      if(reestimate){(a <- a+1)}
      else {a <- 2}
    }
    else if(discard==2){
      if(a==0){in.sample <- (pscore>=minp1 & pscore<=maxp1)}
      if(reestimate){a <- a+1}
      else {a <- 2}
    }
    else if(discard==3){
      if(a==0){in.sample <- (pscore>=minp0 & pscore<=maxp0)}
      if(reestimate){a <- a+1}
      else {a <- 2}
    }
    else(stop("Invalid discard option",call.=FALSE))
  }
  ## return
  if(counter){cat("Done\n")}
  list(in.sample=in.sample, pscore=pscore, treat=treat, covariates=covariates,
                     assign.model=assign.model)
}
