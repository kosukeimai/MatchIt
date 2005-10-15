##
## Carpenter's Model 1 without all the other oversight variables
## 

rm(list=ls())
data <- read.table("ajps2002-full1b.txt", header=T)
data$treat <- data$demsnmaj

analysis <- TRUE
model <- "weibull"

## rescaling
data$hospdisc <- data$hospdisc/100000
data$natreg <- data$natreg/100
data$stafcder <- data$stafcder/100
data$prevgenx <- data$prevgenx/100
data$hhosleng <- data$hhosleng/10
data$condavg3 <- data$condavg3/10
data$orderent <- data$orderent/10
data$vandavg3 <- data$vandavg3/10
data$wpnoavg3 <- data$wpnoavg3/100

library(MatchIt)
spec <- treat ~ orderent  + prevgenx  + lethal +
  deathrt1 + hosp01 + hospdisc + hhosleng + femdiz01 + mandiz01 +
  peddiz01 + acutediz + orphdum + natreg  +  wpnoavg3 +
  vandavg3 + condavg3 + stafcder + sqrt(hospdisc) 

m.out <- matchit(spec, method = "nearest", data = data,
                 discard="both", exact=c("acutediz", "peddiz01", "hosp01",
                                   "lethal", "femdiz01", "mandiz01"))
summary(m.out)

#sink("matchfda.out")
#print(summary(m.out))
#sink()

mdata <- match.data(m.out)

if(analysis) {
  library(Zelig)
  ## the original full model
  fullmodel <- as.formula(Surv(acttime, d) ~ treat + orderent  +
                          prevgenx  + lethal + 
                          deathrt1 + hosp01 + hospdisc + hhosleng +
                          femdiz01 + mandiz01 + peddiz01 + acutediz +
                          orphdum + natreg  +  wpnoavg3 +
                          vandavg3 + condavg3 + stafcder
                          )
  ## full model for matched data
  mfullmodel <- as.formula(Surv(acttime, d) ~ treat + orderent  +
                           prevgenx + deathrt1 + hospdisc + hhosleng +
                           orphdum + natreg  +  wpnoavg3 +
                           vandavg3 + condavg3 + stafcder 
                           )
  
  res <- zelig(fullmodel, data=data, model=model)
  print(summary(res))
  mres <- zelig(mfullmodel, data=mdata, model=model)
  print(summary(mres))
  xvars <- names(res$coefficients)
  xvars <- xvars[3:length(xvars)]
  mxvars <- names(mres$coefficients)
  mxvars <- mxvars[3:length(mxvars)]
  start <- paste("Surv(acttime, d) ~ treat")
  

  library(combinat)

  ## original data
  N <- 1
  total <- 1
  for (i in N:(length(xvars)-1)) 
    total <- total + ncol(combn(xvars, i))
  
  ate <- rep(0,total)
  counter <- 1
  cat("start", date(), "\n")
  for (i in N:(length(xvars)-1)) {
    allsubset <- combn(xvars, i)
    for (j in 1:ncol(allsubset)) {
      ftmp <- start
      for (k in 1:i)
        ftmp <- paste(ftmp, "+",allsubset[k,j])
      ftmp <- as.formula(ftmp)
      tmp <- zelig(ftmp,  data = data, model=model)
      x0 <- x1 <- model.matrix(tmp)
      x0[,"treat"] <- 0
      x1[,"treat"] <- 1
      if (model == "lognorm")
        ate[counter] <- mean(exp(x1%*%tmp$coefficients+0.5*tmp$scale^2) -
                             exp(x0%*%tmp$coefficients+0.5*tmp$scale^2))
      else if (model == "exp")
        ate[counter] <- mean(exp(x1%*%tmp$coefficients) -
                             exp(x0%*%tmp$coefficients))
      else if (model == "weibull")
        ate[counter] <- mean(exp(x1%*%tmp$coefficients)*gamma(1+tmp$scale) -
                             exp(x0%*%tmp$coefficients)*gamma(1+tmp$scale))
      else
        stop("invalid model")
      counter <- counter + 1
    }
    cat(i,"covariates:",date(),"\n")
  }
  
  ## matched data
  N <- 1
  total <- 1
  for (i in N:(length(mxvars)-1))
    total <- total + ncol(combn(mxvars, i))
  
  mate <- rep(0,total)
  counter <- 1
  cat("start", date(), "\n")
  for (i in N:(length(mxvars)-1)) {
    allsubset <- combn(mxvars, i)
    for (j in 1:ncol(allsubset)) {
      ftmp <- start
      for (k in 1:i)
        ftmp <- paste(ftmp, "+", allsubset[k,j])
      ftmp <- as.formula(ftmp)
      tmp <- zelig(ftmp,  data = mdata, model=model)
      x0 <- x1 <- model.matrix(tmp)
      x0[,"treat"] <- 0
      x1[,"treat"] <- 1
      if (model == "lognorm")
        mate[counter] <- mean(exp(x1%*%tmp$coefficients+0.5*tmp$scale^2) -
                              exp(x0%*%tmp$coefficients+0.5*tmp$scale^2))
      else if (model == "exp")
        mate[counter] <- mean(exp(x1%*%tmp$coefficients) -
                              exp(x0%*%tmp$coefficients))
      else if (model == "weibull")
        mate[counter] <- mean(exp(x1%*%tmp$coefficients)*gamma(1+tmp$scale) -
                              exp(x0%*%tmp$coefficients)*gamma(1+tmp$scale))
      else
        stop("invalid model")
      counter <- counter + 1
    }
    cat(i,"covariates:",date(),"\n")
  }
  
  sims <- 5000
  library(MASS)

  ## using the full model
  x0 <- x1 <- model.matrix(res)
  x0[,"treat"] <- 0
  x1[,"treat"] <- 1
  if (model == "lognorm") {
    ate[length(ate)] <- mean(exp(x1%*%res$coefficients+0.5*res$scale^2) -
                             exp(x0%*%res$coefficients+0.5*res$scale^2))
    param <- mvrnorm(sims, mu=c(res$coefficients, log(res$scale)),
                     Sigma=vcov(res)) 
    sate <- exp(param[,-ncol(param)]%*%t(x1) +
                0.5*exp(param[,ncol(param)])^2) -
                  exp(param[,-ncol(param)]%*%t(x0) +
                      0.5*exp(param[,ncol(param)])^2)
  }
  else if (model == "exp") {
    ate[length(ate)] <- mean(exp(x1%*%res$coefficients) -
                             exp(x0%*%res$coefficients))
    param <- mvrnorm(sims, mu=res$coefficients, Sigma=vcov(res)) 
    sate <- exp(param%*%t(x1)) - exp(param%*%t(x0)) 
  }
  else if (model == "weibull") {
    ate[length(ate)] <- mean(exp(x1%*%res$coefficients)*gamma(1+res$scale) -
                             exp(x0%*%res$coefficients)*gamma(1+res$scale))
    sate <- exp(param[,-ncol(param)]%*%t(x1))*gamma(1+exp(param[,ncol(param)]))-
      exp(param[,-ncol(param)]%*%t(x0))*gamma(1+exp(param[,ncol(param)])) 
  }
  else
    stop("invalid model")
  est <- c(mean(apply(sate, 1, mean)), sd(apply(sate, 1, mean)))

  ## matched data
  x0 <- x1 <- model.matrix(mres)
  x0[,"treat"] <- 0
  x1[,"treat"] <- 1
  if (model == "lognorm") {
    mate[length(mate)] <- mean(exp(x1%*%mres$coefficients+0.5*mres$scale^2) -
                             exp(x0%*%mres$coefficients+0.5*mres$scale^2))
    param <- mvrnorm(sims, mu=c(mres$coefficients, log(mres$scale)),
                     Sigma=vcov(mres)) 
    msate <- exp(param[,-ncol(param)]%*%t(x1) +
                0.5*exp(param[,ncol(param)])^2) -
                  exp(param[,-ncol(param)]%*%t(x0) +
                      0.5*exp(param[,ncol(param)])^2)
  }
  else if (model == "exp") {
    mate[length(mate)] <- mean(exp(x1%*%mres$coefficients) -
                             exp(x0%*%mres$coefficients))
    param <- mvrnorm(sims, mu=mres$coefficients, Sigma=vcov(mres)) 
    msate <- exp(param%*%t(x1)) - exp(param%*%t(x0)) 
  }
  else if (model == "weibull") {
    mate[length(mate)] <- mean(exp(x1%*%mres$coefficients)*gamma(1+mres$scale) -
                             exp(x0%*%mres$coefficients)*gamma(1+mres$scale))
    msate <- exp(param[,-ncol(param)]%*%t(x1))*gamma(1+exp(param[,ncol(param)]))-
      exp(param[,-ncol(param)]%*%t(x0))*gamma(1+exp(param[,ncol(param)])) 
  }
  else
    stop("invalid model")
  
  mest <- c(mean(apply(msate, 1, mean)), sd(apply(msate, 1, mean)))

  save.image("matchfdaBaseWei.RData")
  sink()
}
