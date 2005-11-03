## Zelig: 2.4-6, MatchIt: 2.1-4
## Carpenter's Original Model
##

rm(list=ls())
## libraries
library(Zelig)
library(MatchIt)
library(combinat)

## inputs
analysis <- TRUE
model <- "lognorm"
qoi <- "ATT"
insample <- TRUE
se <- TRUE
file <- "matchfdaCI.RData"
sims <- 0

## data
data <- read.table("ajps2002-full1b.txt", header=T)

## democratic senate as the treatment variable
data$treat <- data$demsnmaj

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

## democratic senate
spec <- treat ~ orderent  + prevgenx  + lethal +
  deathrt1 + hosp01 + hospdisc + hhosleng + femdiz01 + mandiz01 +
  peddiz01 + acutediz + orphdum + natreg  +  wpnoavg3 +
  vandavg3 + condavg3 + stafcder + sqrt(hospdisc) 

m.out <- matchit(spec, method = "nearest", data = data,
                 discard="both", exact=c("acutediz", "peddiz01", "hosp01",
                                   "lethal", "femdiz01", "mandiz01"))

print(summary(m.out))
print(summary(m.out, standardize=TRUE))
mdata <- match.data(m.out)

#m.out1 <- matchit(spec, method = "genetic", wait.generations = 1000,
#                  data = data, discard="both") 

#sink("matchfda.out")
#print(summary(m.out))
#print(summary(m.out1))
#sink()

## functions to calculate QOI
## only lognormal with insample ATT has variance
qoical <- function(object, model, what = "ATT", insample = FALSE, se=FALSE) {
  if (se) {
    if (model != "lognorm")
      stop("se is not available for this model.")
    if (what != "ATT")
    if (!insample)
      stop("se is not available for this qoi.")
  }
  x <- model.matrix(object)
  y <- object$y
  beta <- coef(object)
  scale <- object$scale
  if (what == "ATE") {
    x0 <- x1 <- x
    y0 <- y1 <- y
    x0[,"treat"] <- 0
    x1[,"treat"] <- 1
    y1[x[,"treat"] == 0,2] <- 0
    y0[x[,"treat"] == 1,2] <- 0
  }
  if (what == "ATT") {
    x0 <- x1 <- x[x[,"treat"] == 1,] # treated only
    x0[,"treat"] <- 0  # counterfactual for the treated
    y <- y[x[,"treat"] == 1,] # observed Y(1) for the treated
  }
  if (insample) {  ## insample qoi
    if (what == "ATT") {
      if (model == "lognorm") {
        res <- mean(ifelse(y[,2]>0.5, exp(y[,1]),
                           exp(x1%*%beta+0.5*scale^2)) -
                    exp(x0%*%beta+0.5*scale^2))
        if (se) {
          delta1 <- (y[,2]<0.5)*cbind(c(exp(x1%*%beta+0.5*scale^2))*x1,
                                      c(exp(x1%*%beta+0.5*scale^2)*scale^2))
          delta0 <- cbind(c(exp(x0%*%beta+0.5*scale^2))*x0,
                          c(exp(x0%*%beta+0.5*scale^2)*scale^2))
          grad <- apply(delta1-delta0, 2, mean)
          var <- t(grad)%*%vcov(object)%*%grad
          res <- c(res, var, 2*1.96*sqrt(var))
        }
      }
      else if (model == "exp")
        res <- mean(ifelse(y[,2]>0.5, exp(y[,1]),
                           exp(x1%*%beta))-
                    exp(x0%*%beta))
      else if (model == "weibull")
        res <- mean(ifelse(y[,2]>0.5, exp(y[,1]),
                           exp(x1%*%beta)*gamma(1+scale))- 
                    exp(x0%*%beta)*gamma(1+scale))
      else
        stop("invalid model.")
    }
    if (what == "ATE") {
      if (model == "lognorm")
        res <- mean(ifelse(y1[,2]>0.5, exp(y1[,1]),
                           exp(x1%*%beta+0.5*scale^2)) -
                    ifelse(y0[,2]>0.5, exp(y0[,1]),
                           exp(x0%*%beta+0.5*scale^2)))
      else if (model == "exp")
        res <- mean(ifelse(y1[,2]>0.5, exp(y1[,1]),
                           exp(x1%*%beta)) -
                    ifelse(y0[,2]>0.5, exp(y0[,1]),
                           exp(x0%*%beta)))
      else if (model == "weibull")
        res <- mean(ifelse(y1[,2]>0.5, exp(y1[,1]),
                    exp(x1%*%beta)*gamma(1+scale)) -
                    ifelse(y0[,2]>0.5, y0[,1],
                           exp(x0%*%beta)*gamma(1+scale)))
      else
        stop("invalid model.")
    }
  }
  else { ## population qoi
    if (model == "lognorm")
      res <- mean(exp(x1%*%beta+0.5*scale^2) -
                  exp(x0%*%beta+0.5*scale^2))
    else if (model == "exp")
      res <- mean(exp(x1%*%beta) -
                  exp(x0%*%beta))
    else if (model == "weibull")
      res <- mean(exp(x1%*%beta)*gamma(1+scale) -
                  exp(x0%*%beta)*gamma(1+scale))
    else
      stop("invalid model.")
  }  
  return(res)
}

if(analysis) {
  ## the original full model
  fullmodel <- as.formula(Surv(acttime, d) ~ treat + orderent  +
                          prevgenx  + lethal + 
                          deathrt1 + hosp01 + hospdisc + hhosleng +
                          femdiz01 + mandiz01 + peddiz01 + acutediz +
                          orphdum + natreg + wpnoavg3 + condavg3 +
                          vandavg3 + stafcder + I(natreg^2)
                          )
  res <- zelig(fullmodel, data=data, model=model)
  print(summary(res))
  xvars <- names(res$coefficients)
  xvars <- xvars[3:length(xvars)]
  ## the simplest model
  start <- paste("Surv(acttime, d) ~ treat")

  ## counting the number of regressions to run
  total <- rep(0,length(xvars))
  for (i in 1:length(xvars)) 
    total[i] <- nCm(length(xvars), i)
  if (sims > 0)
    run <- round(total/sum(total)*sims)
  else
    run <- total
  if (se)
    ate <- mate <- matrix(NA, nrow=sum(run)+1, ncol=3)
  else
    ate <- mate <- rep(NA, sum(run)+1)
  counter <- 1
  
  ## original data
  tmp <- zelig(as.formula(start), data = data, model=model)
  if (se)
    ate[counter,] <- qoical(tmp, model, qoi, insample, se)
  else
    ate[counter] <- qoical(tmp, model, qoi, insample, se)    
  ## matched data
  tmp <- zelig(as.formula(start), data = mdata, model=model)
  if (se)
    mate[counter,] <- qoical(tmp, model, qoi, insample, se)
  else
    mate[counter] <- qoical(tmp, model, qoi, insample, se)
    
  ## looping
  counter <- 2
  cat("start", date(), "\n")
  for (i in 1:length(xvars)) {
    if (run[i]>0) {
      for (j in 1:run[i]) {
        ftmp <- start
        bk <- sample(xvars,i,replace=F)
        ftmp <- as.formula(paste(ftmp, "+",paste(bk,collapse=" + ")))
        ## original data
        tmp <- zelig(ftmp,  data = data, model=model)
        if (se)
          ate[counter,] <- qoical(tmp, model, qoi, insample, se)
        else
          ate[counter] <- qoical(tmp, model, qoi, insample, se)          
        ## matched data
        xtmp <- model.matrix(ftmp, model.frame(terms(ftmp), mdata))
        if (qr(xtmp)$rank >= ncol(xtmp)) {
          tmp <- zelig(ftmp,  data = mdata, model=model)
          if (se)
            mate[counter,] <- qoical(tmp, model, qoi, insample, se)
          else
            mate[counter] <- qoical(tmp, model, qoi, insample, se)
        }
        else
          if (se)
            mate[counter,] <- mate[counter-1,]
          else
            mate[counter] <- mate[counter-1]
        counter <- counter + 1
      }
      cat(i,"covariates:",date(),"\n")
    }
  }
}
save.image(file, compress=TRUE)
