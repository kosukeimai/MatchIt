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
file <- "matchfdaLin.RData"
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
mdata <- match.data(m.out)

##m.out1 <- matchit(spec, method = "genetic", wait.generations = 1000,
#                  data = data, discard="both") 

#sink("matchfda.out")
#print(summary(m.out))
#print(summary(m.out1))
#sink()

## functions to calculate QOI
qoical <- function(object, model, what = "ATT", insample = FALSE) {
  x <- model.matrix(object)
  y <- object$y
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
      if (model == "lognorm")
        res <- mean(ifelse(y[,2], exp(y[,1]),
                           exp(x1%*%object$coefficients+0.5*object$scale^2)) -
                    exp(x0%*%object$coefficients+0.5*object$scale^2))
      else if (model == "exp")
        res <- mean(ifelse(y[,2], y[,1], exp(x1%*%object$coefficients))-
                    exp(x0%*%object$coefficients))
      else if (model == "weibull")
        res <- mean(ifelse(y[,2], y[,1],
                           exp(x1%*%object$coefficients)*gamma(1+object$scale))- 
                    exp(x0%*%object$coefficients)*gamma(1+object$scale))
      else
        stop("invalid model.")
    }
    if (what == "ATE") {
      if (model == "lognorm")
        res <- mean(ifelse(y1[,2], exp(y1[,1]),
                           exp(x1%*%object$coefficients+0.5*object$scale^2)) -
                    ifelse(exp(y0[,2]), y0[,1],
                           exp(x0%*%object$coefficients+0.5*object$scale^2)))
      else if (model == "exp")
        res <- mean(ifelse(y1[,2], y1[,1], exp(x1%*%object$coefficients)) -
                    ifelse(y0[,2], y0[,1], exp(x0%*%object$coefficients)))
      else if (model == "weibull")
        res <- mean(ifelse(y1[,2], y1[,1],
                    exp(x1%*%object$coefficients)*gamma(1+object$scale)) -
                    ifelse(y0[,2], y0[,1],
                           exp(x0%*%object$coefficients)*gamma(1+object$scale)))
      else
        stop("invalid model.")
    }
  }
  else { ## population qoi
    if (model == "lognorm")
      res <- mean(exp(x1%*%object$coefficients+0.5*object$scale^2) -
                  exp(x0%*%object$coefficients+0.5*object$scale^2))
    else if (model == "exp")
      res <- mean(exp(x1%*%object$coefficients) -
                  exp(x0%*%object$coefficients))
    else if (model == "weibull")
      res <- mean(exp(x1%*%object$coefficients)*gamma(1+object$scale) -
                  exp(x0%*%object$coefficients)*gamma(1+object$scale))
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
                          vandavg3 + stafcder
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
  ate <- mate <- rep(NA,sum(run)+1)
  counter <- 1
  
  ## original data
  tmp <- zelig(as.formula(start), data = data, model=model)
  ate[counter] <- qoical(tmp, model, qoi, insample)

  ## matched data
  tmp <- zelig(as.formula(start), data = mdata, model=model)
  mate[counter] <- qoical(tmp, model, qoi, insample)

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
        ate[counter] <- qoical(tmp, model, qoi, insample)
        ## matched data
        xtmp <- model.matrix(ftmp, model.frame(terms(ftmp), mdata))
        if (qr(xtmp)$rank >= ncol(xtmp)) {
          tmp <- zelig(ftmp,  data = mdata, model=model)
          mate[counter] <- qoical(tmp, model, qoi, insample)
        }
        else
          mate[counter] <- mate[counter-1]
        counter <- counter + 1
      }
      cat(i,"covariates:",date(),"\n")
    }
  }
}
save.image(file)
