data <- read.table("ajps2002-full1b.txt", header=T)
data$treat <- data$demsnmaj

sink("matchfdaW.out")
library(Matchit)
mout <- matchit(treat ~ orderent + stafcder + prevgenx + lethal +
                deathrt1 + hosp01 + I(hospdisc/1000000) +
                hhosleng + femdiz01 + mandiz01 + peddiz01 + acutediz +
                orphdum + natreg + I(natregsq/1000) +  wpnoavg3 +
                sqrt(wpnoavg3) + vandavg3 + condavg3, data=data,
                discard=1)   

print(summary(mout))
mdata <- match.data(mout)

library(Zelig)
fullmodel <- as.formula(Surv(acttime, d) ~ treat + prevgenx + lethal +
                        deathrt1 + acutediz + hosp01 + I(hospdisc/10000) +
                        hhosleng + femdiz01 + mandiz01 + peddiz01 +
                        orphdum + natreg + I(natregsq/1000) + vandavg3 + 
                        wpnoavg3 + condavg3 +
                        orderent + stafcder)

res <- zelig(fullmodel, data=data, model="weibull")
mres <- zelig(fullmodel, data=mdata, model="weibull")
xvars <- names(res$coefficients)
xvars <- xvars[3:length(xvars)]
start <- paste("Surv(acttime, d) ~ treat")


library(combinat)
N <- 1
total <- 1
for (i in N:(length(xvars)-1))
  total <- total + ncol(combn(xvars, i)) 
cat("\n I'm going to run", total, "regressions!\n")
    
ate <- mate <- rep(0,total)
counter <- 1
cat("start", date(), "\n")
for (i in N:(length(xvars)-1)) {
  allsubset <- combn(xvars, i)
  for (j in 1:ncol(allsubset)) {
    ftmp <- start
    for (k in 1:i) 
      ftmp <- paste(ftmp, "+", allsubset[k,j])
    ftmp <- as.formula(ftmp)
    tmp <- survreg(ftmp,  data = data, dist="weibull")
    x0 <- x1 <- model.matrix(tmp)
    x0[,"treat"] <- 0
    x1[,"treat"] <- 1
    ate[counter] <- mean(exp(x1%*%tmp$coefficients)*gamma(1+1/tmp$scale) -
                         exp(x0%*%tmp$coefficients)*gamma(1+1/tmp$scale))
    tmp <- survreg(ftmp, data = mdata, dist="weibull") 
    x0 <- x1 <- model.matrix(tmp)
    x0[,"treat"] <- 0
    x1[,"treat"] <- 1
    mate[counter] <- mean(exp(x1%*%tmp$coefficients)*gamma(1+1/tmp$scale) -
                          exp(x0%*%tmp$coefficients)*gamma(1+1/tmp$scale))
    counter <- counter + 1
  }
  cat(i,"covariates:",date(),"\n")
}

## using the full model
x0 <- x1 <- model.matrix(res)
x0[,"treat"] <- 0
x1[,"treat"] <- 1
ate[length(ate)] <- mean(exp(x1%*%res$coefficients)*gamma(1+1/tmp$scale) -
                         exp(x0%*%res$coefficients)*gamma(1+1/tmp$scale))
sims <- 5000
library(MASS)
param <- mvrnorm(sims, mu=c(res$coefficients, log(res$scale)),
                 Sigma=res$var) 
sate <- exp(param[,-ncol(param)]%*%t(x1))*gamma(1+1/exp(param[,ncol(param)])) -
  exp(param[,-ncol(param)]%*%t(x0))*gamma(1+1/exp(param[,ncol(param)]))

est <- c(mean(apply(sate, 1, mean)), sd(apply(sate, 1, mean)))

x0 <- x1 <- model.matrix(mres)
x0[,"treat"] <- 0
x1[,"treat"] <- 1
mate[length(mate)] <- mean(exp(x1%*%mres$coefficients)*gamma(1+1/tmp$scale) -
                           exp(x0%*%mres$coefficients)*gamma(1+1/tmp$scale))

mparam <- mvrnorm(sims, mu=c(mres$coefficients, log(mres$scale)),
                  Sigma=mres$var) 
msate <- exp(mparam[,-ncol(mparam)]%*%t(x1))*gamma(1+1/exp(mparam[,ncol(mparam)])) -
  exp(mparam[,-ncol(mparam)]%*%t(x0))*gamma(1+1/exp(mparam[,ncol(mparam)]))

mest <- c(mean(apply(msate, 1, mean)), sd(apply(msate, 1, mean)))

save.image("matchfdaW.RData")
sink()
