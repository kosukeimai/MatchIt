data <- read.table("ajps2002-full1b.txt", header=T)
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

sink("matchfda.out")
library(MatchIt)
mout <- matchit(treat ~ orderent + stafcder + prevgenx + lethal +
                deathrt1 + hosp01 + hospdisc +
                hhosleng + femdiz01 + mandiz01 + peddiz01 + acutediz +
                orphdum + natreg + I(natreg^2) +  wpnoavg3 +
                sqrt(wpnoavg3) + vandavg3 + condavg3, data=data,
                discard=1)   

print(summary(mout))
mdata <- match.data(mout)

library(Zelig)
fullmodel <- as.formula(Surv(acttime, d) ~ treat + prevgenx + lethal +
                        deathrt1 + acutediz + hosp01 + hospdisc +
                        hhosleng + femdiz01 + mandiz01 + peddiz01 +
                        orphdum + natreg + I(natreg^2) + vandavg3 + 
                        wpnoavg3 + condavg3 + orderent + stafcder +
                        I(prevgenx^2) + I(deathrt1^2) + I(hospdisc^2) + 
                        I(hhosleng^2) + I(vandavg3^2) + I(wpnoavg3^2) +
                        I(condavg3^2) + I(orderent^2) + I(stafcder^2)
                        )

res <- zelig(fullmodel, data=data, model="lognorm")
print(summary(res))
mres <- zelig(fullmodel, data=mdata, model="lognorm")
print(summary(res))
xvars <- names(res$coefficients)
xvars <- xvars[3:length(xvars)]
start <- paste("Surv(acttime, d) ~ treat")


library(combinat)
sims <- 1000000
N <- 1
total <- rep(0,length(xvars)-1)
for (i in N:(length(xvars)-1)) 
  total[i] <- nCm(length(xvars), i) 
run <- round(total/sum(total)*sims)

ate <- mate <- rep(0,sum(run))
counter <- 1
cat("start", date(), "\n")
for (i in N:(length(xvars)-1)) {
  if (run[i]>0) {
    for (j in 1:ncol(run[i])) {
      ftmp <- start
      bk <- sample(xvars,i,replace=F)
      ftmp <- as.formula(paste(ftmp, "+",paste(bk,collapse=" + ")))
      tmp <- survreg(ftmp,  data = data, dist="lognormal")
      x0 <- x1 <- model.matrix(tmp)
      x0[,"treat"] <- 0
      x1[,"treat"] <- 1
      ate[counter] <- mean(exp(x1%*%tmp$coefficients+0.5*tmp$scale^2) -
                           exp(x0%*%tmp$coefficients+0.5*tmp$scale^2))
      tmp <- survreg(ftmp, data = mdata, dist="lognormal") 
      x0 <- x1 <- model.matrix(tmp)
      x0[,"treat"] <- 0
      x1[,"treat"] <- 1
      mate[counter] <- mean(exp(x1%*%tmp$coefficients+0.5*tmp$scale^2) -
                            exp(x0%*%tmp$coefficients+0.5*tmp$scale^2))
      counter <- counter + 1
    }
    cat(i,"covariates:",date(),"\n")
  }
}

## using the full model
x0 <- x1 <- model.matrix(res)
x0[,"treat"] <- 0
x1[,"treat"] <- 1
ate[length(ate)] <- mean(exp(x1%*%res$coefficients+0.5*tmp$scale^2) -
                         exp(x0%*%res$coefficients+0.5*tmp$scale^2))
sims <- 5000
library(MASS)
param <- mvrnorm(sims, mu=c(res$coefficients, log(res$scale)),
                 Sigma=res$var) 
sate <- exp(param[,-ncol(param)]%*%t(x1) +
            0.5*exp(param[,ncol(param)])^2) -
  exp(param[,-ncol(param)]%*%t(x0) +
      0.5*exp(param[,ncol(param)])^2)

est <- c(mean(apply(sate, 1, mean)), sd(apply(sate, 1, mean)))

x0 <- x1 <- model.matrix(mres)
x0[,"treat"] <- 0
x1[,"treat"] <- 1
mate[length(mate)] <- mean(exp(x1%*%mres$coefficients+0.5*tmp$scale^2) -
                           exp(x0%*%mres$coefficients+0.5*tmp$scale^2))

mparam <- mvrnorm(sims, mu=c(mres$coefficients, log(mres$scale)),
                  Sigma=mres$var) 
msate <- exp(mparam[,-ncol(mparam)]%*%t(x1) +
             0.5*exp(mparam[,ncol(mparam)])^2) -
  exp(mparam[,-ncol(mparam)]%*%t(x0) +
      0.5*exp(mparam[,ncol(mparam)])^2)

mest <- c(mean(apply(msate, 1, mean)), sd(apply(msate, 1, mean)))

save.image("matchfda.RData")
sink()
