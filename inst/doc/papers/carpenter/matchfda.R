data <- read.table("ajps2002-full1b.txt", header=T)
data$treat <- data$demsnmaj

sink("matchfda.out")
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

res <- zelig(fullmodel, data=data, model="lognorm")
mres <- zelig(fullmodel, data=mdata, model="lognorm")
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
    tmp <- survreg(ftmp,  data = data, dist="lognormal")
    x0 <- x1 <- model.matrix(tmp)
    x0[,"treat"] <- 0
    x1[,"treat"] <- 1
    ate[counter] <- mean(exp(x1%*%tmp$coefficient+0.5*tmp$scale^2) -
                         exp(x0%*%tmp$coefficient+0.5*tmp$scale^2))
    tmp <- survreg(ftmp, data = mdata, dist="lognormal") 
    x0 <- x1 <- model.matrix(tmp)
    x0[,"treat"] <- 0
    x1[,"treat"] <- 1
    mate[counter] <- mean(exp(x1%*%tmp$coefficient+0.5*tmp$scale^2) -
                          exp(x0%*%tmp$coefficient+0.5*tmp$scale^2))
    counter <- counter + 1
  }
  cat(i,"covariates:",date(),"\n")
}
x0 <- x1 <- model.matrix(res)
x0[,"treat"] <- 0
x1[,"treat"] <- 1
ate[length(ate)] <- mean(exp(x1%*%res$coefficient+0.5*tmp$scale^2) -
                         exp(x0%*%res$coefficient+0.5*tmp$scale^2))
x0 <- x1 <- model.matrix(mres)
x0[,"treat"] <- 0
x1[,"treat"] <- 1
mate[length(mate)] <- mean(exp(x1%*%mres$coefficient+0.5*tmp$scale^2) -
                           exp(x0%*%mres$coefficient+0.5*tmp$scale^2))

save.image("matchfda.RData")
sink()
