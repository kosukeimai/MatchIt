

lognorm.llik <- function(param, Y, X, C){
  sigma <- exp(param[1])
  beta <- param[2:length(param)]
  mu <- X%*%beta
  value <- sum((1-C)*(-log(sigma)+dnorm((log(Y)-mu)/sigma, log=TRUE) -
                      pnorm((log(Y)-mu)/sigma, lower.tail= FALSE,
                            log.p=TRUE)))
  value <- value + sum(pnorm((log(Y)-mu)/sigma, lower.tail= FALSE,
                             log.p=TRUE))  
#  value <-
#    sum((1-C)*(-log(sigma)+dnorm((log(Y)-mu)/sigma, log=TRUE)))+
#      sum(C*pnorm((log(Y)-mu)/sigma, lower.tail= FALSE, log.p=TRUE))
  return(-value)
}

lognorm.inv.gauss.llik <- function(param, Y, X, C){
  sigma <- exp(param[1])
  theta <- exp(param[2])
  beta <- param[3:length(param)]
  mu <- X%*%beta

  lnS <- log(1 - pnorm((log(Y)-mu)/sigma))
  lnS <- pnorm((log(Y)-mu)/sigma, lower.tail = FALSE, log = TRUE)
  lnh <- dnorm(log(Y), mean=mu, sd=sigma, log = TRUE) - lnS
  value <- sum(-sqrt(1-2*theta*lnS)/theta + (1-C)*lnh -
    0.5*(1-C)*log(1-2*theta*lnS)) 
  
  #value <- sum((1-C)*(-log(sigma) + dnorm((log(Y)-mu)/sigma, log=TRUE) -
  #                    pnorm((log(Y)-mu)/sigma, lower.tail= FALSE,
  #                          log.p=TRUE)) -
  #             0.5*log(1-2*theta*pnorm((log(Y)-mu)/sigma,
  #                                     lower.tail= FALSE, 
  #                                     log.p=TRUE)))
  #value <- value + sum((1-sqrt(1-2*theta*pnorm((log(Y)-mu)/sigma,
  #                                             lower.tail= FALSE, 
  #                                             log.p=TRUE)))/theta) 
  return(-value)
}


## testing the functions
library(survival)
TEST <- FALSE
if (TEST) {
  data(ovarian)
  res1 <- survreg(Surv(futime, fustat) ~ ecog.ps, data = ovarian,
                  dist='lognormal', scale=0) 
  print(summary(res1))
  
  start <- c(0, lm(log(futime) ~ ecog.ps, ovarian)$coef)
  res2 <- optim(start, lognorm.llik, Y = ovarian$futime,
                X=cbind(rep(1, nrow(ovarian)), ovarian$ecog.ps),
                C=(ovarian$fustat==0)*1, method="BFGS")  
  print(res2)
  
  start <- c(1, 1, lm(log(futime) ~ ecog.ps, ovarian)$coef)
  res3 <- optim(start, lognorm.inv.gauss.llik, Y = ovarian$futime,
                X=cbind(rep(1, nrow(ovarian)), ovarian$ecog.ps),
                C=(ovarian$fustat==0)*1, method="BFGS")  
  print(res3)
} else {
  data <- read.table("ajps2002-full1b.txt", header=T)

  toy <- survreg(Surv(acttime, d) ~ hcomm, dist='lognormal', scale=0, data=data) 
  
  res <- survreg(Surv(acttime, d) ~ hcomm + scomm + hfloor + sfloor + demsnmaj
                 + demhsmaj + prespart + orderent + stafcder +
                 prevgenx + lethal + 
                 deathrt1 + hosp01 + I(hospdisc/10000) + hhosleng + femdiz01 +
                 mandiz01 + peddiz01 + acutediz + orphdum + natreg +
                 I(natregsq/10000) + wpnoavg3 + vandavg3 + condavg3,
                 dist='lognormal', scale=0, data=data) 
  res <- toy
  #print(summary(res))

  #start <- c(res$scale, res$coef)
  #res1 <- optim(start, lognorm.llik, Y = data$acttime,
  #              X=model.matrix(res), C=(data$d==0)*1, method="BFGS")
  #print(res1)
  
  start <- c(-0.21, 0.54, 2.92, -2.50)
  print(lognorm.inv.gauss.llik(start, Y=data$acttime,
                               X=model.matrix(res), C=(data$d==0)*1)) 
  res2 <- optim(start, lognorm.inv.gauss.llik, Y = data$acttime,
                X=model.matrix(res), C=(data$d==0)*1, method="BFGS")
  print(res2)
}
