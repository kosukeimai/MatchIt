ft <- function(t, mu, sigma) {
        #return(-log(t)-log(sigma)-0.5*log(2*pi)-(log(t)-mu)^2/(2*sigma^2))
	return((1/(t*sigma*sqrt(2*pi)))*exp((-1/(2*sigma^2))*(log(t)-mu)^2))
}

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

  lnS <- pnorm((log(Y)-mu)/sigma, lower.tail = FALSE, log = TRUE)
  #lnh <- ft(Y, mu, sigma) - lnS
  lnh <- log(ft(Y, mu, sigma)) - lnS	
  lnhold <- dnorm(log(Y), mean=mu, sd=sigma, log = TRUE) - lnS

  #value <- sum(-sqrt(1-2*theta*lnS)/theta + (1-C)*lnh -
  #  0.5*(1-C)*log(1-2*theta*lnS)) 

  # New form: Include first term (1/theta)
  value <- sum((1/theta) - sqrt(1-2*theta*lnS)/theta + (1-C)*lnh -
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

  #toy <- survreg(Surv(acttime, d) ~ hcomm, dist='lognormal', scale=0, data=data) 
  
  res <- survreg(Surv(acttime, d) ~ hcomm + scomm + hfloor + sfloor + 
	         demhsmaj + demsnmaj + prespart + 
                 prevgenx + lethal + 
                 deathrt1 + acutediz + hosp01 + I(hospdisc/1000000) + hhosleng + 
                 femdiz01 + mandiz01 + peddiz01 + orphdum + natreg +
                 I(natregsq/10000) + vandavg3 +
                 wpnoavg3 + condavg3 + orderent + stafcder,
                 dist='lognormal', scale=0, data=data) 
  #res <- toy
  print(summary(res))
  #start <- c(res$scale, res$coef)
  #res1 <- optim(start, lognorm.llik, Y = data$acttime,
  #              X=model.matrix(res), C=(data$d==0)*1, method="BFGS")
  #print(res1)
  
  start <- c(-0.5, 0.5, res$coef)
  stata <- c(-0.47, 0.59, 2.90, -1.98, -0.93, 1.05, 0.64, 0.44, -0.91, -0.35, 
             0.001, -0.047, -0.426, -0.153, -0.145, 1.45, -0.01, -0.38, 0.04, -0.22, 
             -0.26, 0.015, -0.7, 0.0155, -0.0038, 0.0102, 0.0186, -0.0001)  
  res1 <- optim(start, lognorm.inv.gauss.llik, Y=data$acttime, X=model.matrix(res), 
                C=(data$d==0)*1, method="BFGS", control=list(maxit=1000000)) 
  print(res1)
  print(lognorm.inv.gauss.llik(stata, Y=data$acttime, X=model.matrix(res), C=(data$d==0)*1) )

  #start <- c(-0.21, 0.54, 2.92, -2.50)
  #res2 <- optim(start, lognorm.inv.gauss.llik, Y = data$acttime,
  #              X=model.matrix(toy), C=(data$d==0)*1, method="BFGS")
  #print(res2)
}
