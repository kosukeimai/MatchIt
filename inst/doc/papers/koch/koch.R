#Koch replication
rm(list=ls())
library(Matchit)
library(foreign)
library(Zelig)
library(mvtnorm)
setwd("c:/R/match/docs/papers/koch/")
dta <- read.dta("genharvard.dta")
draw.ols <- function(ols1,sims=1000){
  rinvchisq <- function (n, df, scale = 1/df){
    if ((length(scale) != 1) & (length(scale) != n))
      stop("scale should be a scalar or a vector of the same length as x")
    if (df <= 0)
      stop("df must be greater than zero")
    if (any(scale <= 0))
      stop("scale must be greater than zero")
    return((df * scale)/rchisq(n, df = df)) }
  s2 <- summary(ols1)$sigma^2 
  bhat <- ols1$coefficients 
  k <- length(bhat)
  n <- nrow(model.matrix(ols1))
  vhat <- summary(ols1)$cov.unscaled 
  sim.sig2 <- rinvchisq(sims,n-k,s2)
  sim.beta <- matrix(0,sims,k)
  for(i in 1:sims){
    sim.beta[i,] <- rmvnorm(1,mean=bhat,sigma=vhat*sim.sig2[i])
  }
  sim.beta <- as.data.frame(cbind(sim.beta,sim.sig2))
  names(sim.beta) <- c(names(bhat),"sig2")
  return(sim.beta)
}
#overlaying densities
doverlay <- function(x1, x0, xlab = "", main = "", lines = FALSE, ...)
{
  x <- c(x1,x0)
  minobs <- min(x)
  maxobs <- max(x)
  dx1 <- density(x1, from = minobs, to = maxobs)
  dx0 <- density(x0, from = minobs, to = maxobs)
    matplot(dx0$x, cbind(dx1$y, dx0$y), type = "l", ylab = "Density", 
            xlab = xlab, main = main, ...)
  legend(0.7, max(c(dx1$y, dx0$y)), lty = 1:2, col = 1:2, 
         legend = c("Treatment", "Control"), cex=0.5)
}

#ols1 <- lm(pdcanid ~ dviswom + dvisman + demcan1 + dempty + rideo +
#           dproj + demft + aware, data=dta)
# Table 2, Model 1
ols1 <- lm(pdcanid ~ dviswom + demcan1 + dempty + rideo +
           dproj + demft + aware, data=dta)
dta2 <- subset(dta,select=c(pdcanid, dviswom , dvisman , demcan1 ,
                     dempty , rideo, dproj , demft , aware))
dta2 <- na.omit(dta2)
#now matching (but here we're using the full dataset)
m1 <- matchit(dviswom ~ demcan1 + dempty + rideo +
           dproj + demft + aware, data=dta2, nearest=F)
neyman(pdcanid,m1)  #neyman estimate
#should treatment indicator be included on RHS??
#z.out <- zelig(pdcanid ~ dviswom + demcan1 + dempty + rideo +
#           dproj + demft + aware, data = match.data(m1),
#               model="ls")
#model adjustment
ols.m <- lm(pdcanid ~  demcan1 + dempty + rideo +
           dproj + demft + aware, data=match.data(m1,"control"))
b.sim <- as.matrix(draw.ols(ols.m))
b.sim <- b.sim[,1:(ncol(b.sim)-1)]
XX <- as.matrix(cbind(1,match.data(m1,"treat")[,m1$covariates]))
y.hat <- XX%*%t(b.sim)
ate <- matrix(0,1000,1)
for(i in 1:1000){
  ate[i] <- mean(dta2$pdcanid[dta2$dviswom==1]-y.hat[,i])
}
z.out <- zelig(pdcanid ~ dviswom + demcan1 + dempty + rideo +
           dproj + demft + aware, data = match.data(m1),
               model="ls")
x.out <- setx(z.out, data = match.data(m1),cond=T)
x.out1 <- x.out0 <- x.out
x.out1[,dimnames(x.out)[[2]]=="dviswom"] <- 1
x.out0[,dimnames(x.out)[[2]]=="dviswom"] <- 0
s.out <- sim(z.out, x = x.out0, x1 = x.out1)
summary(s.out)
