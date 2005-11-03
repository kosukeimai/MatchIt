rm(list=ls())
library(lattice)
load("matchfdaLin.RData")

##
## calculation of qoi
##
library(MASS)
library(survival)
library(Zelig)

## matched data
mcall <- res$call
mcall$data <- mdata
mres <- eval(mcall)
## original data
#mres <- res


x <- model.matrix(mres)
y <- mres$y
x0 <- x1 <- x[x[,"treat"] == 1,] # treated only
x0[,"treat"] <- 0  # counterfactual for the treated
y <- y[x[,"treat"] == 1,] # observed Y(1) for the treated

beta <- coef(mres)
sigma2 <- mres$scale^2

## delta method
est <- mean(ifelse(y[,2]>0.5, exp(y[,1]), exp(x1%*%beta+0.5*sigma2)) -
                    exp(x0%*%beta+0.5*sigma2))

grad <- matrix(0, nrow=length(beta)+1, ncol=1)
for (i in 1:nrow(x1)) {
  grad <- grad + ifelse(y[i,2]>0.5, matrix(0, nrow=length(beta)+1, ncol=1),
                        matrix(c(exp(x1[i,]%*%beta+0.5*sigma2)*x1[i,],
                                 exp(x1[i,]%*%beta+0.5*sigma2)*sigma2), ncol=1)
                        -matrix(c(exp(x0[i,]%*%beta+0.5*sigma2)*x0[i,],
                                  exp(x0[i,]%*%beta+0.5*sigma2)*sigma2), ncol=1))
}
grad <- grad/nrow(x1)

cat("ATT: ", est, " (", sqrt(t(grad)%*%vcov(mres)%*%grad), ")", sep="")

##
## density figure
##
doverlay <- function(x1, x0, xlab = "", main = "", lines = FALSE,
                     leg=T, ...)
{
  x <- c(x1,x0)
  dx1 <- density(x1)#, from = minobs, to = maxobs)
  dx0 <- density(x0)#, from = minobs, to = maxobs)
  plot(dx1, type = "l", ylab = "Density", 
       xlab = xlab, main = main, xlim=range(c(dx1$x,dx0$x)),
       ylim = range(c(dx1$y,dx0$y)),...)
  lines(dx0, lty=2, ...)
  if(leg){
    legend(minobs, max(c(dx1$y, dx0$y)), lty = 1:2, 
           legend = c("Raw data", "Matched data"), ...)
  }
}

pdf("fdadens.pdf", paper="special", height=3.5, width=6)
par(mar=c(2.5, 2.5, 2, 2) + 0.1, cex.lab=0.8, cex.axis=0.8,
    mgp=c(1.5,0.5,0), cex.main=0.5, cex=0.8, bg="white")
doverlay(mate,ate,lwd=2,
         xlab="Estimated in-sample average treatment effect", leg=F)
arrows(ate[length(ate)], 0.11, ate[length(ate)],0, length=0.1)
text(-75,0.05,"Raw data")
text(-42,0.05,"Matched\ndata")
text(-52.5,0.14,"Point estimate of \n Carpenter's specification \n using raw data")
dev.off()

##
## eCDF figure
##
sm <- summary(m.out, standardize=TRUE)
pdf("eCDF.pdf")
plot(x=sm$sum.all[1:19,"eCDF Max"], y=sm$sum.matched[1:19,"eCDF Max"],
     xlab="before matching", ylab="after matching", pch=19)
points(x=sm$sum.all[1:19,"eCDF Mean"], y=sm$sum.matched[1:19,"eCDF Mean"], pch=22)
points(x=sm$sum.all[1:19,"eCDF Med"], y=sm$sum.matched[1:19,"eCDF Med"], pch=24)
abline(0,1)
dev.off()
