rm(list=ls())
library(lattice)
load("matchfdaLin.RData")

##
## calculation of qoi
##
library(MASS)
library(survival)

## matched data
mcall <- res$call
mcall$data <- mdata
mres <- eval(mcall)
## original data
mres <- res

coef <- coef(mres)
var <- vcov(mres)
nsims <- 100000
par <- mvrnorm(nsims, mu=c(coef(mres), log(mres$scale)), Sigma=vcov(mres))
scale2 <- exp(par[,ncol(par)])^2
par <- par[,-ncol(par)]

x <- model.matrix(mres)
y <- mres$y
x0 <- x1 <- x[x[,"treat"] == 1,] # treated only
x0[,"treat"] <- 0  # counterfactual for the treated
y <- y[x[,"treat"] == 1,] # observed Y(1) for the treated

att <- apply(ifelse(y[,2]>0.5, matrix(rep(exp(y[,1]), nsims), nrow=nsims, byrow=TRUE),
                    exp(par%*%t(x1)+0.5*scale2)) -
                    exp(par%*%t(x0)+0.5*scale2), 1, mean)
cat("ATT: ", mean(att), " (", sd(att), ")", sep="")

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
par(mar=c(2.5, 2, 2, 2) + 0.1, cex.lab=0.8, cex.axis=0.8,
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
