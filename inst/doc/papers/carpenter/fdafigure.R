setwd("~/research/matchit/docs/papers/carpenter")
library(lattice)
load("matchfda.RData")

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
           legend = c("Raw Data", "Matched Data"), ...)
  }
}

trellis.device(device="pdf",file="fdadens.pdf",color=FALSE,width=6,height=4)
par(mar=c(2, 2, 2, 2) + 0.1, cex.lab=0.6, cex.axis=0.6,
    mgp=c(1,0.5,0), cex.main=0.5, cex=0.8)
doverlay(mate,ate,lwd=2,
         xlab="Estimated Average Treatment Effect", leg=F)
arrows(ate[length(ate)], 0.12, ate[length(ate)],0, length=0.1)
text(-20,0.04,"Raw Data")
text(-40,0.08,"Matched\nData")
text(-28,0.1275,"Point Estimate \n of Raw Data")
dev.off()

data <- read.table("ajps2002-full1b.txt", header=T)
data$treat <- data$demsnmaj

library(Matchit)
mout <- matchit(treat ~ orderent + stafcder + prevgenx + lethal +
                deathrt1 + hosp01 + I(hospdisc/1000000) +
                hhosleng + femdiz01 + mandiz01 +
                peddiz01 + acutediz + orphdum + natreg +
                I(natregsq/1000) + wpnoavg3 + sqrt(wpnoavg3) +
                vandavg3 + condavg3, data=data, discard=1)
smout <- summary(mout)
pscore <- mout$data$pscore[mout$data$treat==0]

trellis.device(device="pdf",file="fdabal.pdf",color=FALSE,width=6,height=3)
par(mar=c(2, 2, 2, 2) + 0.1, cex.lab=0.6, cex.axis=0.6,
    mgp=c(1,0.5,0), cex.main=0.5, cex=0.8, mfrow=c(1,2))

plot(abs(smout$sum.all[-18,"T-stat"]),
      abs(smout$sum.matched[-18,"T-stat"]), pch=".", 
      xlim=c(0, max(abs(smout$sum.all[,"T-stat"]))),
      ylim=c(0, max(abs(smout$sum.all[,"T-stat"]))), 
      main="Covariate balance before and after matching",
      xlab="absolute t-statistics before matching",
      ylab="absolute t-statistics after matching")
polygon(c(-0.2,-0.2,1.96,1.96,-0.2), c(-0.2,1.96,1.96,-0.2,-0.2), density=-1,
        col="darkgray", border=NA)
points(abs(smout$sum.all[-18,"T-stat"]),
       abs(smout$sum.matched[-18,"T-stat"]), pch=19, cex=0.6)
abline(0,1)
text(5.5, 1.5, "propensity score", cex=0.6)
text(4.5, 0.1, "FDA review staff", cex=0.6)
text(4, 2, "Washington Post stories", cex=0.6)
text(4, 0.7, "Nightly TV News stories", cex=0.6)

plot(density(mout$data$pscore[mout$data$treat==1 & mout$psweights>0]),
     type="h", col="darkgray", xlim=c(0,1), ylim=c(0,4),
     main="Propensity score before and after matching",
     xlab="Estimated propensity score")
lines(density(pscore[mout$psweights[mout$data$treat==0]>0]), lwd=2)
lines(density(pscore), lwd=2, lty=2)
dev.off()
