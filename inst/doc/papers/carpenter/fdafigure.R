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
           legend = c("Full Data", "Matched Data"), ...)
  }
}

trellis.device(device="pdf",file="fdadens.pdf",color=FALSE,width=6,height=4)
par(mar=c(2, 2, 2, 2) + 0.1, cex.lab=0.6, cex.axis=0.6,
    mgp=c(1,0.5,0), cex.main=0.5, cex=0.8)
doverlay(mate,ate,lwd=2,
         xlab="Estimated Average Treatment Effect", leg=F)
arrows(ate[length(ate)], 0.12, ate[length(ate)],0, length=0.1)
text(-20,0.04,"Full Data")
text(-40,0.08,"Matched\nData")
text(-28,0.1275,"Point Estimate \n of Full Data")
dev.off()
