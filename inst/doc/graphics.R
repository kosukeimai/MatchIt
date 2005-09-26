# File to create sample graphics and output for documentation

library(MatchIt)
library(lattice)
data(lalonde)
m.out <- matchit(treat ~ re74+re75+educ+black+hispan+age, data=lalonde, method="nearest")

print(summary(m.out))
#trellis.device(device="pdf",file="qqplotnn.pdf",color=FALSE,width=8,height=4)
postscript(file="figs/qqplotnn1.eps")
par(mar=c(2, 2, 2, 2) + 0.1, cex.lab=0.7, cex.axis=0.5,
    mgp=c(1,0.5,0), cex.main=0.8, cex=1, bg="white", mfrow=c(1,2))
plot(m.out, which.xs=c("re74", "re75", "educ"), interactive=FALSE)
dev.off()

postscript(file="figs/qqplotnn2.eps")
par(mar=c(2, 2, 2, 2) + 0.1, cex.lab=0.7, cex.axis=0.5,
    mgp=c(1,0.5,0), cex.main=0.8, cex=1, bg="white", mfrow=c(1,2))
plot(m.out, which.xs=c("black", "hispan", "educ"), interactive=FALSE)
dev.off()

#trellis.device(device="pdf",file="jitterplotnn.pdf",color=FALSE,width=8,height=4)
postscript(file="figs/jitterplotnn.eps")
par(mar=c(2, 2, 2, 2) + 0.1, cex.lab=0.7, cex.axis=0.5,
    mgp=c(1,0.5,0), cex.main=0.8, cex=1, bg="white")
plot(m.out, type="jitter", interactive=FALSE)
dev.off()

m.outs <- matchit(treat ~ re74+re75+educ+black+hispan+age, data=lalonde, method="subclass")
print(summary(m.outs))
#trellis.device(device="pdf",file="qqplotsub.pdf",color=FALSE,width=8,height=4)
pdf(file="figs/qqplotsub.pdf")
par(mar=c(2, 2, 2, 2) + 0.1, cex.lab=0.7, cex.axis=0.5,
    mgp=c(1,0.5,0), cex.main=0.8, cex=1, bg="white", mfrow=c(1,2))
plot(m.outs)
dev.off()

#trellis.device(device="pdf",file="jitterplotsub.pdf",color=FALSE,width=8,height=4)
pdf(file="figs/jitterplotsub.pdf")
par(mar=c(2, 2, 2, 2) + 0.1, cex.lab=0.7, cex.axis=0.5,
    mgp=c(1,0.5,0), cex.main=0.8, cex=1, bg="white")
plot(m.outs, type="jitter")
dev.off()



