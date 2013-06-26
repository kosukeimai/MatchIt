# File to create sample graphics and output for documentation

library(MatchIt)
library(lattice)
data(lalonde)
m.out <- matchit(treat ~ re74+re75+educ+black+hispan+age, data=lalonde, method="nearest")

print(summary(m.out))

ps.options(family = c("Times"), pointsize = 8)

postscript(file="figs/qqplotnn1.eps", horizontal=FALSE, paper="special", width=2.75, height=2.75)
par(mar=c(2, 2, 2, 2) + 0.1, cex.lab=0.7, cex.axis=0.5,
    mgp=c(1,0.5,0), cex.main=0.8, cex=1, bg="white", mfrow=c(1,2))
plot(m.out, which.xs=c("re74", "re75", "educ"), interactive=FALSE)
dev.off()

pdf(file="figs/qqplotnn1.pdf", width=2.75, height=2.75, pointsize=8, family="Times")
par(mar=c(2, 2, 2, 2) + 0.1, cex.lab=0.7, cex.axis=0.5,
    mgp=c(1,0.5,0), cex.main=0.8, cex=1, bg="white", mfrow=c(1,2))
plot(m.out, which.xs=c("re74", "re75", "educ"), interactive=FALSE)
dev.off()

postscript(file="figs/qqplotnn2.eps", horizontal=FALSE, paper="special", width=2.75, height=2.75)
par(mar=c(2, 2, 2, 2) + 0.1, cex.lab=0.7, cex.axis=0.5,
    mgp=c(1,0.5,0), cex.main=0.8, cex=1, bg="white", mfrow=c(1,2))
plot(m.out, which.xs=c("black", "hispan", "educ"), interactive=FALSE)
dev.off()

pdf(file="figs/qqplotnn2.pdf", width=2.75, height=2.75, pointsize=8, family="Times")
par(mar=c(2, 2, 2, 2) + 0.1, cex.lab=0.7, cex.axis=0.5,
    mgp=c(1,0.5,0), cex.main=0.8, cex=1, bg="white", mfrow=c(1,2))
plot(m.out, which.xs=c("black", "hispan", "educ"), interactive=FALSE)
dev.off()

postscript(file="figs/jitterplotnn.eps", horizontal=FALSE, paper="special", width=5.5, height=3.5)
par(mar=c(2, 2, 2, 2) + 0.1, cex.lab=0.7, cex.axis=0.5,
    mgp=c(1,0.5,0), cex.main=0.8, cex=1, bg="white")
plot(m.out, type="jitter", interactive=FALSE)
dev.off()

pdf(file="figs/jitterplotnn.pdf", width=5.5, height=3.5, pointsize=8, family="Times")
par(mar=c(2, 2, 2, 2) + 0.1, cex.lab=0.7, cex.axis=0.5,
    mgp=c(1,0.5,0), cex.main=0.8, cex=1, bg="white")
plot(m.out, type="jitter", interactive=FALSE)
dev.off()


