# File to create sample graphics and output for documentation

library(MatchIt)
data(lalonde)
m.out <- matchit(treat ~ re74+re75+educ+black+hispan+age, data=lalonde, method="nearest")

print(summary(m.out))
#postscript("docs/figs/qqplotnn.ps")
#par(mfrow=c(1,2))
#plot(m.out)
#dev.off()

#postscript("docs/figs/jitterplotnn.ps")
#plot(m.out, type="jitter")
#dev.off()

m.outs <- matchit(treat ~ re74+re75+educ+black+hispan+age, data=lalonde, method="subclass")
print(summary(m.outs))
#postscript("docs/figs/qqplotsub.ps")
#par(mfrow=c(1,2))
#plot(m.outs)
#dev.off()

#postscript("docs/figs/jitterplotsubs.ps")
#plot(m.outs, type="jitter")
#dev.off()



