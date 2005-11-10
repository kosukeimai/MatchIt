#Koch replication
rm(list=ls())
library(MatchIt)
library(foreign)
library(Zelig)
library(mvtnorm)
library(combinat)
library(lattice)
#setwd("c:/R/match/docs/papers/koch/")
#setwd("c:/match/docs/papers/koch/")
setwd("matchit/docs/papers/koch/")
source("fn.R")
dta <- read.dta("genharvard.dta")

#======================================
# Model 3 (Republican Male Candidates)
#======================================

#creating full dataset 
dta.full <- subset(dta,select=c(prcanid, rviswom, rvisman, repcan1 , goppty ,
              rideo , rproj , repft , aware, repwom), subset=c(repman==1 & voter==1))
dta.full <- na.omit(dta.full)

# Matching and creating matched dataset
fml <- as.formula(prcanid ~ I(rvisman==0) + repcan1 + goppty + rideo +
                  rproj + repft + aware)
res <- zelig(fml, data=dta.full, model="ls")
xvars <- names(res$coefficients)
xvars <- xvars[3:length(xvars)]
tt <- attr(terms(res),"term.labels")[1]
yy <- attr(terms(res),"variables")[[2]]
mfml <- as.formula(paste(tt,"~",paste(xvars,collapse=" + ")))
m1 <- matchit(mfml, method="optimal", data=dta.full)
#m1 <- matchit(mfml, method="optimal", data=dta.full, ratio=2)
#m1 <- matchit(mfml, method="full", data=dta.full)
#m1 <- matchit(mfml, data=dta.full, ratio=2,replace=T)
#m1 <- matchit(mfml, data=dta.full, exact=c("goppty","rideo","repft"))
#m1 <- matchit(mfml, data=dta.full, exact=c("goppty"))
#m1 <- matchit(mfml, data=dta.full, exact=c("goppty","rideo"))
#m1 <- matchit(mfml, data=dta.full, exact=c("goppty","repft"))
#m1 <- matchit(mfml, data=dta.full, exact=c("rideo","repft"))
#m1 <- matchit(mfml, data=dta.full, exact=c("repft"))
#m1 <- matchit(mfml, data=dta.full, exact=c("rideo"))
#m1 <- matchit(mfml, data=dta.full,
#              discard="hull.both")
#m1 <- matchit(mfml, data=dta.full,
#              discard="hull.treat")
#m1 <- matchit(mfml, data=dta.full,
#              discard="hull.control")
#m1 <- matchit(mfml, data=dta.full,
#              discard="both")
#m1 <- matchit(mfml, data=dta.full,
#              discard="treat")
#m1 <- matchit(mfml, data=dta.full,
#              discard="control")
#m1 <- matchit(mfml, data=dta.full, ratio=2, discard="both")
#m1 <- matchit(mfml, data=dta.full, ratio=2, discard="both",replace=T)
#m1 <- matchit(mfml, data=dta.full, ratio=2, discard="control", replace=T)
#m1 <- matchit(mfml, data=dta.full, ratio=2, discard="both",
#              replace=T, exact="goppty")
#m1 <- matchit(mfml, data=dta.full, ratio=3, discard="both",
#              replace=T, exact=c("goppty"))
#m1 <- matchit(mfml, data=dta.full, ratio=4, discard="both",
#              replace=T, exact=c("goppty","repft"))
#m1 <- matchit(mfml, data=dta.full, ratio=4, discard="both",
#              replace=T, exact=c("repft"))
#m1 <- matchit(mfml, data=dta.full, ratio=4, discard="both",
#              replace=T, exact=c("rideo"))
#m1 <- matchit(mfml, data=dta.full, ratio=4, discard="both",
#              replace=T, exact=c("goppty"))
dta.match <- match.data(m1)
dta.full.controls <- subset(dta.full,subset=c(dta.full$rvisman!=0))
dta.full.treated <- subset(dta.full,subset=c(dta.full$rvisman==0))
dta.match.controls <- dta.match[dta.match$rvisman!=0,]
dta.match.treated <- dta.match[dta.match$rvisman==0,]
summary(m1)

# QQ Plot
pdf("kochqq.pdf", paper="special", height=4, width=4)
par(mar=c(2.5, 2.5, 2, 2) + 0.1, cex.lab=0.8, cex.axis=0.8,
    mgp=c(1.5,0.5,0), cex.main=0.8, cex=0.8, bg="white")
eqqplot(m1$distance[m1$treat==1],m1$distance[m1$treat==0],
        xlim=range(m1$distance),ylim=range(m1$distance),
        main="QQ Plot of Propensity Score",
        xlab="Treated Units", ylab="Control Units",
        pch=16, cex=0.5)
abline(a=0,b=1)
eqqplot(m1$distance[m1$treat==1 & m1$weights==1],
        m1$distance[m1$treat==0 & m1$weights==1], addit=T)
#eqqplot(m1$distance[m1$treat==1 & m1$weights==1],
#        sample(m1$distance[m1$treat==0], 5000, replace=TRUE, prob=m1$weights[m1$treat==0]),        
#        addit=T)
text(0.27,0.42,"Matched Data",cex=0.8)
text(0.45,0.24,"Raw Data",cex=0.8)
dev.off()

# Outcome Regression
summary(lm(fml, data=dta.full))
summary(lm(fml, data=dta.match))

# Total number of cells
cc <- apply(dta.full[,xvars],2,table)
tc <- 2 #2 treatment values
for(i in 1:length(cc)){
#  print(length(cc[[i]]))
  tc <- tc*length(cc[[i]])
}

# Sensitivity
# unconditional effects =========
#start <- paste(yy,"~",tt)
# conditional effects ==========
start <- paste(yy,"~") #conditional effects
coef <- mcoef <- NULL
N <- 1
total <- 1
for (i in N:(length(xvars)-1))
  total <- total + ncol(combn(xvars, i)) 
cat("\n I'm going to run", total, "regressions!\n")

cil <- mcil <- coef <- mcoef <- rep(0,total)
counter <- 1
for (i in N:(length(xvars)-1)) {
  allsubset <- combn(xvars, i)
  for (j in 1:ncol(allsubset)) {
    ftmp <- start
    for (k in 1:i)
      ftmp <- paste(ftmp, "+", allsubset[k,j])
    ftmp <- as.formula(ftmp)
    # unconditional effects =============
#    tmp0 <- lm(ftmp,  data = dta.full)
#    coef[counter] <- tmp0$coefficient[2]
#    cil[counter] <- 2*1.96*sqrt(vcov(tmp0)[2,2])
#    tmp1 <- lm(ftmp, data = dta.match, weights=m1$weights[m1$weights!=0])
#    tmp1 <- lm(ftmp, data = dta.match)
#    mcoef[counter] <- tmp1$coefficient[2]
#    mcil[counter] <- 2*1.96*sqrt(vcov(tmp1)[2,2])
    #conditional effects ============
    res0 <- zelig(ftmp,
                  data=dta.full.controls, model="ls")
    xout0 <- setx(res0,data=dta.full.treated,cond=T)
    out0 <- sim(res0,x=xout0)
    quant0 <-  quantile(out0$qi$ate.ev,c(0.025,0.95))
    cil[counter] <- quant0[2]-quant0[1]
    coef[counter] <- mean(out0$qi$ate.ev)
    res1 <- zelig(ftmp,
                  data=dta.match.controls, model="ls")
    xout1 <- setx(res1,data=dta.match.treated,cond=T)
    out1 <- sim(res1,x=xout1)
    quant1 <-  quantile(out1$qi$ate.ev,c(0.025,0.95))
    mcil[counter] <- quant1[2]-quant1[1]
    mcoef[counter] <- mean(out1$qi$ate.ev)
    counter <- counter + 1
  }
  cat(i,"covariates:",counter,"\n")
  cat(date(), "\n\n")
  # unconditional effects ==============
#  coef[length(coef)] <- res$coefficients[2]
#  cil[counter] <- 2*1.96*sqrt(vcov(res)[2,2])
#  res1 <- lm(fml, data=dta.match)
#  mcoef[length(mcoef)] <- res1$coefficients[2]
#  mcil[counter] <- 2*1.96*sqrt(vcov(res1)[2,2])
  # conditional effects =============
  res0 <- zelig(prcanid ~ repcan1 + goppty + rideo +
                  rproj + repft + aware,
                data=dta.full.controls, model="ls")
  xout0 <- setx(res0,data=dta.full.treated,cond=T)
  out0 <- sim(res0,x=xout0)
  quant0 <-  quantile(out0$qi$ate.ev,c(0.025,0.95))
  cil[counter] <- quant0[2]-quant0[1]
  coef[counter] <- mean(out0$qi$ate.ev)
  res1 <- zelig(prcanid ~ repcan1 + goppty + rideo +
                  rproj + repft + aware,
                data=dta.match.controls, model="ls")
  xout1 <- setx(res1,data=dta.match.treated,cond=T)
  out1 <- sim(res1,x=xout1)
  quant1 <-  quantile(out1$qi$ate.ev,c(0.025,0.95))
  mcil[counter] <- quant1[2]-quant1[1]
  mcoef[counter] <- mean(out1$qi$ate.ev)  
}
#CI lengths
mean(cil)
mean(mcil)

# Plotting Densities
#setwd("c:/match/docs/papers/koch/writeup")
#setwd("c:/R/match/docs/papers/koch/writeup")

pdf("kochdens.pdf", paper="special", height=3.5, width=6)
par(mar=c(2.5, 2.5, 2, 2) + 0.1, cex.lab=0.8, cex.axis=0.8,
    mgp=c(1.5,0.5,0), cex.main=0.5, cex=0.8, bg="white")
doverlay(mcoef,coef,lwd=2, bw=0.004,
         xlab="Estimated average treatment effect", leg=F)
arrows(0.007,22, coefficients(res)[2],0, length=0.1)
text(0.07,9,"Raw data")
text(-0.055,50,"Matched\ndata")
text(0.007,27,"Point estimate \n of raw data")
dev.off()

#saving image
save.image("koch-11-7-05.Rdata")
