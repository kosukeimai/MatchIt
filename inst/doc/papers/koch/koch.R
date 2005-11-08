#Koch replication
rm(list=ls())
library(MatchIt)
library(foreign)
library(Zelig)
library(mvtnorm)
library(combinat)
library(lattice)
setwd("c:/R/match/docs/papers/koch/")
#setwd("c:/match/docs/papers/koch/")
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
dta.match <- match.data(m1)
summary(m1)

# QQ Plot
trellis.device(device="pdf",file="kochQQ.pdf",color=FALSE,width=4,height=4)
par(mar=c(2, 2, 2, 2) + 0.1, cex.lab=0.6, cex.axis=0.6,
    mgp=c(1,0.5,0), cex.main=0.8, cex=0.8, bg="white")
eqqplot(m1$distance[m1$treat==1],m1$distance[m1$treat==0],
        xlim=range(m1$distance),ylim=range(m1$distance),
        main="QQ Plot of Propensity Score",
        xlab="Treated Units", ylab="Control Units",
        pch=16, cex=0.5)
abline(a=0,b=1)
eqqplot(m1$distance[m1$treat==1 & m1$weights==1],
        m1$distance[m1$treat==0 & m1$weights==1], addit=T)
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
  tc <- tc*length(cc[[i]])
}

# Sensitivity
start <- paste(yy,"~",tt)
coef <- mcoef <- NULL
N <- 1
total <- 1
for (i in N:(length(xvars)-1))
  total <- total + ncol(combn(xvars, i)) 
cat("\n I'm going to run", total, "regressions!\n")

coef <- mcoef <- rep(0,total)
counter <- 1
for (i in N:(length(xvars)-1)) {
  allsubset <- combn(xvars, i)
  for (j in 1:ncol(allsubset)) {
    ftmp <- start
    for (k in 1:i)
      ftmp <- paste(ftmp, "+", allsubset[k,j])
    ftmp <- as.formula(ftmp)
    tmp0 <- lm(ftmp,  data = dta.full)
    coef[counter] <- tmp0$coefficient[2]
    tmp1 <- lm(ftmp, data = dta.match)
    mcoef[counter] <- tmp1$coefficient[2]
    counter <- counter + 1
  }
  cat(i,"covariates:",counter,"\n")
  cat(date(), "\n\n")
}
save.image("koch-11-7-05.Rdata")

# Plotting Densities
#setwd("c:/match/docs/papers/koch/writeup")
#setwd("c:/R/match/docs/papers/koch/writeup")
trellis.device(device="pdf",file="kochdens.pdf",color=FALSE,width=6,height=3)
par(mar=c(2, 2, 2, 2) + 0.1, cex.lab=0.6, cex.axis=0.6,
    mgp=c(1,0.5,0), cex.main=0.5, cex=0.8, bg="white")
doverlay(mcoef,coef,lwd=2, bw=0.004,
         xlab="Estimated average treatment effect", leg=F)
arrows(0.007,22, coefficients(res)[2],0, length=0.1)
text(0.07,9,"Raw data")
text(-0.055,50,"Matched\ndata")
text(0.007,27,"Point estimate \n of raw data")
dev.off()

