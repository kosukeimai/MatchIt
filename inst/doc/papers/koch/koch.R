#Koch replication
rm(list=ls())
library(Matchit)
library(foreign)
library(Zelig)
library(mvtnorm)
library(combinat)
#setwd("c:/R/match/docs/papers/koch/")
setwd("c:/match/docs/papers/koch/")
source("fn.R")
dta <- read.dta("genharvard.dta")

#==================================
# Model 3 (Republicans)
#==================================

#creating full dataset 
dta.full <- subset(dta,select=c(prcanid, rviswom, rvisman, repcan1 , goppty ,
              rideo , rproj , repft , aware), subset=c(rvisman==0 & voter==1))
dta.full <- na.omit(dta.full)

#matching and creating matched dataset
#now testing sensitivity
#fml <- as.formula(prcanid ~ I(rviswom==1 | rvisman==1) + repcan1 + goppty + rideo +rproj + repft + aware)
#fml <- as.formula(prcanid ~ rviswom + rvisman + repcan1 + goppty + rideo +rproj + repft + aware)
fml <- as.formula(prcanid ~ rviswom + repcan1 + goppty + rideo + rproj + repft + aware)
res <- zelig(fml, data=dta.full, model="ls")
xvars <- names(res$coefficients)
xvars <- xvars[3:length(xvars)]
tt <- attr(terms(res),"term.labels")[1]
yy <- attr(terms(res),"variables")[[2]]
mfml <- as.formula(paste(tt,"~",paste(xvars,collapse=" + ")))
#m1 <- matchit(mfml, caliper=0.1,data=dta.full)
m1 <- matchit(mfml,data=dta.full)
dta.match <- match.data(m1)
mres <- zelig(fml, data=dta.match, model="ls")

#original models
summary(ols1 <- lm(fml,data=dta.full))
summary(ols2 <- lm(fml,data=dta.match))
#sampling from posterior
f1 <- zsim(fml,dta.full,num=1000,zz=F,treat="rviswom")
f2 <- zsim(fml,dta.match,num=1000,zz=F,treat="rviswom")
#overlaying densities
doverlay(f1,f2)

#sensitivity
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
    tmp <- lm(ftmp,  data = dta.full) 
    coef[counter] <- tmp$coefficient[tt]
    tmp <- lm(ftmp, data = dta.match) 
    mcoef[counter] <- tmp$coefficient[tt]
    counter <- counter + 1
  }
  cat(i,"covariates:",length(coef),"\n")
  cat(date(), "\n\n")
}
coef[length(coef)] <- res$coefficients[tt]
mcoef[length(mcoef)] <- mres$coefficients[tt] 

#plotting
doverlay(coef,mcoef)
abline(v=coefficients(res)[tt])

#visualizing extrapolation
lm.full <- lm(prcanid~repcan1+rviswom, data=dta.full)
lm.match <- lm(prcanid~repcan1+rviswom, data=dta.match)
par(mfrow=c(2,1))
plot(dta.full$repcan1, jitter(dta.full$prcanid),
     col=(dta.full$rviswom+2*(1-dta.full$rviswom)), cex=0.5,
     pch=16, main="Full Data", ylab="Outcome (jittered)", xlab="Candidate Ideology")
abline(coef=c(sum(coefficients(lm.full)[c(1,3)]),coefficients(lm.full)[2]))
abline(coef=c(coefficients(lm.full)[1],coefficients(lm.full)[2]))
plot(dta.match$repcan1,jitter(dta.match$prcanid),
     col=(dta.match$rviswom+2*(1-dta.match$rviswom)), cex=0.5,
     pch=16, main="Matched Data", ylab="Outcome (jittered)", xlab="Candidate Ideology")
abline(coef=c(sum(coefficients(lm.match)[c(1,3)]),coefficients(lm.match)[2]))
abline(coef=c(coefficients(lm.match)[1],coefficients(lm.match)[2]))

tot <- 0
for(i in 1:6){
  tot <- tot+choose(6,i)
}
