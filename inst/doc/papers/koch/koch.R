#Koch replication
rm(list=ls())
library(Matchit)
library(foreign)
library(Zelig)
library(mvtnorm)
#setwd("c:/R/match/docs/papers/koch/")
setwd("c:/match/docs/papers/koch/")
source("fn.R")
dta <- read.dta("genharvard.dta")

#==================================
# Model 3 (Republicans)
#==================================

#creating full dataset 
dta.full <- subset(dta,select=c(prcanid, rviswom, rvisman, repcan1 , goppty ,
              rideo , rproj , repft , aware))
dta.full <- na.omit(dta.full)
#matching and creating matched dataset
m1 <- matchit(rviswom ~ rvisman + repcan1 + goppty + rideo + rproj +
              repft + aware, caliper=0.1,data=dta.full)
#m1 <- matchit(rviswom ~ rvisman + repcan1 + goppty + rideo + rproj +
#              repft + aware,data=dta.full)
dta.match <- match.data(m1)

#original models
fml <- as.formula(prcanid ~ rviswom + rvisman + repcan1 + goppty + rideo +rproj + repft + aware)
summary(ols1 <- lm(fml,data=dta.full))
summary(ols2 <- lm(fml,data=dta.match))
#sampling from posterior
f1 <- zsim(fml,dta.full,num=10000,zz=F,treat="rviswom")
f2 <- zsim(fml,dta.match,num=10000,zz=F,treat="rviswom")
#overlaying densities
doverlay(f1,f2)


