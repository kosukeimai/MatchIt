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
# Model 1 (Democrats)
#==================================

#creating full dataset
dta.full <- subset(dta,select=c(pdcanid, dviswom , dvisman , demcan1 ,
                     dempty , rideo, dproj , demft , aware))
dta.full <- na.omit(dta.full)
#now creating matched data
m1 <- matchit(dviswom ~ demcan1 + dempty + rideo +
           dproj + demft + aware, data=dta.full)
dta.match <- match.data(m1)

#creating possible RHS specifications
xvar <- names(dta.full)[2:length(names(dta.full))]
xvar <- xvar[xvar!="dvisman"]

#another way to show the curse of dimensionality
foo<-apply(dta.full[,xvar],2,table)
kk <- 1
for(i in 1:length(foo)){
  kk <- kk*length(foo[[i]])
}
# really hard to summarize 380108484 possible covariate cells!!

k <- length(xvar)
for(i in 1:k){
  for(j in i:k){
    xvar <- c(xvar,paste("I(",xvar[i],"*",xvar[j],")",sep=""))
  }
}

totk <- length(xvar)
tot <- 0
for(i in 1:totk){
  tot <- tot+choose(totk,i)
}
# 3.436e+10 total specifications

sims <- 1000
ate.full <- ate.match <- NULL
simseq <- floor(seq(sims/10,sims,sims/10))
for(i in 1:sims){
  sp <- FALSE
  while(!sp){ #only sampling specns with ttt on RHS
    xs <- xvar[sample(1:totk,size=sample(1:totk,size=1),replace=F)]
    if(length(grep("dviswom",xs))!=0){sp <- TRUE}
  }
  fml <- as.formula(paste("pdcanid ~", paste(xs, collapse= " + ")))
  ff <- zsim(fml,dta.full,num=10,zz=F)
  ate.full <- c(ate.full,ff)
  mm <- zsim(fml,dta.match,num=10,zz=F)
  ate.match <- c(ate.match,mm)
  if(i%in%simseq){
    cat(10*which(simseq==i),"%...",sep="")
    if(i==sims){cat("Done\n")}
  }
}
doverlay(ate.full,ate.match)



#==================================
# Model 3 (Republicans)
#==================================

dta.full <- subset(dta,select=c(prcanid, rviswom, rvisman, repcan1 , goppty ,
              rideo , rproj , repft , aware))
dta.full <- na.omit(dta.full)
m1 <- matchit(rviswom ~ rvisman + repcan1 + goppty + rideo + rproj +
              repft + aware, data=dta.full, caliper=0.1)
dta.match <- match.data(m1)

#creating possible RHS specifications
xvar <- names(dta.full)[2:length(names(dta.full))]
#xvar <- xvar[xvar!="rvisman"]

k <- length(xvar)
for(i in 1:k){
  for(j in i:k){
    xvar <- c(xvar,paste("I(",xvar[i],"*",xvar[j],")",sep=""))
  }
}
totk <- length(xvar)
tot <- 0
for(i in 1:totk){
  tot <- tot+choose(totk,i)
}

sims <- 1000
ate.full <- ate.match <- NULL
simseq <- floor(seq(sims/10,sims,sims/10))
for(i in 1:sims){
  sp <- FALSE
  while(!sp){ #only sampling specns with ttt on RHS
    xs <- xvar[sample(1:totk,size=sample(1:totk,size=1),replace=F)]
    if(length(grep("rviswom",xs))!=0){sp <- TRUE}
  }
  fml <- as.formula(paste("prcanid ~", paste(xs, collapse= " + ")))
  ff <- zsim(fml,dta.full,num=1,zz=F,treat="rviswom")
  ate.full <- c(ate.full,ff)
  mm <- zsim(fml,dta.match,num=1,zz=F,treat="rviswom")
  ate.match <- c(ate.match,mm)
  if(i%in%simseq){
    cat(10*which(simseq==i),"%...",sep="")
    if(i==sims){cat("Done\n")}
  }
}
doverlay(ate.full,ate.match)
