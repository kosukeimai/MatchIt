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
                     dempty , rideo, dproj , demft , aware), subset=c(demwom==0))
dta.full <- subset(dta,select=c(pdcanid, dviswom , dvisman , demcan1 ,
                     dempty , rideo, dproj , demft , aware), subset=c(demwom==1))
#dta.full <- subset(dta,select=c(pdcanid, dviswom , dvisman , demcan1 ,
#                     dempty , rideo, dproj , demft , aware))
dta.full <- na.omit(dta.full)
#now creating matched data
m1 <- matchit(dviswom ~ demcan1 + dempty + rideo +
           dproj + demft + aware, data=dta.full)
m1 <- matchit(dviswom ~ demcan1 + I(demcan1^2) + I(demcan1^3) + dempty + rideo +
              dproj + demft + aware, data=dta.full,
              exact=c("dempty","rideo","demft","aware"))
m1 <- matchit(dviswom ~ demcan1 + I(demcan1^2) + I(demcan1^3) + dproj, data=dta.full,
              exact=c("dempty","rideo","demft","aware"))
m1 <- matchit(dviswom ~ demcan1 + dproj, data=dta.full,
              exact=c("dempty","rideo","demft","aware"))
m1 <- matchit(dviswom ~ demcan1 + dempty + rideo +
              dproj + demft + aware, data=dta.full,
              nearest=F, subclass=8, discard=1)
m1 <- matchit(as.numeric(dvisman==0) ~ demcan1 + dempty + rideo +
              dproj + demft + aware, data=dta.full,
              nearest=F, subclass=8, discard=1)
m1 <- matchit(as.numeric(dvisman==0) ~ demcan1 + dempty + rideo +
              dproj + demft + aware, data=dta.full,
              nearest=T)
m1 <- matchit(as.numeric(dviswom==0) ~ demcan1 + dempty + rideo +
           dproj + demft + aware, data=dta.full)
dta.match <- match.data(m1)
#summary(neyman(pdcanid,m1))

#model sensitivity
fml <- as.formula(pdcanid ~ dvisman + demcan1 + dempty + rideo + dproj + demft + aware)


#woman
fml <- as.formula(pdcanid ~ dviswom + demcan1 + dempty + rideo + dproj + demft + aware)

ccode <- function(tt){
  cc <- tt
  cc[cc==1] <- "red"
  cc[cc==0] <- "yellow"
  cc
}

ols1 <- lm(fml,data=dta.full)
plot(dta.full$demcan1,residuals(ols1),pch=16,col=ccode(dta.full$dviswom))
abline(coef(lm(residuals(ols1)~dta.full$demcan1,subset=c(dta.full$dviswom==1))))
abline(coef(lm(residuals(ols1)~dta.full$demcan1,subset=c(dta.full$dviswom==0))))

plot(dta.full$dempty,residuals(ols1),pch=16,col=ccode(dta.full$dviswom))
plot(dta.full$dproj,residuals(ols1),pch=16,col=ccode(dta.full$dviswom))
abline(coef(lm(residuals(ols1)~dta.full$dproj,subset=c(dta.full$dviswom==1))))
abline(coef(lm(residuals(ols1)~dta.full$dproj,subset=c(dta.full$dviswom==0))))
abline(h=0)

z.out <- zelig(fml, data = dta.full, model="ls")
x.out<- setx(z.out, data = dta.full,cond=T)
s.out <- sim(z.out, x = x.out, fn=NULL, num=1000)
summary(s.out)
hist(s.out$qi$ate.ev)

ols2 <- lm(fml,data=dta.match)
z.out <- zelig(fml, data = dta.match, model="ls")
x.out<- setx(z.out, data = dta.match,cond=T)
s.out <- sim(z.out, x = x.out, fn=NULL, num=1000)
summary(s.out)
hist(s.out$qi$ate.ev)

hist(s.out$qi$ate.ev)

f1 <- zsim(fml,dta.full,num=1000,zz=F,treat="dvisman")
f2 <- zsim(fml,dta.match,num=1000,zz=F,treat="dvisman")
doverlay(f1,f2)

f1 <- zsim(fml,dta.full,num=1000,zz=F,treat="dviswom")
f2 <- zsim(fml,dta.match,num=1000,zz=F,treat="dviswom")
ols1 <- lm(fml,data=dta.full)
ols2 <- lm(fml,data=dta.match)
doverlay(f1,f2)


#creating possible RHS specifications
xvar <- names(dta.full)[2:length(names(dta.full))]
xvar <- xvar[xvar!="dvisman"]
#xvar <- xvar[xvar!="dviswom"]

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
#    xs <- xvar[sample(1:totk,size=sample(1:totk,size=1),replace=F)]
#    xs <- xvar[c(1,sample(2:7,size=sample(1:6,size=1),replace=F),
#    sample(8:35,size=sample(1:28,size=1),replace=F))]
    xs <- xvar[c(1,sample(2:7,size=sample(1:6,size=1),replace=F),
                 sample(8:35,size=sample(1:28,size=1,
                               prob=rev(pnorm(seq(from=-4,to=0,length=28)))),
                        replace=F))]
#    dx <- xvar[sample(2:7,size=sample(1:6,size=1),replace=F)]
#    sx <- xvar[8:35]
#    for(i in 1:length(dx)){
#      sx <- grep(dx[i],sx,value=T)
#    }
#    sx <- xvar[xvar%in%dx]
#    sx <- xvar[sample(1:length(sx),size=sample(1:length(sx),size=1),replace=F)]
    if(length(grep("dviswom",xs))!=0){sp <- TRUE}
                                        #    if(length(grep("dvisman",xs))!=0){sp <- TRUE}
  }
  fml <- as.formula(paste("pdcanid ~", paste(xs, collapse= " + ")))
  ff <- mean(zsim(fml,dta.full,num=1,zz=F))
  ate.full <- c(ate.full,ff)
  mm <- mean(zsim(fml,dta.match,num=1,zz=F))
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
dta.full <- subset(dta,select=c(prcanid, rviswom, rvisman, repcan1 , goppty ,
              rideo , rproj , repft , aware),subset=c(rvisman==0))
dta.full <- na.omit(dta.full)
#m1 <- matchit(rviswom ~ rvisman + repcan1 + goppty + rideo + rproj +
#              repft + aware, data=dta.full, subclass=5, nearest=T)
#m1 <- matchit(rviswom ~ rvisman + repcan1 + goppty + rideo + rproj +
#              repft + aware, data=dta.full, subclass=5, nearest=T, sub.by="treat")
m1 <- matchit(rviswom ~ rvisman + repcan1 + goppty + rideo + rproj +
              repft + aware, data=dta.full)
m1 <- matchit(rviswom ~ rvisman + repcan1 + I(repcan1^2) + goppty + rideo + rproj +
              repft + aware + I((repcan1^2)*aware), data=dta.full)
m1 <- matchit(rviswom ~ repcan1 + goppty + rideo + rproj +
              repft + aware, data=dta.full, exact=T)
m1 <- matchit(rviswom ~ repcan1 + goppty + rideo + rproj +
              repft + aware, data=dta.full, exact=T)
dta.match <- match.data(m1)

#original models
fml <- as.formula(prcanid ~ rviswom + rvisman + repcan1 + goppty + rideo +rproj + repft + aware)
ols1 <- lm(fml,data=dta.full)
summary(ols2 <- lm(fml,data=dta.match))
#testing sensitivity
# (the following models provide estimates of ttt effects that do not meet conventional alpha levels)
fml <- as.formula(prcanid ~ rviswom + repcan1 + I(rviswom*repcan1))
fml <- as.formula(prcanid ~ rviswom + repcan1 + I(rviswom*repcan1) + rvisman)
fml <- as.formula(prcanid ~ rviswom + rvisman + repcan1 + I(rviswom*repcan1) + I(rvisman*repcan1))
fml <- as.formula(prcanid ~ rviswom + rvisman + repcan1 + I(rviswom*repcan1) + I(repcan1^2) + goppty)
fml <- as.formula(prcanid ~ rviswom + rvisman + repcan1 + I(rviswom*repcan1) + I(repcan1^2) + goppty + rideo)
fml <- as.formula(prcanid ~ rviswom + rvisman + repcan1 + I(rviswom*repcan1) + I(repcan1^2) + goppty + rideo +rproj)
fml <- as.formula(prcanid ~ rviswom + rvisman + repcan1 + I(rviswom*repcan1) + goppty + rideo +rproj)
fml <- as.formula(prcanid ~ rviswom + rvisman + repcan1 + I(rviswom*repcan1) + goppty + rideo +rproj + repft + aware)
f1 <- zsim(fml,dta.full,num=1000,zz=F,treat="rviswom")
f2 <- zsim(fml,dta.match,num=1000,zz=F,treat="rviswom")
doverlay(f1,f2)

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
