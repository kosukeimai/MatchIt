#Koch replication
rm(list=ls())
library(Matchit)
library(foreign)
library(Zelig)
setwd("c:/R/match/docs/papers/koch/")
dta <- read.dta("genharvard.dta")
#ols1 <- lm(pdcanid ~ dviswom + dvisman + demcan1 + dempty + rideo +
#           dproj + demft + aware, data=dta)
# Table 2, Model 1
ols1 <- lm(pdcanid ~ dviswom + demcan1 + dempty + rideo +
           dproj + demft + aware, data=dta)
dta2 <- subset(dta,select=c(pdcanid, dviswom , dvisman , demcan1 ,
                     dempty , rideo, dproj , demft , aware))
dta2 <- na.omit(dta2)
#now matching (but here we're using the full dataset)
m1 <- matchit(dviswom ~ demcan1 + dempty + rideo +
           dproj + demft + aware, data=dta2, nearest=F)
neyman(pdcanid,m1)  #neyman estimate
#model adjustment 
ols.m <- lm(pdcanid ~ dviswom + demcan1 + dempty + rideo +
           dproj + demft + aware, data=match.data(m1))
summary(ols.m)
#should treatment indicator be included on RHS??
#z.out <- zelig(pdcanid ~ dviswom + demcan1 + dempty + rideo +
#           dproj + demft + aware, data = match.data(m1),
#               model="ls")
z.out <- zelig(pdcanid ~ demcan1 + dempty + rideo +
           dproj + demft + aware, data = match.data(m1),
               model="ls")
x.out <- setx(z.out, data = match.data(m1), cond=T)
s.out <- sim(z.out, x = x.out)
summary(s.out)
