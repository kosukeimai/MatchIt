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

#creating full dataset
dta.full <- subset(dta,select=c(pdcanid, dviswom , dvisman , demcan1 ,
                     dempty , rideo, dproj , demft , aware))
dta.full <- na.omit(dta.full)
#now matching (but here we're using the full dataset)
m1 <- matchit(dviswom ~ demcan1 + dempty + rideo +
           dproj + demft + aware, data=dta.full)
dta.match <- match.data(m1)

fml <- as.formula(pdcanid ~ dviswom + demcan1 +
           dempty + rideo + dproj + demft + aware)
fml <- as.formula(pdcanid ~ dviswom + I(dviswom*rideo) + demcan1 +
           dempty + rideo + dproj + demft + aware)
foo <- zsim(fml,dta.full)
foo <- zsim(fml,dta.match)

