#Koch replication
rm(list=ls())
library(Matchit)
library(foreign)
library(Zelig)
library(mvtnorm)
library(combinat)
library(lattice)
setwd("c:/R/match/docs/papers/koch/")
#setwd("c:/match/docs/papers/koch/")
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
fml <- as.formula(prcanid ~ rviswom + repcan1 + goppty + rideo + rproj + repft + aware)
res <- zelig(fml, data=dta.full, model="ls")
xvars <- names(res$coefficients)
xvars <- xvars[3:length(xvars)]
tt <- attr(terms(res),"term.labels")[1]
yy <- attr(terms(res),"variables")[[2]]
mfml <- as.formula(paste(tt,"~",paste(xvars,collapse=" + ")))
m1 <- matchit(mfml, data=dta.full, subclass=6, nearest=F, discard=1)
dta.match <- match.data(m1)

#total number of cells
cc <- apply(dta.full[,xvars],2,table)
tc <- 1
for(i in 1:length(cc)){
  tc <- tc*length(cc[[i]])
}

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
    tmp <- zelig(ftmp, data = dta.match, model="ls", by="psclass")
    wate <- 0
    for(l in 1:length(tmp)){
      wate <- wate + tmp[[l]]$coefficient[tt]*sum(dta.match$psclass==l & dta.match[,tt]==1)/sum(dta.match$psclass!=0 & dta.match[,tt]==1)
    }
    mcoef[counter] <- wate
    counter <- counter + 1
  }
  cat(i,"covariates:",length(coef),"\n")
  cat(date(), "\n\n")
}
coef[length(coef)] <- res$coefficients[tt]
mres <- zelig(fml, data=dta.match, model="ls",by="psclass")
wate <- 0
for(l in 1:length(tmp)){
  wate <- wate + tmp[[l]]$coefficient[tt]*sum(dta.match$psclass==l & dta.match[,tt]==1)/sum(dta.match$psclass!=0 & dta.match[,tt]==1)      
}
mcoef[length(mcoef)] <- wate

#tables
library(xtable)
sm1 <- summary(m1)
tab1 <- cbind(sm1$sum.all[1:7,c(1,2,4)],sm1$q.table[1:7,c(1,2,4) , 4])
row.names(tab1) <- c("Propensity Score", "Candidate Ideology",
                 "Perception of Party Ideology", "Respondent Ideology",
                 "Respondent Ideology * Feeling Thermometer",
                 "Feeling Thermometer", "Political Awareness")
xtable(tab1)


number.t <- sm1$q.table[nrow(sm1$q.table), , ][1, ]
number.c <- sm1$q.table[nrow(sm1$q.table), , ][2, ]
number.all <- sm1$q.table[nrow(sm1$q.table), , ][3, ]
tab2 <- rbind(sm1$q.table[2:7,4,],number.t, number.c, 
              number.all)
row.names(tab2) <- c("Candidate Ideology",
                     "Perception of Party Ideology", "Respondent Ideology",
                     "Respondent Ideology * Feeling Thermometer",
                     "Feeling Thermometer", "Political Awareness",
                     "No. treated", "No. control", "N")
xtable(tab2)

#plotting figure
#setwd("c:/match/docs/papers/koch/writeup")
setwd("c:/R/match/docs/papers/koch/writeup")
trellis.device(device="pdf",file="kochdens.pdf",color=FALSE,width=6,height=4)
par(mar=c(2, 2, 2, 2) + 0.1, cex.lab=0.6, cex.axis=0.6,
    mgp=c(1,0.5,0), cex.main=0.5, cex=0.8)
doverlay(mcoef,coef,lwd=2,
         xlab="Estimated Average Treatment Effect", leg=F)
arrows(coefficients(res)[tt], 4.7, coefficients(res)[tt],0, length=0.1)
text(-0.532,3,"Full Data")
text(-0.25,8.85,"Matched\nData")
text(-0.3,5.4,"Point Estimate \n of Full Data")
dev.off()

nn <- 0
for(i in 1:6){
  nn <- nn + choose(6,i)}
