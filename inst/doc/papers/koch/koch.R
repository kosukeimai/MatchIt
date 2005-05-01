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

#==================================
# Model 3 (Republicans)
#==================================

#creating full dataset 
dta.full <- subset(dta,select=c(prcanid, rviswom, rvisman, repcan1 , goppty ,
              rideo , rproj , repft , aware), subset=c(rvisman==0 & voter==1))
dta.full <- na.omit(dta.full)

#matching and creating matched dataset
#now testing sensitivity
#fml <- as.formula(prcanid ~ rviswom + repcan1 + goppty + rideo + rproj + repft + aware)
fml <- as.formula(prcanid ~ rviswom + repcan1 + goppty + rideo + rproj
                  + repft + aware + I(repcan1^2) + I(goppty^2) +
                  I(rideo^2) + I(rproj^2) + I(repft^2) + I(aware^2))
res <- zelig(fml, data=dta.full, model="ls")
xvars <- names(res$coefficients)
xvars <- xvars[3:length(xvars)]
tt <- attr(terms(res),"term.labels")[1]
yy <- attr(terms(res),"variables")[[2]]
mfml <- as.formula(paste(tt,"~",paste(xvars,collapse=" + ")))
m1 <- matchit(mfml, data=dta.full, nearest=T, discard=1,
              caliper=0.3, reestimate=T)
#m1 <- matchit(mfml, data=dta.full, subclass=6, nearest=F, discard=1)
dta.match <- match.data(m1)

#new diagnostics
X1 <- dta.full$rideo[dta.full$rviswom==1]
X0 <- dta.full$rideo[dta.full$rviswom==0]
sdta <- dta.match[dta.match$psclass==1,]
x1 <- sdta$rideo[sdta$rviswom==1]
x0 <- sdta$rideo[sdta$rviswom==0]

X1 <- dta.full$repcan1[dta.full$rviswom==1]
X0 <- dta.full$repcan1[dta.full$rviswom==0]
sdta <- dta.match[dta.match$psclass==1,]
x1 <- sdta$repcan1[sdta$rviswom==1]
x0 <- sdta$repcan1[sdta$rviswom==0]

pdf("kochsub1.pdf",width=8,height=4)
par(mfrow=c(1,2))
plot(ecdf(X1),verticals=T, 
     main="ECDF in Raw Data")
lines(ecdf(X0),verticals=T,lty=2,lwd=2.5)
plot(ecdf(x1),verticals=T, 
     main="ECDF in Subclass 1")
lines(ecdf(x0),verticals=T,lty=2,lwd=2.5)
dev.off()

pdf("kochqq.pdf",width=8,height=4)
par(mfrow=c(1,2))
foo0 <- qqplot(X1,X0,type="l",
       ylab="Treated Respondent Ideology",
       xlab="Control Respondent Ideology",
       main="Raw Data QQ-Plot")
abline(a=0,b=1,col="darkgrey")
foo1 <- qqplot(x1,x0,type="l",
       ylab="Treated Respondent Ideology",
       xlab="Control Respondent Ideology",
       main="Matched Subclass 1 QQ-Plot")
abline(a=0,b=1,col="darkgrey")
dev.off()

mean(abs(foo0$x-foo0$y))
mean(abs(foo1$x-foo1$y))


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
    #for subclass
#    tmp <- zelig(ftmp, data = dta.match, model="ls", by="psclass")
#    wate <- 0
#    for(l in 1:length(tmp)){
#      wate <- wate + tmp[[l]]$coefficient[tt]*sum(dta.match$psclass==l & dta.match[,tt]==1)/sum(dta.match$psclass!=0 & dta.match[,tt]==1)
#    }
#    mcoef[counter] <- wate
    #for caliper matching
    mcoef[counter] <- lm(ftmp, data = dta.match)$coefficient[tt]
    counter <- counter + 1
  }
  cat(i,"covariates:",counter,"\n")
  cat(date(), "\n\n")
}
save.image("koch.Rdata")
#coef[length(coef)] <- res$coefficients[tt]
#mres <- zelig(fml, data=dta.match, model="ls",by="psclass")
#wate <- 0
#for(l in 1:length(tmp)){
#  wate <- wate + tmp[[l]]$coefficient[tt]*sum(dta.match$psclass==l & dta.match[,tt]==1)/sum(dta.match$psclass!=0 & dta.match[,tt]==1)  
}
#mcoef[length(mcoef)] <- wate

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
#setwd("c:/R/match/docs/papers/koch/writeup")
trellis.device(device="pdf",file="kochdens.pdf",color=FALSE,width=6,height=3)
par(mar=c(2, 2, 2, 2) + 0.1, cex.lab=0.6, cex.axis=0.6,
    mgp=c(1,0.5,0), cex.main=0.5, cex=0.8, bg="white")
doverlay(mcoef,coef,lwd=2,
         xlab="Estimated average treatment effect", leg=F)
arrows(coefficients(res)[tt], 4.7, coefficients(res)[tt],0, length=0.1)
text(-0.532,3,"Raw data")
text(-0.25,8.85,"Matched\ndata")
text(0.5,0.5,"Matched Data")
text(-0.3,5.4,"Point estimate \n of raw data")
dev.off()

#traceplots
trellis.device(device="pdf",file="kochtrace.pdf",color=FALSE,width=6,height=3)
par(mar=c(2, 2, 2, 2) + 0.1, cex.lab=0.6, cex.axis=0.6,
    mgp=c(1,0.5,0), cex.main=0.5, cex=0.8, bg="white",par(mfrow=c(2,1)))
ts.plot(coef,ylim=range(c(coef,mcoef)),
        ylab="Treatment Effect",main="Raw Data")
ts.plot(mcoef,ylim=range(c(coef,mcoef)),
        ylab="Treatment Effect",main="Matched Data")
dev.off()

nn <- 0
for(i in 1:6){
  nn <- nn + choose(6,i)}

#other functions
#x is in solid steps
#y is in grey
#y is in dashed steps
histbackback2 <- function (x, y, z, brks = NULL, xlab = NULL, axes = TRUE, probability = FALSE, 
                          xlim = NULL, ylab = "", main="", isf = FALSE, ...) 
{
  if(!isf){
    if (!length(brks)) 
      brks <- hist(c(x, y), plot = FALSE)$breaks
    ll <- hist(x, breaks = brks, plot = FALSE, probability =
               probability)
    rr <- hist(y, breaks = brks, plot = FALSE, probability =
               probability)
    rr2 <- hist(z, breaks = brks, plot = FALSE, probability =
                probability)
    if (probability) {
      ll$counts <- ll$density
      rr$counts <- rr$density
      rr2$counts <- rr2$density
    }
    xl <- pretty(range(c(ll$counts,rr$counts,rr2$counts)))
    xl <- c(xl[1], xl[length(xl)])
    barplot(rr$counts, ylim = xl,  horiz = F, space=0,
              axes = FALSE, col = "darkgrey",border="darkgrey")      
    lines(0:length(brks), c(ll$counts,0,0),type="s",lwd=1)
    lines(0:length(brks), c(rr2$counts,0,0),type="s",lwd=1,lty=2)
  } else {
    xl <- pretty(range(c(table(x)/sum(table(x)),table(y)/sum(table(y)),table(z)/sum(table(z)))))
    xl <- c(xl[1], xl[length(xl)])
    barplot((t(table(y)/sum(table(y)))), ylim = xl, space = 0, horiz = F, 
            axes = FALSE, col = "darkgrey",border="darkgrey", cex.names=0.6)  
    lines(c(0,0:length(table(x))), c(0,as.numeric(t(table(x)/sum(table(x)))),0), type="s",lwd=1)
    lines(c(0,0:length(table(z))), c(0,as.numeric(t(table(z)/sum(table(z)))),0),type="s",lwd=1,lty=2)
  }
  if (axes) {
    axis(2, at = pretty(xl), labels = format(abs(pretty(xl))))
    if(!isf){
      del <- (brks[2] - brks[1] - (brks[3] - brks[2]))/2
      brks[1] <- brks[1] + del
      brks[-1] <- brks[-1] - del
      axis(1, at = 0:(length(brks) - 1), labels = formatC(brks, 
                                               format = "f", digits =
                                               1))
    } 
    if (ylab != "" | main!="") 
      title(ylab = ylab, main=main)
  }
  box()
}
