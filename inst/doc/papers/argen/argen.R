# Argentinian Independence Data

#processing data 
rm(list=ls())
setwd("C:/R/match/docs/papers/argen")
library(matchit)
library(foreign)
dta <- read.dta("Hoaggregate.dta")
dta <- dta[dta$nacional==1,]

# functions
calc.ate <- function(x){
  covars2 <- eval(x$data, sys.parent())
  ttt <- eval(x$treat,covars2)
  weighted.var <- function(x, w) sum(w * (x - weighted.mean(x,w))^2)/(sum(w) - 1)
  n1 <- sum(x$psweights!=0 & ttt==1)
  n0 <- sum(x$psweights!=0 & ttt==0)
  wy1 <- weighted.mean(covars2$supercon[ttt==1],x$psweights[ttt==1])
  wy0 <- weighted.mean(covars2$supercon[ttt==0],x$psweights[ttt==0])
  ate <- wy1-wy0
  vary1 <- weighted.var(covars2$supercon[ttt==1],x$psweights[ttt==1])/n1
  vary0 <- weighted.var(covars2$supercon[ttt==0],x$psweights[ttt==0])/n0
  sdate <- sqrt(vary1+vary0)
  cat("ATE is:",round(ate,2))
  cat("\nSD(ATE) is:",round(sdate,2))
  cat("\n")
  xx <- c(ate,sdate,n1,n0)
}

sub.anal <- function(ww1){
  covars2 <- eval(ww1$data, sys.parent())
  ttt <- eval(ww1$treat,covars2)
  qq <- max(ww1$qindex)
  # ATE, SD, N1, N0, covariates
  xn <- length(attr(delete.response(terms(ww1$formula)),"term.labels"))
  sub.ate <- matrix(0,qq,4+xn)
  for(i in 1:qq){
    sub.ate[i,5:ncol(sub.ate)] <-  ww1$q.table[2:(xn+1),1,i]
    y1 <- covars2$supercon[ww1$qindex==i & ttt==1]
    y0 <- covars2$supercon[ww1$qindex==i & ttt==0]
    sub.ate[i,3] <- n1 <- length(y1)
    if(length(y0)==0){
      sub.ate[i,4] <- 0
      sub.ate[i,1] <-  sub.ate[i,2] <- NA
    } else {
      sub.ate[i,4] <- n0 <- length(y0)
      sub.ate[i,1] <- mean(y1)-mean(y0)
      sub.ate[i,2] <- sqrt(var(y1)/n1+var(y0)/n0)
    }
  }
  sub.ate <- na.omit(sub.ate)
  sub.ate <- as.data.frame(sub.ate)
  names(sub.ate)[5:ncol(sub.ate)] <- row.names(ww1$q.table)[2:(xn+1)]
  names(sub.ate)[1:4] <- c("ATE","SD(ATE)","N_1","N_0")
  return(sub.ate)
}

#======================================================================
#Replication of Table 2, Eq. 1a
t2eq1a <- glm(supercon~unisuper+unisimp+gobdiv+mc50+ley,
              family=binomial(link=logit), data=dta)
summary(t2eq1a)

#Matchit model for divided government
dta2 <- subset.data.frame(dta, select=c(supercon,unisuper,unisimp,gobdiv,mc50,ley))
dta2 <- na.omit(dta2)
foo <- matchit(gobdiv ~ unisuper+unisimp+mc50+ley,
              exact=T, data=dta2)
calc.ate(foo)
round(sub.anal(foo),2)

#======================================================================
#Replication of Table 2, Eq. 2a
t2eq2a <- glm(supercon~unisuper+unisimp+gobdiv+ptenp0+ley,
              family=binomial(link=logit), data=dta)
summary(t2eq2a)

#Matchit model for divided government
dta2 <- subset.data.frame(dta, select=c(supercon,unisuper,unisimp,gobdiv,ptenp0,ley))
dta2 <- na.omit(dta2)
foo <- matchit(gobdiv ~ unisuper+unisimp+ptenp0+ley,
              exact=T, data=dta2)
calc.ate(foo)
round(sub.anal(foo),2)

#=====================================================================
#Replication of Table 2, Eq. 3a
t2eq3a <- glm(supercon~unisuper+unisimp+gobdiv+mc50+ley+lp,
              family=binomial(link=logit), data=dta)
summary(t2eq3a)

#Matchit model for divided government
dta2 <- subset.data.frame(dta, select=c(supercon,unisuper,unisimp,gobdiv,mc50,ley,lp))
dta2 <- na.omit(dta2)
foo <- matchit(gobdiv ~ unisuper+unisimp+mc50+ley+lp,
              exact=T, data=dta2)
calc.ate(foo)
round(sub.anal(foo),2)

#=====================================================================
#Replication of Table 2, Eq. 4a
t2eq4a <- glm(supercon~unisuper+unisimp+gobdiv+mc50+procsupe+procdf+ley,
              family=binomial(link=logit), data=dta)
summary(t2eq4a)

#Matchit model for divided government
dta2 <- subset.data.frame(dta,select=c(supercon,unisuper,unisimp,gobdiv,
                          mc50,ley,procsupe,procdf))
dta2 <- na.omit(dta2)
foo <- matchit(gobdiv ~ unisuper+unisimp+mc50+procsupe+procdf+ley,
              exact=T, data=dta2)
calc.ate(foo)
round(sub.anal(foo),2)

#=====================================================================
#Replication of Table 2, Eq. 5a
t2eq5a <- glm(supercon~unisuper+unisimp+gobdiv+ptenp0+procsupe+procdf+ley,
              family=binomial(link=logit), data=dta)
summary(t2eq5a)

#Matchit model for divided government
dta2 <- subset.data.frame(dta,select=c(supercon,unisuper,unisimp,gobdiv,
                          ptenp0,ley,procsupe,procdf))
dta2 <- na.omit(dta2)
foo <- matchit(gobdiv ~ unisuper+unisimp+ptenp0+procsupe+procdf+ley,
              exact=T, data=dta2)
#no treated units matched when exact matching
#complete separation of data!!
foo <- matchit(gobdiv ~ unisuper+unisimp+ptenp0+procsupe+procdf+ley,
              data=dta2)
plot(foo)

