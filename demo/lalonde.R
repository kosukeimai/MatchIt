# MatchIt Demo script with Lalonde data

data(lalonde)

#Exact matching
foo1 <- matchit(treat ~ black + hispan, exact=TRUE, data=lalonde)
print(foo1)
summary(foo1)
summary(foo1,verbose=T)
row.names(foo1$data)[foo1$psclass==1][1:20]
foo1$data[foo1$psclass==2,c("black", "hispan")][1:20,]
silent <- readline("\nPress <return> to continue: ")

#Propensity score matching
foo2 <- matchit(treat ~ re74 + re75, data=lalonde)
print(foo2)
summary(foo2)
summary(foo2, verbose=T)
plot(foo2)

# Greedy vs. optimal matching
if ("optmatch"%in%.packages(all=T)) {
foo <- matchit(treat ~ re74+re75+age+educ, data=lalonde)
print(summary(foo))
foo <- matchit(treat ~ re74+re75+age+educ, data=lalonde, m.order=1)
print(summary(foo))
}

#Propensity score with exact restriction
foo3 <- matchit(treat ~ re74 + re75, exact=c("black","hispan"), data=lalonde)
print(foo3)
sum.foo3 <- summary(foo3)
sum.foo3$sum.all[c("hispan","black"),]
sum.foo3$sum.matched[c("hispan","black"),]
dta.matched <- match.data(foo3)
dta.treated <- match.data(foo3, group = "treat")
dta.control <- match.data(foo3, group = "control")

#Replication of Dehejia & Wahba
foo <- matchit(treat ~ age + I(age^2) + educ + I(educ^2) + black +
               hispan + married + nodegree + re74 + I(re74^2) +
               re75 + I(re75^2) + I(as.numeric(re74==0)) +
               I(as.numeric(re75==0)), data=lalonde, replace=TRUE,
               discard=1)

#Non-matched and discarded units
foo$in.sample[!foo$in.sample]
foo$psweights[foo$in.sample & foo$psweights==0]
sdisc <- apply(lalonde[!foo$in.sample,],2,mean) #taking column-wise means
streat <- apply(lalonde[lalonde$treat==1,],2,mean)
round(rbind(sdisc,streat),2) #rounding

#Look at weights
print(foo$psweights[foo$data$treat==1][1:10])
print(foo$psweights[foo$data$treat==0][1:20])
print((foo$psweights[foo$data$treat==0]*sum(foo$psweights[foo$data$treat==1])/sum(foo$psweights[foo$data$treat==0]))[1:20])

#Subclassification
foo1 <- matchit(treat ~ age + educ + black + hispan + married +
                nodegree + re74 + re75, data=lalonde, replace=TRUE,
                subclass=6)
foo2 <- matchit(treat ~ age + educ + black + hispan + married +
                nodegree + re74 + re75, data=lalonde, nearest=FALSE,
                replace=TRUE, subclass=6)
print(foo2)
plot(foo2)

#Exact matching on all covariates
foo <- matchit(treat ~ age + educ + black + hispan + married +
               nodegree + re74 + re75, data=lalonde, exact=TRUE)
summary(foo)

# Full matching
if ("optmatch"%in%.packages(all=T)) {
foo1 <- matchit(treat ~ age + educ + black + hispan + married +
               nodegree + re74 + re75, data=lalonde, full=T)
print(summary(foo1))
}

# Caliper matching
foo <- matchit(treat ~ age + educ + black + hispan + married +
               nodegree + re74 + re75, data=lalonde, caliper=0.25,
               replace=TRUE)

#Mahalanobis matching
foo <- matchit(treat ~ age + educ, data=lalonde, mahvars = c("age", "educ"),
               caliper=Inf, replace=TRUE)

#Assignment Model specification
foo1 <- matchit(treat~educ+re74+re75,data=lalonde,model="probit")
foo2 <- matchit(treat~educ+re74+re75,data=lalonde,model="nnet")
foo3 <- matchit(treat~educ+re74+re75,data=lalonde,model="GAM")
foo4 <- matchit(treat~educ+re74+re75,data=lalonde,model="cart")

#Using Observation Names
test <- lalonde[1:4,] #taking a lalonde subset
row.names(test) <- c("Dan","Kosuke","Liz","Gary") #assigning row names
print(test)

#Matching on One Covariate
index <- rep(TRUE,nrow(lalonde)) #keeping all units
names(index) <- row.names(lalonde) #assigning obs names
treat <- lalonde$treat
names(treat) <- row.names(lalonde)
re75 <- lalonde$re75
names(re75) <- row.names(lalonde)
foo <- matchdef(treat~re75,index,pscore=re75,replace=TRUE,data=lalonde)
cbind(lalonde[row.names(foo$match.matrix),]$re75,
      lalonde[as.vector(foo$match.matrix[,1]),]$re75) #matched re75
plot(lalonde[row.names(foo$match.matrix),]$re75,
     lalonde[as.vector(foo$match.matrix[,1]),]$re75, pch=16,
     xlab="Treated", ylab="Controls") #plotting matched data
