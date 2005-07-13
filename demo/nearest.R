###
### An Example Script for Nearest Neighbor Matching
###
data(lalonde)
m.out <- matchit(treat ~ re74+re75+age+educ, data=lalonde,
                 method = "nearest", distance = "logit", ratio=2)
cat("DEMO 2:1 Nearest neighbor matching \n")
print(m.out)
summary(m.out)

m.out2 <- matchit(treat ~ re74+re75+age+educ, data=lalonde,
                 method = "nearest", distance = "logit", mahvars=c("re74", "re75"), exact=c("married"))
cat("DEMO 1:1 Nearest neighbor matching with Mahalanobis matching on re74 and re75 and exact matching on married \n")
print(m.out2)
summary(m.out2)

m.out3 <- matchit(treat ~ re74+re75+age+educ, data=lalonde,
                 method = "nearest", distance = "logit", discard="both")
cat("DEMO 1:1 Nearest neighbor matching with units outside the common support discarded \n")
print(m.out3)
summary(m.out3)

m.out4 <- matchit(treat ~ re74+re75+age+educ, data=lalonde,
                 method = "nearest", distance = "logit", replace=TRUE, ratio=2)
cat("DEMO 2:1 Nearest neighbor matching with replacement \n")
print(m.out4)
summary(m.out4)

m.out5 <- matchit(treat ~ re74+re75+age+educ, data=lalonde,
                 method = "nearest", distance = "logit", subclass=5)
cat("DEMO 1:1 Nearest neighbor matching followed by subclassification \n")
print(m.out5)
summary(m.out5)


