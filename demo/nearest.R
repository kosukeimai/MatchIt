###
### An Example Script for Nearest Neighbor Matching
###
data(lalonde)
cat("DEMO 2:1 Nearest neighbor matching \n")
m.out <- matchit(treat ~ re74+re75+age+educ, data=lalonde,
                 method = "nearest", distance = "logit", ratio=2)
print(m.out)
summary(m.out)

user.prompt()
cat("DEMO 1:1 Nearest neighbor matching with Mahalanobis matching on re74 and re75 and exact matching on married \n")
m.out2 <- matchit(treat ~ re74+re75+age+educ, data=lalonde,
                 method = "nearest", distance = "logit", mahvars=c("re74", "re75"), exact=c("married"), caliper=.25)
print(m.out2)
summary(m.out2)

user.prompt()
cat("DEMO 1:1 Nearest neighbor matching with units outside the common support discarded \n")
m.out3 <- matchit(treat ~ re74+re75+age+educ, data=lalonde,
                 method = "nearest", distance = "logit", discard="both")
print(m.out3)
summary(m.out3)

user.prompt()
cat("DEMO 2:1 Nearest neighbor matching with replacement \n")
m.out4 <- matchit(treat ~ re74+re75+age+educ, data=lalonde,
                 method = "nearest", distance = "logit", replace=TRUE, ratio=2)
print(m.out4)
summary(m.out4)

user.prompt()
cat("DEMO 1:1 Nearest neighbor matching followed by subclassification \n")
m.out5 <- matchit(treat ~ re74+re75+age+educ, data=lalonde,
                 method = "nearest", distance = "logit", subclass=5)
print(m.out5)
summary(m.out5)


