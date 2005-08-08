###
### An Example Script for Nearest Neighbor Matching
###
data(lalonde)
user.prompt()

## 1:1 Nearest neighbor matching
m.out <- matchit(treat ~ re74 + re75 + educ + black + hispan + age,
                 data = lalonde, method = "nearest")
user.prompt()

## print a short summary
print(m.out)
user.prompt()

## balance diagnostics through statistics
summary(m.out)
user.prompt()

## balance diagnostics through graphics
plot(m.out)

## 2:1 Nearest neighbor matching
m.out <- matchit(treat ~ re74+re75+age+educ, data=lalonde,
                 method = "nearest", distance = "logit", ratio=2)
user.prompt()

## print a short summary
print(m.out)
user.prompt()

## balance diagnostics through statistics
summary(m.out)
user.prompt()

## balance diagnostics through graphics
plot(m.out)

## 1:1 Nearest neighbor matching with Mahalanobis matching on re74 and re75 and exact matching on married
m.out2 <- matchit(treat ~ re74+re75+age+educ, data=lalonde,
                 method = "nearest", distance = "logit", mahvars=c("re74", "re75"), exact=c("married"), caliper=.25)
user.prompt()

## print a short summary
print(m.out2)
user.prompt()

## balance diagnostics through statistics
summary(m.out2)
user.prompt()

## balance diagnostics through graphics
plot(m.out2)

## 1:1 Nearest neighbor matching with units outside the common support discarded
m.out3 <- matchit(treat ~ re74+re75+age+educ, data=lalonde,
                 method = "nearest", distance = "logit", discard="both")
user.prompt()

## print a short summary
print(m.out3)
user.prompt()

## balance diagnostics through statistics
summary(m.out3)
user.prompt()

## balance diagnostics through graphics
plot(m.out3)

## 2:1 Nearest neighbor matching with replacement
m.out4 <- matchit(treat ~ re74+re75+age+educ, data=lalonde,
                 method = "nearest", distance = "logit", replace=TRUE, ratio=2)
user.prompt()

## print a short summary
print(m.out4)
user.prompt()

## balance diagnostics through statistics
summary(m.out4)
user.prompt()

## balance diagnostics through graphics
plot(m.out4)

## 1:1 Nearest neighbor matching followed by subclassification
m.out5 <- matchit(treat ~ re74+re75+age+educ, data=lalonde,
                 method = "nearest", distance = "logit", subclass=5)
user.prompt()

## print a short summary
print(m.out5)
user.prompt()

## balance diagnostics through statistics
summary(m.out5)
user.prompt()

## balance diagnostics through graphics
plot(m.out5)
