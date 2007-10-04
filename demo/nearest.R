###
### An Example Script for Nearest Neighbor Matching
###
data(lalonde)
user.prompt()

## 1:1 Nearest neighbor matching
m.out <- matchit(treat ~ re74 + re75 + educ + black + hispan + age,
                 data = lalonde, method = "nearest")

## print a short summary
print(m.out)
user.prompt()

## balance diagnostics through statistics
s.out <- summary(m.out, standardize=TRUE)
print(s.out)
user.prompt()

## balance diagnostics through graphics
plot(m.out)
user.prompt()
plot(m.out, type="jitter")
user.prompt()
plot(m.out, type="hist")
user.prompt()
plot(s.out)

## 2:1 Nearest neighbor matching
m.out1 <- matchit(treat ~ re74+re75+age+educ, data=lalonde,
                  method = "nearest", distance = "logit", ratio=2)
user.prompt()

## print a short summary
print(m.out1)
user.prompt()

## balance diagnostics through statistics
print(summary(m.out1))
user.prompt()

## balance diagnostics through graphics
plot(m.out)
user.prompt()
plot(m.out, type="jitter")

## 1:1 Nearest neighbor matching with Mahalanobis matching on re74 and re75 and exact matching on married
m.out2 <- matchit(treat ~ re74+re75+age+educ, data=lalonde,
                 method = "nearest", distance = "logit", mahvars=c("re74", "re75"), exact=c("married"), caliper=.25)
user.prompt()

## print a short summary
print(m.out2)
user.prompt()

## balance diagnostics through statistics
s.out2 <- summary(m.out2, standardize=TRUE)
print(s.out2)
user.prompt()

## balance diagnostics through graphics
plot(m.out2)
user.prompt()
plot(m.out2, type="jitter")
user.prompt()
plot(s.out2)

## 1:1 Nearest neighbor matching with units outside the common support discarded
m.out3 <- matchit(treat ~ re74+re75+age+educ, data=lalonde,
                 method = "nearest", distance = "logit", discard= "both")
user.prompt()

## print a short summary
print(m.out3)
user.prompt()

## balance diagnostics through statistics
print(summary(m.out3))
user.prompt()

## balance diagnostics through graphics
plot(m.out3)
plot(m.out3, type="jitter")

## 2:1 Nearest neighbor matching with replacement
m.out4 <- matchit(treat ~ re74+re75+age+educ, data=lalonde,
                 method = "nearest", distance = "logit", replace=TRUE, ratio=2)
user.prompt()

## print a short summary
print(m.out4)
user.prompt()

## balance diagnostics through statistics
print(summary(m.out4))
user.prompt()

## balance diagnostics through graphics
plot(m.out4)
plot(m.out4, type="jitter")
plot(m.out4, type="hist")

## 1:1 Nearest neighbor matching followed by subclassification
m.out5 <- matchit(treat ~ re74+re75+age+educ, data=lalonde,
                 method = "nearest", distance = "logit", subclass=5)
user.prompt()

## print a short summary
print(m.out5)
user.prompt()

## balance diagnostics through statistics
print(summary(m.out5))
user.prompt()

## balance diagnostics through graphics
plot(m.out5)
user.prompt()

s.out5 <- summary(m.out5, standardize=TRUE)
plot(s.out5)
