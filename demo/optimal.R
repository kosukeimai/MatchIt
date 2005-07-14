###
### An Example Script for Optimal Matching
###

## load the Lalonde data
data(lalonde)
## optimal ratio matching using the propensity score based on logistic regression
m.out <- matchit(treat ~ re74 + re75 + age + educ, data = lalonde,
                 method = "optimal", distance = "logit", ratio = 2)
## a short summary
print(m.out)
## balance diagnostics through statistics
summary(m.out)
## balance diagnostics through graphics
plot(m.out)
