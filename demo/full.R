###
### An Example Script for Full Matching
###

## load the Lalonde data
data(lalonde)
## conduct full matching using the propensity score based on logistic regression
m.out <- matchit(treat ~ age + educ + black + hispan + married +
                 nodegree + re74 + re75, data = lalonde,
                 method = "full", distance = "logit")
## print a short summary
print(m.out)
## balance diagnostics through statistics
summary(m.out)
