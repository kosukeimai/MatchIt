###
### An Example Script for Full Matching
###

## load the Lalonde data
data(lalonde)
user.prompt()

## conduct full matching using the propensity score based on logistic regression
m.out <- matchit(treat ~ age + educ + black + hispan + married +
                 nodegree + re74 + re75, data = lalonde,
                 method = "full", distance = "logit")
user.prompt()

## print a short summary
print(m.out)
user.prompt()

## balance diagnostics through statistics
summary(m.out)
user.prompt()

## balance diagnostics through graphics
plot(m.out)
