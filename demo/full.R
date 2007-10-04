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
user.prompt()

## balance diagnostics through statistics
s.out <- summary(m.out)
print(s.out)
s.out <- summary(m.out, standardize=TRUE)
print(s.out)
user.prompt()

## balance diagnostics through graphics
plot(m.out)
plot(s.out)
