###
### An Example Script for Optimal Matching
###

## load the Lalonde data
data(lalonde)

## optimal ratio matching using the propensity score based on logistic regression
m.out <- matchit(treat ~ re74 + re75 + age + educ, data = lalonde,
                 method = "optimal", distance = "logit", ratio = 2)
user.prompt()

## a short summary
print(m.out)
user.prompt()

## balance diagnostics through statistics
print(summary(m.out))
user.prompt()

## balance diagnostics through graphics
plot(m.out)
user.prompt()
plot(m.out, type="jitter")

s.out <- summary(m.out, standardize=TRUE)
plot(s.out)

