####
#### demo file for Genetic Matching
#### 

## loading the lalonde data
library(lalonde)

## using logistic propensity score as one of the covariates
m.out <- matchit(treat ~ age + educ + black + hispan + married + nodegree + re74 + re75, data = lalonde, method = "genetic", distance = "logit")

## printing a short summary
print(m.out)

## numerical balance diagonstics
summary(m.out)

## graphical balance diagnostics
plot(m.out)
