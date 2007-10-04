####
#### demo file for Genetic Matching
#### 

## loading the lalonde data
data(lalonde)

## using logistic propensity score as one of the covariates
m.out <- matchit(treat ~ age + educ + black + hispan + married + nodegree + re74 + re75, data = lalonde, method = "genetic", distance = "logit")
user.prompt()

## printing a short summary
print(m.out)
user.prompt()

## numerical balance diagonstics
s.out <- summary(m.out, standardize=TRUE)
print(s.out)
user.prompt()

## graphical balance diagnostics
plot(m.out)
plot(s.out)
