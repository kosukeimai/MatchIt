# Demo of propensity score specification diagnostics
# Initial model
data(lalonde)
match.out1 <- matchit(treat ~ age + married + black + hispan +
                      nodegree + educ + re74 + re75, data=lalonde,                      
                      subclass=c(0, .5, 1))

summary(match.out1, verbose=F)

# Create more blocks if t-stat of pscore > 2.5
match.out1 <- matchit(treat ~ age + educ + married + nodegree + hispan
                      + black + re74 + re75, data=lalonde,                      
                      subclass=c(0,.25, .5,.75,1))

summary(match.out1, verbose=T)

# Add terms that are significant in > 1 block
match.out2 <- matchit(treat ~ age + educ + married + nodegree + hispan
                      + black + re74 + re75 + I(educ*black) +
                      I(black*re74) + I(age*married) + I(educ*married)
                      + I(educ*nodegree) + I(married*nodegree) +
                      I(married*black) + I(age*nodegree),
                      data=lalonde, subclass=c(0,.5,1))

summary(match.out2, verbose=F)

match.out2 <- matchit(treat ~ age + educ + married + nodegree + hispan
                      + black + re74 + re75 + I(educ*black) +
                      I(black*re74) + I(age*married) + I(educ*married)
                      + I(educ*nodegree) + I(married*nodegree) +
                      I(married*black) + I(age*nodegree),
                      data=lalonde, subclass=c(0,.25,.5,1))

summary(match.out2, verbose=F)
 
