###
### An Example Script for Full Matching
###
data(lalonde)
m.out <- matchit(treat ~ age + educ + black + hispan + married +
                 nodegree + re74 + re75, data = lalonde,
                 method = "full", distance = "logit")
print(m.out)
summary(m.out)
