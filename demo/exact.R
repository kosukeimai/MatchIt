###
### An Example Script for Exact Matching
###

data(lalonde)

m.out <- matchit(treat ~ age + educ + black + hispan + married, data = lalonde,
                 method = "exact")

print(m.out)

summary(m.out)
