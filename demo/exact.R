###
### An Example Script for Exact Matching
###

data(lalonde)
user.prompt()

m.out <- matchit(treat ~ educ + black + hispan, data = lalonde,
                 method = "exact")
user.prompt()

print(m.out)
user.prompt()

summary(m.out)
