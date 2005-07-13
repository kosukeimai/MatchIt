###
### An Example Script for Exact Matching
###

data(lalonde)

m.out <- matchit(treat ~ educ + black + hispan, data = lalonde,
                 method = "exact")

print(m.out)

summary(m.out)
