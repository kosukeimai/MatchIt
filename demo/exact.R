###
### An Example Script for Exact Matching
###

## laod the Lalonde data
data(lalonde)
user.prompt()

## exact matching
m.out <- matchit(treat ~ educ + black + hispan, data = lalonde,
                 method = "exact")
user.prompt()

## a short summary
print(m.out)
user.prompt()

## a full summary
summary(m.out)
