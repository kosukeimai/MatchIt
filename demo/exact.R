###
### An Example Script for Exact Matching
###

## laod the Lalonde data
data(lalonde)

## exact matching
m.out <- matchit(treat ~ educ + black + hispan, data = lalonde,
                 method = "exact")
user.prompt()

## print a short summary
print(m.out)
user.prompt()

## balance diagnostics through statistics
print(summary(m.out, covariates = T))

