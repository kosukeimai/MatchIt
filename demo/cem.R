###
### An Example Script for Coarsened Exact Matching
###

## load the Lalonde data
data(lalonde)

## coarsened exact matching with automatic coarsening
m.out <- matchit(treat ~ age + educ + black + hispan + married + nodegree 
                 + re74 + re75, data = lalonde, method = "cem")
user.prompt()

## print a short summary
print(m.out)
user.prompt()

## balance diagnostics through statistics; standardize = T for plotting
s.out <- summary(m.out, covariates = T, standardize = T)
print(s.out)
user.prompt()

## graphical balance checks
plot(m.out)
plot(s.out)


## create some cutpoints for continuous variables
re74cut <- hist(lalonde$re74, br=seq(0,max(lalonde$re74)+1000, by=1000),plot=FALSE)$breaks
re75cut <- hist(lalonde$re75, br=seq(0,max(lalonde$re75)+1000, by=1000),plot=FALSE)$breaks
agecut <- hist(lalonde$age, br=seq(15,55, length=14),plot=FALSE)$breaks
mycp <- list(re75=re75cut, re74=re74cut, age=agecut)

## coarsened exact matching with user-given cutpoints
m.out2 <- matchit(treat ~ age + educ + black + hispan + married + nodegree 
                  + re74 + re75, data = lalonde, method = "cem",
                  cutpoints = mycp)
user.prompt()

## print a short summary
print(m.out2)
user.prompt()

## balance diagnostics through statistics
s.out2 <- summary(m.out2, covariates = T)
print(s.out2)
user.prompt()

