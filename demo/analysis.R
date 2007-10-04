###
### Example 1: calculating the average treatment effect for the treated
###

## load the Lalonde data
data(lalonde)

## load Zelig package: if not already installed, try install.package("Zelig")
library(Zelig)

## propensity score matching
m.out1 <- matchit(treat ~ age + educ + black + hispan + nodegree + married + re74 + re75, 
                  method = "nearest", data = lalonde)
user.prompt()

## fit the linear model to the control group controlling for propensity score and 
## other covariates
z.out1 <- zelig(re78 ~ age + educ + black + hispan + nodegree + married + re74 + re75 +
                       distance, data = match.data(m.out1, "control"), model = "ls")
user.prompt()

## set the covariates to the covariates of matched treated units
## use conditional prediction by setting cond = TRUE.
x.out1 <- setx(z.out1, data = match.data(m.out1, "treat"), fn = NULL, cond = TRUE)
user.prompt()

## simulate quantities of interest
s.out1 <- sim(z.out1, x = x.out1)
user.prompt()

## obtain a summary
print(summary(s.out1))
user.prompt()


###
### Example 2: calculating the average treatment effect for the entire sample
###

## fit the linear model to the treatment group controlling for propensity score and 
## other covariates
z.out2 <- zelig(re78 ~ age + educ + black + hispan + nodegree + married + re74 + re75 +
                       distance, data = match.data(m.out1, "treat"), model = "ls")
user.prompt()

## conducting the simulation procedure for the control group
x.out2 <- setx(z.out2, data = match.data(m.out1, "control"), fn = NULL, cond = TRUE)
user.prompt()

s.out2 <- sim(z.out2, x = x.out2)
user.prompt()

##  Note that Zelig calculates the difference between observed and
##  either predicted or expected values.  This means that the treatment
##  effect for the control units is actually the effect of control
##  (observed control outcome minus the imputed outcome under treatment
##  from the model).  Hence, to combine treatment effects just reverse
##  the signs of the estimated treatment effect of controls.
ate.all <- c(s.out1$qi$att.ev, -s.out2$qi$att.ev)
user.prompt()

## some summaries
## point estimate
print(mean(ate.all))
user.prompt()
## standard error
print(sd(ate.all))
user.prompt()
## 95% confidence interval
print(quantile(ate.all, c(0.025, 0.975)))
user.prompt()


###
### Example 3: subclassification
###

## subclassification with 4 subclasses
m.out2 <- matchit(treat ~ age + educ + black + hispan + nodegree + married + re74 + re75,  
                  data = lalonde, method = "subclass", subclass = 4)
user.prompt()

## controlling only for the estimated prpensity score and lagged Y within each subclass
## one can potentially control for more
z.out3 <- zelig(re78 ~ re74 + re75 + distance, data = match.data(m.out2, "control"), 
                model = "ls", by = "subclass")
user.prompt()

## conducting simulations
x.out3 <- setx(z.out3, data = match.data(m.out2, "treat"), fn = NULL, cond = TRUE)
user.prompt()

## for the demonstration purpose, we set the number of simulations to be 100
s.out3 <- sim(z.out3, x = x.out3, num = 100)
user.prompt()

## overall results
print(summary(s.out3)) 
user.prompt()

## summary for each subclass
print(summary(s.out3, subset = 1))
user.prompt()

print(summary(s.out3, subset = 2))
user.prompt()

print(summary(s.out3, subset = 3))


