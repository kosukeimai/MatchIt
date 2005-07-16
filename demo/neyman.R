###
### Example Script for Estimating the Average Treatment Effect through
### Neyman's Method
###

## load the Lalonde data
data(lalonde)
user.prompt()

## conduct the default propensity score matching
m.out1 <- matchit(treat ~ age + educ + black + hispan + nodegree +
                  married + re74 + re75, data = lalonde)
user.prompt()

## calculate ATE on re78
n.out1 <- neyman(re78, m.out1)
user.prompt()

## a short summary
print(n.out1)
user.prompt()

## a more detailed summary
summary(n.out1) 
user.prompt()


## subclassification
m.out2 <- matchit(treat ~ age + educ + black + hispan + nodegree +
                  married + re74 + re75, data = lalonde, subclass=4)
user.prompt()

n.out2 <- neyman(re78, m.out2)
user.prompt()

print(n.out2)
user.prompt()

summary(n.out2)
user.prompt()


## exact matching
m.out3 <- matchit(treat ~ black + hispan, data = lalonde, method = "exact")
user.prompt()

n.out3 <- neyman(re78, m.out3)
user.prompt()

print(n.out3)
user.prompt()

summary(n.out3)
