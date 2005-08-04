###
### An Example Script for Obtaining Mathced Data
###

## load the Lalonde data
data(lalonde)
user.prompt()

## perform nearest neighbor matching
m.out1 <- matchit(treat ~ re74 + re75 + age + educ, data = lalonde,
                 method = "nearest", distance = "logit")
user.prompt()

## obtain matched data
m.data1 <- match.data(m.out1)
user.prompt()

## summarize the resulting matched data
summary(m.data1)
user.prompt()

## obtain matched data for the treatment group
m.data2 <- match.data(m.out1, group = "treat")
user.prompt()

summary(m.data2)
user.prompt()

## obtain matched data for the control group
m.data3 <- match.data(m.out1, group = "control")
user.prompt()

summary(m.data3)
user.prompt()

## run a subclassification method
m.out2 <- matchit(treat ~ re74 + re75 + age + educ, data=lalonde, method = "subclass")
user.prompt()

## specify different names
m.data4 <- match.data(m.out2, subclass = "block", weights = "w",
                      distance = "pscore")
user.prompt()

## print the variable names of the matched data
names(m.data4)
