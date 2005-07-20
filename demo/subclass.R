###
### An Example Script for Subclassification
###

## load the Lalonde data
data(lalonde)
user.prompt()

## sublcalssification
m.out <- matchit(treat ~ re74+re75+age+educ, data=lalonde, method = "subclass")
user.prompt()

## a short summary
print(m.out)
user.prompt()

## balance diagnostics
summary(m.out)

