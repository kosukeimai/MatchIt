###
### An Example Script for Subclassification
###
data(lalonde)
user.prompt()

m.out <- matchit(treat ~ re74+re75+age+educ, data=lalonde, method = "subclass")
user.prompt()

print(m.out)
user.prompt()

summary(m.out)

