###
### An Example Script for Subclassification
###
data(lalonde)
m.out <- matchit(treat ~ re74+re75+age+educ, data=lalonde, method = "subclass")
print(m.out)
summary(m.out)

