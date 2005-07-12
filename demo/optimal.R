###
### An Example Script for Optimal Matching
###
data(lalonde)
m.out <- matchit(treat ~ re74+re75+age+educ, data=lalonde,
                 method = "optimal", distance = "logit", ratio=2)
print(m.out)
summary(m.out)
