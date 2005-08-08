###
### An Example Script for Subclassification
###

## load the Lalonde data
data(lalonde)
user.prompt()

## sublclassification
m.out <-  matchit(treat ~ re74 + re75 + educ + black + hispan + age,
                  data = lalonde, method = "subclass")
user.prompt()

## a short summary
print(m.out)
user.prompt()

## balance diagnostics
summary(m.out)
user.prompt()

## balance diagnostics through plots
plot(m.out)
plot(m.out, type="jitter")
