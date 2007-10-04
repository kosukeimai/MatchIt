###
### An Example Script for Subclassification
###

## load the Lalonde data
data(lalonde)

## sublclassification
m.out <-  matchit(treat ~ re74 + re75 + educ + black + hispan + age,
                  data = lalonde, method = "subclass")
user.prompt()

## a short summary
print(m.out)
user.prompt()

## balance diagnostics
print(summary(m.out))
user.prompt()

## balance diagnostics through plots
plot(m.out)
user.prompt()
plot(m.out, type="jitter")

s.out <- summary(m.out, standardize=TRUE)
plot(s.out)
