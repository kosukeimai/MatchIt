# Examples of analyses, with and without Zelig
# Without zelig

data(lalonde)

# First just neyman on matched sample
print("Neyman estimates using a matched sample")
match.out1 <- matchit(treat ~ age + educ + black + hispan + nodegree +
                      married + re74 + re75, data = lalonde)
neyman.out <- neyman(re78, match.out1)
print(neyman.out)
print(summary(neyman.out))
silent <- readline("\nPress <return> to continue: ")

# Neyman using propensity score subclasses
print("Neyman estimates with propensity score subclassification")
match.out1s <- matchit(treat ~ age + educ + black + hispan + nodegree
                       + married + re74 + re75, data = lalonde, subclass=4)
neyman.outs <- neyman(re78, match.out1s)
print(neyman.outs)
print(summary(neyman.outs))
silent <- readline("\nPress <return> to continue: ")

# Neyman using exact matching
print("Neyman estimates using exact matching")
match.out1e <- matchit(treat ~ black + hispan, data = lalonde, exact=TRUE)
neyman.oute <- neyman(re78, match.out1e)
print(neyman.oute)
print(summary(neyman.oute))
silent <- readline("\nPress <return> to continue: ")

# Run linear regression on matched data: use match.out1
print("Linear regression on matched data set")
lm.1 <- lm(re78 ~ treat + pscore, data=match.data(match.out1))

print(summary(lm.1)$coef[2,])
silent <- readline("\nPress <return> to continue: ")

# Zelig part from zelig documentation (match.R in demo directory)
foo <- .find.package("Zelig",quiet=T)
if(length(foo)==0){
  cat("Zelig is not installed.  \nFor further illustrations of analyses,\nyou may install Zelig from \nhttp://gking.harvard.edu/zelig/\n")
} else {
  library(Zelig)

  ## an example for propensity score matching
  print("Using Zelig for analysis")
  z.out1 <- zelig(re78 ~ pscore, data = match.data(match.out1,
                                   "control"), model = "ls")
  x.out1 <- setx(z.out1, data = match.data(match.out1, "treat"),
                 fn = NULL, cond = TRUE)
  s.out1 <- sim(z.out1, x = x.out1)
  print(summary(s.out1))
  user.prompt()
  
  ## an example for subclassification
  print("Using Zelig for analysis, with propensity score subclassification")
  z.out2 <- zelig(re78 ~ pscore, data = match.data(match.out1s,
                                   "control"), model="ls", by="psclass")
  x.out2 <- setx(z.out2, data = match.data(match.out1s, "treat"),
                 fn = NULL, cond = TRUE)
  s.out2 <- sim(z.out2, x = x.out2, num = 100)
  print(summary(s.out2)) # overall results
  print(summary(s.out2, subset = 1)) # subclass 1
  print(summary(s.out2, subset = 2)) # subclass 2
  print(summary(s.out2, subset = 3)) # subclass 3
  user.prompt()
}
