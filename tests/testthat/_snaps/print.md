# print.matchit: printed output is stable

    Code
      print(matchit(f_pr, data = lalonde))
    Output
      A `matchit` object
       - method: 1:1 nearest neighbor matching without replacement
       - distance: Propensity score
                   - estimated with logistic regression
       - number of obs.: 614 (original), 370 (matched)
       - target estimand: ATT
       - covariates: age, educ, race, re74

---

    Code
      print(matchit(f_pr, data = lalonde, method = NULL))
    Output
      A `matchit` object
       - method: None (no matching)
       - distance: Propensity score
                   - estimated with logistic regression
       - number of obs.: 614 (original)
       - target estimand: ATT
       - covariates: age, educ, race, re74

---

    Code
      print(matchit(f_pr, data = lalonde, mahvars = ~ age + educ, caliper = c(0.1,
        age = 2), discard = "control", replace = TRUE, ratio = 2))
    Output
      A `matchit` object
       - method: 2:1 nearest neighbor matching with replacement
       - distance: Mahalanobis [matching]
                   Propensity score [caliper, common support]
                   - estimated with logistic regression
       - caliper: <distance> (0.028), age (18.316)
       - common support: control units dropped
       - number of obs.: 614 (original), 303 (matched)
       - target estimand: ATT
       - covariates: age, educ, race, re74

---

    Code
      print(matchit(treat ~ race + married, data = lalonde, method = "exact"))
    Output
      A `matchit` object
       - method: Exact matching
       - number of obs.: 614 (original), 614 (matched)
       - target estimand: ATT
       - covariates: race, married

---

    Code
      print(matchit(f_pr, data = lalonde, method = "subclass", subclass = 4))
    Output
      A `matchit` object
       - method: Subclassification (4 subclasses)
       - distance: Propensity score
                   - estimated with logistic regression
       - number of obs.: 614 (original), 614 (matched)
       - target estimand: ATT
       - covariates: age, educ, race, re74

---

    Code
      print(matchit(f_pr, data = lalonde, method = "cardinality"))
    Output
      A `matchit` object
       - method: Cardinality matching
       - number of obs.: 614 (original), 236 (matched)
       - target estimand: ATT
       - covariates: age, educ, race, re74

---

    Code
      print(matchit(f_pr, data = lalonde, distance = "robust_mahalanobis"))
    Output
      A `matchit` object
       - method: 1:1 nearest neighbor matching without replacement
       - distance: Robust Mahalanobis
       - number of obs.: 614 (original), 370 (matched)
       - target estimand: ATT
       - covariates: age, educ, race, re74

---

    Code
      print(matchit(f_pr, data = lalonde, s.weights = rep(c(1, 2), length.out = nrow(
        lalonde))))
    Output
      A `matchit` object
       - method: 1:1 nearest neighbor matching without replacement
       - distance: Propensity score
                   - estimated with logistic regression
                   - sampling weights included in estimation
       - number of obs.: 614 (original), 370 (matched)
       - sampling weights: present
       - target estimand: ATT
       - covariates: age, educ, race, re74

