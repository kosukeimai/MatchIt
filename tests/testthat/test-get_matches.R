
library(testthat)
library(MatchIt)

data("lalonde")

context("get_matches")

test_that("can get correct matches, exact", {  
  # preliminaries
  m.out <- matchit(treat ~ educ + black + hispan, data = lalonde,
                   method = "exact")
  n_matched <- sum(m.out$nn[2,])
  nms_matched <- names(m.out$subclass[!is.na(m.out$subclass)])
  
  lalonde$id <- 1:nrow(lalonde)
  lalonde2 <- lalonde[sample.int(500, 3000, replace=TRUE),]
  matches <- get_matches(m.out, lalonde)
  matches2 <- get_matches(m.out, model_frame= lalonde, id_cols= "id", newdata= lalonde2)
  
  # tests -- no newdata
  expect_equal(n_matched, nrow(matches))
  expect_equal(sum(matches$treat), m.out$nn[2,2])
  expect_equal(nrow(matches) - sum(matches$treat), m.out$nn[2,1])
  expect_equal(sum(is.na(matches)), 0L)
  expect_equal(names(matches), c(names(lalonde), "weight"))
  expect_equal(ncol(matches), ncol(lalonde) + 1)
  expect_true(all(matches[matches$treat == 1,]$weight == 1))
  expect_true(all(matches$weight %in% m.out$weights))
  expect_true(all(rownames(matches) %in% nms_matched))
  expect_true(all(nms_matched %in% rownames(matches)))
  
  # tests -- newdata
  expect_equal(sum(is.na(matches2)), 0L)
  expect_true(all(names(matches2) %in% c(names(lalonde2), "weight")))
  expect_true(all(c(names(lalonde2), "weight") %in% names(matches2)))
  expect_equal(ncol(matches2), ncol(lalonde2) + 1)
  expect_true(all(matches2[matches2$treat == 1,]$weight == 1))
  expect_true(all(matches2$weight %in% m.out$weights))
})

test_that("can get correct matches, nearest", {
  m.out <- matchit(treat ~ re74 + re75 + educ + black + hispan + age,
                   data = lalonde, method = "nearest")
  n_matched <- sum(m.out$nn[2,])
  nms_matched <- c(rownames(m.out$match.matrix), c(m.out$match.matrix))
  
  lalonde$id <- 1:nrow(lalonde)
  lalonde2 <- lalonde[sample.int(500, 3000, replace=TRUE),]
  matches <- get_matches(m.out, lalonde)
  matches2 <- get_matches(m.out, model_frame= lalonde, id_cols= "id", newdata= lalonde2)
  
  # tests -- no newdata
  expect_equal(n_matched, nrow(matches))
  expect_equal(sum(matches$treat), m.out$nn[2,2])
  expect_equal(nrow(matches) - sum(matches$treat), m.out$nn[2,1])
  expect_equal(sum(is.na(matches)), 0L)
  expect_equal(names(matches), c(names(lalonde), "weight"))
  expect_equal(ncol(matches), ncol(lalonde) + 1)
  expect_true(all(matches[matches$treat == 1,]$weight == 1))
  expect_true(all(matches$weight %in% m.out$weights))
  expect_true(all(rownames(matches) %in% nms_matched))
  expect_true(all(nms_matched %in% rownames(matches)))
  
  # tests -- newdata
  expect_equal(sum(is.na(matches2)), 0L)
  expect_true(all(names(matches2) %in% c(names(lalonde2), "weight")))
  expect_true(all(c(names(lalonde2), "weight") %in% names(matches2)))
  expect_equal(ncol(matches2), ncol(lalonde2) + 1)
  expect_true(all(matches2[matches2$treat == 1,]$weight == 1))
  expect_true(all(matches2$weight %in% m.out$weights))
  
})


test_that("can get correct matches, full", {
  # preliminaries
  m.out <- matchit(treat ~ age + educ + black + hispan + married +
                     nodegree + re74 + re75, data = lalonde, method = "full")
  n_matched <- sum(m.out$nn[2,])
  nms_matched <- names(m.out$subclass[!is.na(m.out$subclass)])
  
  lalonde$id <- 1:nrow(lalonde)
  lalonde2 <- lalonde[sample.int(500, 3000, replace=TRUE),]
  matches <- get_matches(m.out, lalonde)
  matches2 <- get_matches(m.out, model_frame= lalonde, id_cols= "id", newdata= lalonde2)
  
  # tests -- no newdata
  expect_equal(n_matched, nrow(matches))
  expect_equal(sum(matches$treat), m.out$nn[2,2])
  expect_equal(nrow(matches) - sum(matches$treat), m.out$nn[2,1])
  expect_equal(sum(is.na(matches)), 0L)
  expect_equal(names(matches), c(names(lalonde), "weight"))
  expect_equal(ncol(matches), ncol(lalonde) + 1)
  expect_true(all(matches[matches$treat == 1,]$weight == 1))
  expect_true(all(matches$weight %in% m.out$weights))
  expect_true(all(rownames(matches) %in% nms_matched))
  expect_true(all(nms_matched %in% rownames(matches)))
  
  # tests -- newdata
  expect_equal(sum(is.na(matches2)), 0L)
  expect_true(all(names(matches2) %in% c(names(lalonde2), "weight")))
  expect_true(all(c(names(lalonde2), "weight") %in% names(matches2)))
  expect_equal(ncol(matches2), ncol(lalonde2) + 1)
  expect_true(all(matches2[matches2$treat == 1,]$weight == 1))
  expect_true(all(matches2$weight %in% m.out$weights))
})


test_that("can get correct matches, optimal", {
  ratio <- 2
  m.out <- matchit(treat ~ re74 + re75 + age + educ, data = lalonde,
                   method = "optimal", ratio = ratio)
  n_matched <- sum(m.out$nn[2,])
  nms_matched <- c(rownames(m.out$match.matrix), c(m.out$match.matrix))
  
  lalonde$id <- 1:nrow(lalonde)
  lalonde2 <- lalonde[sample.int(500, 3000, replace=TRUE),]
  matches <- get_matches(m.out, lalonde)
  matches2 <- get_matches(m.out, model_frame= lalonde, id_cols= "id", newdata= lalonde2)
  
  # tests -- no newdata
  expect_equal(n_matched, nrow(matches))
  expect_equal(sum(matches$treat), m.out$nn[2,2])
  expect_equal(nrow(matches) - sum(matches$treat), m.out$nn[2,1])
  expect_equal(sum(is.na(matches)), 0L)
  expect_equal(names(matches), c(names(lalonde), "weight"))
  expect_equal(ncol(matches), ncol(lalonde) + 1)
  expect_true(all(matches[matches$treat == 1,]$weight == 1))
  expect_true(all(matches$weight %in% c(1, 1 / ratio)))
  expect_true(all(rownames(matches) %in% nms_matched))
  expect_true(all(nms_matched %in% rownames(matches)))
  
  # tests -- newdata
  expect_equal(sum(is.na(matches2)), 0L)
  expect_true(all(names(matches2) %in% c(names(lalonde2), "weight")))
  expect_true(all(c(names(lalonde2), "weight") %in% names(matches2)))
  expect_equal(ncol(matches2), ncol(lalonde2) + 1)
  expect_true(all(matches2[matches2$treat == 1,]$weight == 1))
  expect_true(all(matches2$weight %in% c(1, 1 / ratio)))
})


test_that("can get correct matches, subclass", {
  m.out <- matchit(treat ~ re74 + re75 + educ + black + hispan + age,
                   data = lalonde, method = "subclass")
  n_matched <- sum(m.out$nn[2,])
  nms_matched <- names(m.out$subclass[!is.na(m.out$subclass)])
  
  lalonde$id <- 1:nrow(lalonde)
  lalonde2 <- lalonde[sample.int(500, 3000, replace=TRUE),]
  matches <- get_matches(m.out, lalonde)
  matches2 <- get_matches(m.out, model_frame= lalonde, id_cols= "id", newdata= lalonde2)
  
  # tests -- no newdata
  expect_equal(n_matched, nrow(matches))
  expect_equal(sum(matches$treat), m.out$nn[2,2])
  expect_equal(nrow(matches) - sum(matches$treat), m.out$nn[2,1])
  expect_equal(sum(is.na(matches)), 0L)
  expect_equal(names(matches), c(names(lalonde), "weight"))
  expect_equal(ncol(matches), ncol(lalonde) + 1)
  expect_true(all(matches[matches$treat == 1,]$weight == 1))
  expect_true(all(matches$weight %in% m.out$weights))
  expect_true(all(rownames(matches) %in% nms_matched))
  expect_true(all(nms_matched %in% rownames(matches)))
  
  # tests -- newdata
  expect_equal(sum(is.na(matches2)), 0L)
  expect_true(all(names(matches2) %in% c(names(lalonde2), "weight")))
  expect_true(all(c(names(lalonde2), "weight") %in% names(matches2)))
  expect_equal(ncol(matches2), ncol(lalonde2) + 1)
  expect_true(all(matches2[matches2$treat == 1,]$weight == 1))
  expect_true(all(matches2$weight %in% m.out$weights))
})

test_that("can get correct matches, genetic", {  ## this one is problematic -- problematic returns
  m.out <- matchit(treat ~ age + educ + black + hispan + married + nodegree +
                     re74 + re75, data = lalonde, method = "genetic")
  n_matched <- sum(m.out$nn[2,])
  nms_matched <- c(rownames(m.out$match.matrix), m.out$match.matrix[,1], 
                   m.out$match.matrix[which(!is.na(m.out$match.matrix[,2])), 2])
  
  lalonde$id <- 1:nrow(lalonde)
  lalonde2 <- lalonde[sample.int(500, 3000, replace=TRUE),]
  matches <- get_matches(m.out, lalonde)
  matches2 <- get_matches(m.out, model_frame= lalonde, id_cols= "id", newdata= lalonde2)
  
  # tests -- no newdata
  expect_equal(n_matched, nrow(matches))
  expect_equal(sum(matches$treat), m.out$nn[2,2])
  expect_equal(nrow(matches) - sum(matches$treat), m.out$nn[2,1])
  expect_equal(sum(is.na(matches)), 0L)
  expect_equal(names(matches), c(names(lalonde), "weight"))
  expect_equal(ncol(matches), ncol(lalonde) + 1)
  expect_true(all(matches[matches$treat == 1,]$weight == 1))
  expect_true(all(rownames(matches) %in% nms_matched))
  expect_true(all(nms_matched %in% rownames(matches)))
  
  # tests -- newdata
  expect_equal(sum(is.na(matches2)), 0L)
  expect_true(all(names(matches2) %in% c(names(lalonde2), "weight")))
  expect_true(all(c(names(lalonde2), "weight") %in% names(matches2)))
  expect_equal(ncol(matches2), ncol(lalonde2) + 1)
  expect_true(all(matches2[matches2$treat == 1,]$weight == 1))
})

test_that("can get correct matches, cem", {
  library(cem)
  m.out <- matchit(treat ~ age + educ + black + hispan + married + nodegree
                   + re74 + re75, data = lalonde, method = "cem")
  n_matched <- sum(m.out$nn[2,])
  nms_matched <- names(m.out$subclass[!is.na(m.out$subclass)])
  
  lalonde$id <- 1:nrow(lalonde)
  lalonde2 <- lalonde[sample.int(500, 3000, replace=TRUE),]
  matches <- get_matches(m.out, lalonde)
  matches2 <- get_matches(m.out, model_frame= lalonde, id_cols= "id", newdata= lalonde2)
  
  # tests -- no newdata
  expect_equal(n_matched, nrow(matches))
  expect_equal(sum(matches$treat), m.out$nn[2,2])
  expect_equal(nrow(matches) - sum(matches$treat), m.out$nn[2,1])
  expect_equal(sum(is.na(matches)), 0L)
  expect_equal(names(matches), c(names(lalonde), "weight"))
  expect_equal(ncol(matches), ncol(lalonde) + 1)
  expect_true(all(matches[matches$treat == 1,]$weight == 1))
  expect_true(all(matches$weight %in% m.out$weights))
  expect_true(all(rownames(matches) %in% nms_matched))
  expect_true(all(nms_matched %in% rownames(matches)))
  
  # tests -- newdata
  expect_equal(sum(is.na(matches2)), 0L)
  expect_true(all(names(matches2) %in% c(names(lalonde2), "weight")))
  expect_true(all(c(names(lalonde2), "weight") %in% names(matches2)))
  expect_equal(ncol(matches2), ncol(lalonde2) + 1)
  expect_true(all(matches2[matches2$treat == 1,]$weight == 1))
  expect_true(all(matches2$weight %in% m.out$weights))
  
})

l_treat <- do.call("rbind",replicate(10, lalonde[lalonde$treat == 1, ], simplify = FALSE))
l_contr <- do.call("rbind",replicate(3, lalonde[lalonde$treat == 0, ], simplify = FALSE))
l2 <- do.call("rbind", list(l_treat, l_contr)); rm(l_treat, l_contr)
rownames(l2) <- 1:nrow(l2)

test_that("correct weights for replace= TRUE", {
  m.out <- matchit(treat ~ re74 + re75 + educ + black + hispan + age,
                   data = l2, method = "nearest", replace= TRUE)
  n_matched <- sum(m.out$nn[2,])
  nms_matched <- c(rownames(m.out$match.matrix), c(m.out$match.matrix))
  
  matches <- get_matches(m.out, l2)
  exp_wts <- c(1, unique(as.vector(table(m.out$match.matrix))))
  
  # tests -- no newdata
  expect_equal(n_matched, nrow(matches))
  expect_equal(sum(matches$treat), m.out$nn[2,2])
  expect_equal(nrow(matches) - sum(matches$treat), m.out$nn[2,1])
  expect_equal(sum(is.na(matches)), 0L)
  expect_equal(names(matches), c(names(lalonde), "weight"))
  expect_equal(ncol(matches), ncol(lalonde) + 1)
  expect_true(all(matches[matches$treat == 1,]$weight == 1))
  expect_true(all(matches$weight %in% exp_wts))
  expect_true(all(rownames(matches) %in% nms_matched))
  expect_true(all(nms_matched %in% rownames(matches)))
})


test_that("correct weights for replace= TRUE; test2", {
  m.out <- matchit(treat ~ re74 + re75 + educ + black + hispan + age,
                   data = l2, method = "cem")
  n_matched <- sum(m.out$nn[2,])
  nms_matched <- names(m.out$subclass[!is.na(m.out$subclass)])
  
  matches <- get_matches(m.out, l2)
  exp_wts <- c(1, unique(as.vector(table(m.out$match.matrix))))
  
  # tests -- no newdata
  expect_equal(n_matched, nrow(matches))
  expect_equal(sum(matches$treat), m.out$nn[2,2])
  expect_equal(nrow(matches) - sum(matches$treat), m.out$nn[2,1])
  expect_equal(sum(is.na(matches)), 0L)
  expect_equal(names(matches), c(names(lalonde), "weight"))
  expect_equal(ncol(matches), ncol(lalonde) + 1)
  expect_true(all(matches[matches$treat == 1,]$weight == 1))
  expect_true(all(matches$weight %in% m.out$weights))
  expect_true(all(rownames(matches) %in% nms_matched))
  expect_true(all(nms_matched %in% rownames(matches)))
})
