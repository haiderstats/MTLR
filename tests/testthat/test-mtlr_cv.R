testthat::context("Testing the mtlr_cv function. Most tests are just checking against stored values.")



# mtlr_cv ------------------------------------------------------------
testthat::test_that("mtlr_cv function is consistent for basic survival dataset",{
  formula = survival::Surv(time,status)~.
  data = survival::leukemia
  expect_equal_to_reference(mtlr_cv(formula,data),"mtlrcv_leuk.rds")
})

testthat::test_that("mtlr_cv function is consistent for more complex survival dataset",{
  formula = survival::Surv(time,status)~.
  data = survival::lung
  expect_equal_to_reference(mtlr_cv(formula,data),"mtlrcv_lung.rds", tolerance = 1e-3)
})

testthat::test_that("mtlr_cv function is consistent for all censored survival dataset",{
  formula = survival::Surv(time,status)~.
  data = survival::lung
  data = data[data$status == 1,]
  expect_equal_to_reference(mtlr_cv(formula,data),"mtlrcv_censored.rds", tolerance = 1e-3)
})

testthat::test_that("mtlr_cv function is consistent for all uncensored survival dataset",{
  formula = survival::Surv(time,status)~.
  data = survival::leukemia
  data = data[data$status == 1,]
  expect_equal_to_reference(mtlr_cv(formula,data),"mtlrcv_uncensored.rds")
})

testthat::test_that("mtlr_cv function catches negative C1 survival dataset",{
  formula = survival::Surv(time,status)~.
  data = survival::leukemia
  data = data[data$status == 1,]
  C1_vec = c(-1,1)
  expect_error(mtlr_cv(formula,data,C1_vec = C1_vec),"All values of C1 must be non-negative.")
})


testthat::test_that("mtlr_cv function works with multiple types of censoring",{
  time1 = c(NA, 4, 7, 12, 10, 6, NA, 3,5,9,10,12,NA,4,6,2,NA,16,15,11)
  time2 = c(14, 4, 10, 12, NA, 9, 5, NA, NA, NA, NA, 15,22,4,8,6,2,20,23,11)
  set.seed(42)
  dat = cbind.data.frame(time1, time2, importantfeature1 = rnorm(20),importantfeature2 = rnorm(20),
                         importantfeature3 = rnorm(20),importantfeature4 = rnorm(20),importantfeature5 = rnorm(20),
                         importantfeature6 = rbinom(20,1,.3),importantfeature7 = rbinom(20,1,.3))
  formula = survival::Surv(time1,time2,type = "interval2")~.
  expect_equal_to_reference(mtlr_cv(formula, dat),"mtlrcv_mixed_censoring.rds",tolerance = 1e-3)
})



