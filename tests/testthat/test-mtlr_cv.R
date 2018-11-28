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
  expect_equal_to_reference(mtlr_cv(formula,data),"mtlrcv_lung.rds")
})

testthat::test_that("mtlr_cv function is consistent for all censored survival dataset",{
  formula = survival::Surv(time,status)~.
  data = survival::lung
  data = data[data$status == 1,]
  expect_equal_to_reference(mtlr_cv(formula,data),"mtlrcv_censored.rds")
})

testthat::test_that("mtlr_cv function is consistent for all uncensored survival dataset",{
  formula = survival::Surv(time,status)~.
  data = survival::leukemia
  data = data[data$status == 1,]
  expect_equal_to_reference(mtlr_cv(formula,data),"mtlrcv_uncensored.rds")
})

testthat::test_that("mtlr_cv function is consistent for all uncensored survival dataset",{
  formula = survival::Surv(time,status)~.
  data = survival::leukemia
  data = data[data$status == 1,]
  C1_vec = c(-1,1)
  expect_error(mtlr_cv(formula,data,C1_vec = C1_vec),"All values of C1 must be non-negative.")
})





