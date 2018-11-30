testthat::context("Testing the primary mtlr function.")



# mtlr ------------------------------------------------------------
testthat::test_that("mtlr function is consistent for basic survival dataset",{
  formula = survival::Surv(time,status)~.
  data = survival::leukemia
  expect_equal_to_reference(mtlr(formula,data),"mtlr_leuk.rds")
})

testthat::test_that("mtlr function is consistent for more complex survival dataset",{
  formula = survival::Surv(time,status)~.
  data = survival::lung
  expect_equal_to_reference(mtlr(formula,data),"mtlr_lung.rds")
})

testthat::test_that("mtlr function is consistent for more complex survival dataset - no extra bias training",{
  formula = survival::Surv(time,status)~.
  data = survival::lung
  expect_equal_to_reference(mtlr(formula,data, train_biases = F),"mtlr_lung_nobias.rds")
})

testthat::test_that("mtlr function is consistent for all censored survival dataset",{
  formula = survival::Surv(time,status)~.
  data = survival::leukemia
  data = data[data$status == 0,]
  expect_equal_to_reference(mtlr(formula,data),"mtlr_censored.rds")
})

testthat::test_that("mtlr function is consistent for all uncensored survival dataset",{
  formula = survival::Surv(time,status)~.
  data = survival::leukemia
  data = data[data$status == 1,]
  expect_equal_to_reference(mtlr(formula,data),"mtlr_uncensored.rds")
})

testthat::test_that("mtlr function is consistent for basic survival dataset UNNORMALIZED",{
  formula = survival::Surv(time,status)~.
  data = survival::leukemia
  expect_equal_to_reference(mtlr(formula,data, normalize = F),"mtlr_leukUNNORMALIZED.rds")
})

testthat::test_that("mtlr function is consistent for basic survival dataset for chosen nintervals",{
  formula = survival::Surv(time,status)~.
  data = survival::leukemia
  expect_equal_to_reference(mtlr(formula,data, nintervals = 3),"mtlr_leuk_timepoints.rds")
})



testthat::test_that("mtlr argument specifications are working.",{
  formula = survival::Surv(time,status)~.
  data = survival::leukemia
  data$time[1] = -10
  expect_error(mtlr(formula,data),"All event times must be non-negative")
  data$time[1] = 10

  formula = time~.
  expect_error(mtlr(formula,data),"The response must be a Surv object.")

  formula = survival::Surv(time,status)~.

  expect_error(mtlr(formula,data, C1 = -10),"C1 must be non-negative.")
  expect_error(mtlr(formula,data, C1 = -1e-10),"C1 must be non-negative.")
  expect_error(mtlr(formula,data, threshold = -1e-10),"The threshold must be positive.")
  expect_error(mtlr(formula,data, threshold = 0),"The threshold must be positive.")
})

testthat::test_that("when training mtlr fails optim error is caught",{
  formula = survival::Surv(time,status)~.
  data = survival::lung
  data$meal.cal = data$meal.cal*1e100
  expect_error(mtlr(formula,data,normalize = F))

})















