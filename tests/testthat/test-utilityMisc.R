testthat::context("Testing the extra utility methods (get_param_influence,log_loss).")



# get_param_influence ------------------------------------------------------------
#Basic test for basic function (mostly to increase code coverage).
testthat::test_that("get_influence functionality",{
  formula = survival::Surv(time,status)~.
  data = survival::leukemia
  mod = mtlr(formula,data)
  weights = mod$weight_matrix
  expect_equal(unname(get_param_influence(mod)), sum(abs(weights[,2, drop=FALSE])))
})

testthat::test_that("log_loss functionality all uncensored",{
  formula <- survival::Surv(time,status)~.
  data <- survival::leukemia
  object <- mtlr(formula,data)

  #All uncensored
  newdata <- survival::leukemia[c(1,2,4,5),]
  curves <- predict(object, newdata)

  time <- curves[,1]
  #We only need one curve because everyone is the same.
  curve <- curves[,2]
  curve_diff = -diff(c(curve,0))
  ind = findInterval(newdata$time, time)
  loss = -sum(log(curve_diff[ind]+1e-05))
  expect_equal(log_loss(object,newdata),loss/nrow(newdata))
})

# log_loss ------------------------------------------------------------
testthat::test_that("log_loss functionality all censored",{
  formula <- survival::Surv(time,status)~.
  data <- survival::leukemia
  object <- mtlr(formula,data)

  #All uncensored
  newdata <- survival::leukemia[c(3,6,9,11),]
  curves <- predict(object, newdata)

  time <- curves[,1]
  #We only need one curve because everyone is the same.
  curve <- curves[,2]
  probs = predict_prob(curve,time,newdata$time)

  loss = -sum(log(probs+1e-05))
  expect_equal(log_loss(object,newdata),loss/nrow(newdata))
})

testthat::test_that("log_loss functionality mixed uncensored and censored",{
  formula <- survival::Surv(time,status)~.
  data <- survival::leukemia
  object <- mtlr(formula,data)

  #All uncensored
  newdata <- survival::leukemia[c(1,2,4,5,3,6,9,11),]
  curves <- predict(object, newdata)

  time <- curves[,1]
  #We only need one curve because everyone is the same.
  curve <- curves[,2]
  #Uncensored
  curve_diff = -diff(c(curve,0))
  ind = findInterval(newdata$time[1:4], time)
  loss = -sum(log(curve_diff[ind]+1e-05))

  #Censored
  probs = predict_prob(curve,time,newdata$time[5:8])
  loss = loss -sum(log(probs+1e-05))

  expect_equal(log_loss(object,newdata),loss/nrow(newdata))
})

















