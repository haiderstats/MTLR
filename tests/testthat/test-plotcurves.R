testthat::context("plotcurves is consistent with previous plots.")



# plotcurves ------------------------------------------------------------
testthat::test_that("plotcurves function is consistent",{
  formula <- survival::Surv(time,status)~.
  data <- survival::lung
  mod <- mtlr(formula,data)
  curves <- predict(mod)
  expect_warning(plotcurves(curves, index=  1:3, color = c("red","green")))
  oneCurve <- plotcurves(curves)
  expect_doppelganger("oneCurve", oneCurve)
  onePinkCurve <- plotcurves(curves, color = "pink")
  expect_doppelganger("onePinkCurve", onePinkCurve)
  manyCurves <- plotcurves(curves, 1:10)
  expect_doppelganger("manyCurves", manyCurves)
  shortCurve <- plotcurves(curves, xlim = c(0,4))
  expect_doppelganger("shortCurve", shortCurve)
})
