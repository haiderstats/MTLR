testthat::context("plotcurves is consistent with previous plots.")



# plotcurves ------------------------------------------------------------
testthat::test_that("plotcurves function is consistent",{
  formula <- survival::Surv(time,status)~.
  data <- survival::lung
  mod <- mtlr(formula,data)
  curves <- predict(mod)
  expect_warning(plotcurves(curves, index=  1:3, color = c("red","green")))
  oneCurve <- plotcurves(curves)
  vdiffr::expect_doppelganger("oneCurve", oneCurve)
  onePinkCurve <- plotcurves(curves)
  vdiffr::expect_doppelganger("onePinkCurve", onePinkCurve)
  manyCurves <- plotcurves(curves)
  vdiffr::expect_doppelganger("manyCurves", manyCurves)
  shortCurve <- plotcurves(curves)
  vdiffr::expect_doppelganger("shortCurve", shortCurve)
})
