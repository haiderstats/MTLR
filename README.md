
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MTLR

[![R-CMD-check](https://github.com/haiderstats/MTLR/workflows/R-CMD-check/badge.svg)](https://github.com/haiderstats/MTLR/actions)
[![codecov](https://codecov.io/gh/haiderstats/MTLR/branch/master/graphs/badge.svg)](https://codecov.io/gh/haiderstats/MTLR)
![](http://cranlogs.r-pkg.org/badges/grand-total/MTLR)

The goal of `MTLR` is to provide an R implementation for [Multi-Task
Logistic
Regression](https://papers.nips.cc/paper/4210-learning-patient-specific-cancer-survival-distributions-as-a-sequence-of-dependent-regressors).
In addition to supplying the model provided by [Yu et
al.](https://papers.nips.cc/paper/4210-learning-patient-specific-cancer-survival-distributions-as-a-sequence-of-dependent-regressors)
we have extended the model for left censoring, interval censoring, and a
mixture of censoring types. Functionality includes training an MTLR
model, predicting survival curves for new observations, and plotting
these survival curves and feature weights estimated by MTLR.

## Installation

You can install the version from CRAN or the development version from
GitHub:

``` r
# CRAN:
install.packages("MTLR")

# GitHub:
# install.packages("devtools")
devtools::install_github("haiderstats/MTLR")
```

## Example

Given a survival dataset containing event time and event status
indicator (censored/uncensored), we can produce an MTLR model. For
example, consider the `lung` dataset from the `survival` package:

``` r
# Load survival for the lung dataset and the Surv() function.
library(survival)
#Here we will use 9 intervals (10 time points) just for plotting purposes. 
#The default is sqrt the number of observations.
mod <- mtlr(Surv(time,status)~., data = lung, nintervals = 9)
print(mod)
#> 
#> Call:  mtlr(formula = Surv(time, status) ~ ., data = lung, nintervals = 9) 
#> 
#> Time points:
#>  [1]  62.3 145.4 179.3 210.4 241.4 284.5 308.2 383.5 476.3 642.8
#> 
#> 
#> Weights:
#>           Bias      inst      age      sex ph.ecog ph.karno pat.karno meal.cal  wt.loss
#> 62.27   0.1113 -0.017979  0.04894 -0.01250 0.00465 -0.00652   0.01347 -0.02619 -0.01326
#> 145.36  0.1381 -0.021136  0.03275 -0.00471 0.02310 -0.02226  -0.01751 -0.01105 -0.02975
#> 179.27  0.2095 -0.008194  0.02260 -0.02562 0.02393 -0.02045  -0.03160 -0.02310 -0.02218
#> 210.36  0.0464  0.000371  0.00815 -0.03634 0.04702 -0.02227  -0.04127 -0.01410 -0.03367
#> 241.36 -0.1012  0.009570 -0.01580 -0.04404 0.06693 -0.04142  -0.05454 -0.00808 -0.01286
#> 284.55 -0.2372  0.004866 -0.00474 -0.05181 0.04829 -0.02112  -0.03138  0.00237 -0.02549
#> 308.18 -0.2960 -0.007731 -0.00466 -0.05426 0.04912 -0.02284  -0.03266 -0.01607  0.02130
#> 383.45 -0.0286 -0.019974 -0.01029 -0.03263 0.02870 -0.00682  -0.02323 -0.01458 -0.00236
#> 476.27 -0.1277 -0.010224  0.00106 -0.02226 0.01716  0.01699  -0.01556 -0.02111 -0.01589
#> 642.82 -0.3538 -0.014985  0.02291 -0.02550 0.01846  0.00557  -0.00854  0.00376 -0.00555
#Plot feature weights:
plot(mod)
```

<img src="man/figures/README-example-1.png" style="display: block; margin: auto;" />

``` r
#Get survival curves for the lung dataset:
curves <- predict(mod)
#Plot the first 20 survival curves:
plotcurves(curves, 1:20)
```

<img src="man/figures/README-example-2.png" style="display: block; margin: auto;" />
