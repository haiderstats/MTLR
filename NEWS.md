# MTLR 0.2.1
* Added the option to obtain survival probabilities at specified times (see `predict` with `type = "prob_times"`).
* When using `predict` with  "prob_event" there is no longer a linear extension added if the survival curve never reaches 0. The linear extension is still present in the "mean_time" and "median_time" options however.
* Internally some `sapply` commands were refactored so to not iterate over indexes.
* Tests were changed as to not clash with future ggplot2 updates -- tests were checking against internal structures as opposed to visual differences.

# MTLR 0.2.0
* Added an option to optimize concordance instead of the log-likelihood when selecting a value for C1 in `mtlr_cv`.
* Removed ellipses from `create_folds` since it's easy to use `type` instead of `foldtype` when passing in arguments.
* Fixed a bug where if all data was left or right censored an error was thrown since there was no interval censored data.
* `mtlr_cv` output is not `best_C1` and not `bestC1` to match the documentation.
* Removed usage of Rcpp data structures inside MTLRHelper.cpp. This was not needed and just wasting time when transferring them to armadillo structures.
* A number of printing fixes were made, e.g. printing out all variable names at once instead of one at a 
time for warnings and actually printing results when verbose = T for `mtlr_cv`.
* When biases are trained they now assume the data are uncensored. 

# MTLR 0.1.0
* First CRAN version.
