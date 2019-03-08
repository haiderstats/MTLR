# MTLR 0.2.0
* Added an option to optimize concordance instead of the log-likelihood when selecting a value for C1 in `mtlr_cv`.
* Removed ellipses from `create_folds` since it's easy to use `type` instead of `foldtype` when passing in arguments.
* Fixed a bug where if all data was left or right censored an error was thrown since there was no interval censroed data.
* `mtlr_cv` output is not `best_C1` and not `bestC1` to match the documentation.
* Removed usage of Rcpp data structures inside MTLRHelper.cpp. This was not needed and just wasting time when transferring them to armadillo structures.
* A number of printing fixes were made, e.g. printing out all variable names at once instead of one at a 
time for warnings and actually printing results when verbose = T for `mtlr_cv`.
* When biases are trained they now assume the data are uncensored. 

# MTLR 0.1.0
* First CRAN version.
