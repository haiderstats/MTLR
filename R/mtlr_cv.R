#' Primary MTLR function. Run the mode
#'
#' @inheritParams mtlr
#' @param C1 a vector of regularization parameters to test. All values must be non-negative.
#' @param loss a string indicating the loss to optimize for which to choose the regularization parameter. Currently on the log-likelihood (ll) is supported.
#' @param nfolds the number of internal cross validation folds
#' @param foldtype type of cross validation folds
#' @export
mtlr_cv <- function(formula,
                 data,
                 nintervals = NULL,
                 normalize = T,
                 C1 = c(0.001,0.01,0.1,1,10,100,1000),
                 train_biases = T,
                 loss = c("ll"),
                 nfolds = 5,
                 foldtype = c("fullstrat","censorstrat","random"),
                 threshold = 1e-05,
                 maxit = 5000,
                 lower = -20,
                 upper = 20){
  if(any(C1 < 0)){
    stop("All values of C1 must be non-negative.")
  }
  foldtype = match.arg(foldtype)


  mf <- stats::model.frame(formula = formula, data)
  y <- stats::model.response(mf)
  time <- y[,1]
  delta <- y[,2]
  fold_index <- foldme(time,delta,nfolds,foldtype)

  data <- data[stats::complete.cases(data),]
  res_mat <- matrix(rep(0,nfolds*length(C1)), ncol = length(C1),nrow =nfolds)
  for(fold in 1:nfolds){
    datacv <- data[-fold_index[[fold]],]
    result <- c()
    for(i in C1){
      mod <- mtlr(formula, datacv,nintervals,normalize,i, train_biases, threshold, maxit, lower, upper)
      result <- c(result, log_loss(mod,data[fold_index[[fold]],]))
    }
    res_mat[fold,] <- result
  }
  avg_results <- apply(res_mat, 2, mean)
  names(avg_results) = C1
  best_C1 <- C1[which.min(avg_results)]


  to_return = list(bestC1 = best_C1, avg_loss = avg_results)
  return(to_return)
}







