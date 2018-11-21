#' Primary MTLR function. Run the mode
#'
#' @inheritParams mtlr
#' @param C1 a vector of regularization parameters to test. All values must be non-negative.
#' @param loss a string indicating the loss to optimize for which to choose the regularization parameter. Currently on the log-likelihood (ll) is supported.
#' @export
mtlr_cv <- function(formula,
                 data,
                 nintervals = NULL,
                 normalize = T,
                 C1 = c(0.01,0.1,1,10,100),
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
  time = y[,1]
  delta = y[,2]
  fold_index = foldme(time,delta,foldtype,nfolds)

  data = data[complete.cases(data),]
  resultsMatrix = matrix(rep(0,nfolds*length(C1)), ncol = length(C1),nrow =nfolds)
  for(fold in 1:nfolds){
    datacv = data[-fold_index[[fold]],]
    result = c()
    for(i in C1){
      mod = mtlr(formula, datacv,ninternvals,normalize,i, train_biases, threshold, maxit, lower, upper)
      result = c(result, ll(mod,datacv[fold_index[[fold]],]))
    }
  }
}







