#' MTLR Internal Cross-Validation for Selecting C1.
#'
#' @inheritParams mtlr
#' @param C1_vec a vector of regularization parameters to test. All values must be non-negative. For large datasets you may want to reduce the number
#' of value tried to increase efficiency. Similarly for nfolds.
#' @param loss a string indicating the loss to optimize for which to choose the regularization parameter. Currently only the log-likelihood (ll)
#'  is supported.
#' @param nfolds the number of internal cross validation folds, default is 5.
#' @param foldtype type of cross validation folds. Full stratification, "fullstrat", sorts observations by their event time and their event indicators
#' and numbers them off into folds. This effectively give each fold approximately the same number of uncensored observations as well as keeps the range
#' of time points as equivalent as possible across folds. This type of cross-validation is completely determininstic.
#' Censored stratification, "censorstrat", will put approximately the same number of uncensored observations in each fold but not pay any attention to
#' event time. This is partially stochastic. The totally random cross-validation, "random", randomly assigns observations to folds without considering
#' event time nor event status.
#' @details Currently only the log-likelihood loss is supported for optimizing C1. Here the loss considers censored and uncensored observations differently.
#' For uncensored observations, we assign a loss of the negative log probability assigned to the interval in which the observation had their event, \emph{e.g.}
#' if an observation had a 20% chance of having its event between timepoint1 and timepoint2 and it did have it's event in that interval then the loss
#' is -log(0.2). We want these probabilities to be large so we would normally want to maximize this value (since logs of probabilities are negative)
#' but we take the negative and instead minimize the value, thus we want the lowest loss. For censored observations we take the log of the probability
#' of survival at the time of censoring, \emph{e.g.} if an observation is censored at time = 42 we take the negative log of the survival probability assigned
#' to time 42 as the loss.
#' @return Performing mtlr_cv will return the following:
#' \itemize{
#'   \item best_C1: The value of C1 which achieved the best (lowest) loss.
#'   \item avg_loss: The averaged value of loss across the five folds for each value of C1 tested.
#'   }
#' @examples
#' library(survival)
#' cv_mod = mtlr_cv(Surv(time,status)~., data = lung)
#' #Note the best C1 also corresponds to the lost average loss:
#' cv_mod
#' @seealso \code{\link[MTLR]{mtlr}}
#' @export
mtlr_cv <- function(formula,
                 data,
                 nintervals = NULL,
                 normalize = T,
                 C1_vec = c(0.001,0.01,0.1,1,10,100,1000),
                 train_biases = T,
                 loss = c("ll"),
                 nfolds = 5,
                 foldtype = c("fullstrat","censorstrat","random"),
                 threshold = 1e-05,
                 maxit = 5000,
                 lower = -20,
                 upper = 20){
  if(any(C1_vec < 0)){
    stop("All values of C1 must be non-negative.")
  }
  foldtype = match.arg(foldtype)


  mf <- stats::model.frame(formula = formula, data)
  y <- stats::model.response(mf)
  time <- y[,1]
  delta <- y[,2]
  fold_index <- foldme(time,delta,nfolds,foldtype)

  data <- data[stats::complete.cases(data),]
  res_mat <- matrix(rep(0,nfolds*length(C1_vec)), ncol = length(C1_vec),nrow =nfolds)
  for(fold in 1:nfolds){
    datacv <- data[-fold_index[[fold]],]
    result <- c()
    for(i in C1_vec){
      mod <- mtlr(formula, datacv,nintervals,normalize,i, train_biases, threshold, maxit, lower, upper)
      result <- c(result, loglik_loss(mod,data[fold_index[[fold]],]))
    }
    res_mat[fold,] <- result
  }
  avg_results <- apply(res_mat, 2, mean)
  names(avg_results) = C1_vec
  best_C1 <- C1_vec[which.min(avg_results)]


  to_return = list(bestC1 = best_C1, avg_loss = avg_results)
  return(to_return)
}







