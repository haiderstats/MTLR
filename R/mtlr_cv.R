#' MTLR Internal Cross-Validation for Selecting C1.
#'
#' @inheritParams mtlr
#' @param C1_vec a vector of regularization parameters to test. All values must be non-negative. For large datasets you may want to reduce the number
#' of value tried to increase efficiency. Similarly for nfolds.
#' @param previous_weights a boolean specifying if sequential folds should use the previous fold's parameters as seed_weights. Doing this will likely
#' speed up the computation time for cross-validation as we are providing weights which are (likely) close to the optimal weights. Note that this is
#' done separately for each value of C1 so there is no parameter sharing between different values of C1, and instead only across the same value of C1.
#' @param loss a string indicating the loss to optimize for which to choose the regularization parameter. Currently one can optimize for the log-likelihood ("ll")
#'  or concordance ("concordance"). See details regarding these losses.
#' @param nfolds the number of internal cross validation folds, default is 5.
#' @param foldtype type of cross validation folds. Full stratification, "fullstrat", sorts observations by their event time and their event indicators
#' and numbers them off into folds. This effectively give each fold approximately the same number of uncensored observations as well as keeps the range
#' of time points as equivalent as possible across folds. This type of cross-validation is completely deterministic.
#' Censored stratification, "censorstrat", will put approximately the same number of uncensored observations in each fold but not pay any attention to
#' event time. This is partially stochastic. The totally random cross-validation, "random", randomly assigns observations to folds without considering
#' event time nor event status.
#' @param verbose if TRUE the progress will be printed for every completed value of C1.
#' @details The log-likelihood loss and concordance are supported for optimizing C1. Here the log-likelihood loss considers censored and uncensored observations differently.
#' For uncensored observations, we assign a loss of the negative log probability assigned to the interval in which the observation had their event, \emph{e.g.}
#' if an observation had a 20% chance of having its event between timepoint1 and timepoint2 and it did have it's event in that interval then the loss
#' is -log(0.2). We want these probabilities to be large so we would normally want to maximize this value (since logs of probabilities are negative)
#' but we take the negative and instead minimize the value, thus we want the lowest loss. For censored observations we take the log of the probability
#' of survival at the time of censoring, \emph{e.g.} if an observation is censored at time = 42 we take the negative log of the survival probability assigned
#' to time 42 as the loss.
#'
#' For the concordance loss, C1 is chosen to maximize the overall concordance when using the negative median as the "risk" score. This is completed using survConcordance in the survival package.
#' @return Performing mtlr_cv will return the following:
#' \itemize{
#'   \item best_C1: The value of C1 which achieved the best (lowest) loss.
#'   \item avg_loss: The averaged value of loss across the five folds for each value of C1 tested.
#'   }
#' @examples
#' library(survival)
#' cv_mod <- mtlr_cv(Surv(time,status)~., data = lung)
#' #Note the best C1 also corresponds to the lost average loss:
#' cv_mod
#' @seealso \code{\link[MTLR]{mtlr}}
#' @export
mtlr_cv <- function(formula,
                 data,
                 time_points = NULL,
                 nintervals = NULL,
                 normalize = T,
                 C1_vec = c(0.001,0.01,0.1,1,10,100,1000),
                 train_biases = T,
                 train_uncensored = T,
                 seed_weights = NULL,
                 previous_weights = T,
                 loss = c("ll","concordance"),
                 nfolds = 5,
                 foldtype = c("fullstrat","censorstrat","random"),
                 verbose = FALSE,
                 threshold = 1e-05,
                 maxit = 5000,
                 lower = -15,
                 upper = 15
                 ){
  if(any(C1_vec < 0)){
    stop("All values of C1 must be non-negative.")
  }
  if(any(dim(data)==0)){
    stop("Dimensions of the dataset must be non-zero.")
  }
  foldtype <- match.arg(foldtype)
  loss <- match.arg(loss)
  lossfnc <- switch(loss,
                    ll = loglik_loss,
                    concordance = concordance_loss
  )
  mf <- stats::model.frame(formula = formula, data)
  y <- stats::model.response(mf)
  time <- y[,1]
  delta <- y[,2]
  fold_index <- create_folds(time,delta,nfolds,foldtype)
  time_delta_names = attr(y,"dimnames")[[2]]
  #There may be some NAs in the time if we are in the inteval setting. We do not want to remove these rows since they
  #are valid NAs. However, we only want x values with nonmissing data.
  data <- data[stats::complete.cases(data[which(!names(data)%in%time_delta_names)]),]
  res_mat <- matrix(rep(0,nfolds*length(C1_vec)), ncol = length(C1_vec),nrow =nfolds)
  #Give the same number of intervals to every fold.
  if(is.null(time_points) & is.null(nintervals)){
    nintervals <- ceiling(sqrt(((nfolds-1)/nfolds)*nrow(y))) #Default number of intervals is the sq. root of the number of observations in k-1/k % of the data.
  }
  if(previous_weights){
    parList <- list()
  }
  for(fold in 1:nfolds){
    if(verbose){
      print(paste("Starting fold ",fold," of ",nfolds,".", sep = ""))
    }
    datacv <- data[-fold_index[[fold]],]
    result <- c()
    for(i in seq_along(C1_vec)){
      if(previous_weights){
        if(fold == 1){
          mod <- mtlr(formula,datacv,time_points,nintervals,
                      normalize, C1_vec[i], train_biases,train_uncensored, seed_weights,
                      threshold, maxit, lower, upper)
          parList[[i]] <- c(mod$weight_matrix)
        }else{
          mod <- mtlr(formula,datacv,time_points,nintervals,
                      normalize, C1_vec[i], train_biases,train_uncensored, parList[[i]],
                      threshold, maxit, lower, upper)
        }
      }else{
        mod <- mtlr(formula,datacv,time_points,nintervals,
                    normalize, C1_vec[i], train_biases,train_uncensored, seed_weights,
                    threshold, maxit, lower, upper)
      }
      result <- c(result, lossfnc(mod,data[fold_index[[fold]],]))
    }
    res_mat[fold,] <- result
  }
  avg_results <- apply(res_mat, 2, mean)
  names(avg_results) <- C1_vec
  best_C1 <- C1_vec[which.min(avg_results)]
  if(loss == "concordance"){
    avg_results <- -1*avg_results
  }
  to_return <- list(best_C1 = best_C1, avg_loss = avg_results)
  return(to_return)
}







