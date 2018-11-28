#' @useDynLib MTLR
#' @importFrom Rcpp sourceCpp
NULL

#' Train a Mutli-Task Logistic Regression (MTLR) Model
#'
#' Trains a MTLR model for survival prediction. Right, left, and interval censored data is all supported.
#'
#' @param formula a formula object with the response to the left of the "~" operator. The response must be a survival object returned
#' by the \code{\link[survival]{Surv}} function.
#' @param data a data.frame containing the features for survival prediction. These must be variables corresponding to the formula object.
#' @param nintervals Number of time intervals to use for MTLR. Note the number of time points will be nintervals + 1. If left as NULL
#' a default of sqrt(N) is used where N is the number of observations in the supplied dataset.
#' @param normalize if TRUE, variables will be normalized (mean 0, standard deviation of 1). This is STRONGLY suggested. If normalization
#' does not occur it is much more likely that MTLR will fail to converge. Additionally, if FALSE consider adjusting "lower" and "upper"
#' used for L-BFGS-B optimization.
#' @param C1 The L2 regularization parameter for MTLR. C1 can also be selected via \code{\link[MTLR]{mtlr_cv}}. See "Learning Patient-Specific Cancer Survival Distributions as a Sequence of Dependent
#'  Regressors" by Yu et al. (2011) for details.
#' @param train_biases if TRUE, biases will be trained before feature weights (and again trained while training feature weights). This
#' has shown to speed up total training time.
#' @param threshold The threshold for training the weights. This threshold will be passed to \link[stats]{optim}.
#' @param maxit The maximum interations to run for MTLR. This parameter will be passed to \link[stats]{optim}.
#' @param lower The lower bound for L-BFGS-B optimization. This parameter will be passed to \link[stats]{optim}.
#' @param upper The upper bound for L-BFGS-B optimization. This parameter will be passed to \link[stats]{optim}.
#' @details This function allows one to train an MTLR model given a dataset containing survival data. mtlr uses the Limited-Memory
#' Broyden–Fletcher–Goldfarb–Shanno (L-BFGS-B) approximation method to train feature weights. This training is outsourced to the internal
#' \link[stats]{optim} function in R. Currently only a few parameters (namely threshold, maxit,lower, upper) of optim are supported, more will
#' likely become avaliable in the future.
#'
#' Weights are initialized to 0 prior to training. Under default settings, the bias weights
#' will be trained before considering feature weights. As Yu et al. (2011) specified, the introduction of censored observations creates a non-convex
#' loss function. To address this, weights are first trained assuming all patietns event times are \emph{uncensored}. Once these starting weights have
#' been trained another round of training is performed using the true values of the event indicator (censored/uncensored). Future iterations of this
#' package will add more options to address this non-convexity.
#'
#' Yu et al. (2011) actually suggested two regularization parameters, C1 to control the size of the feature weights and C2 to control the smoothness.
#' In Ping Jin's masters thesis (Using Survival Prediction Techniques to Learn Consumer-Specific Reservation Price Distributions) he showed that C2
#' is not required for smoothness and C1 will suffice (Appendix A.2) so we do not support the C2 parameter in this implemenetation.
#'
#' @return An mtlr object returns the following:
#' \itemize{
#'   \item weight_matrix: The matrix of feature weights determined by MTLR.
#'   \item x: The dataframe of features (response removed). Note observations with missing values will have been removed (this is the dataset on which
#'   MTLR was trained).
#'   \item y: The matrix of response values MTLR uses for training. Each column corresponds to an observation and rows as time points. A value of 1
#'   indicates a observation was either censored or had their event occur by that time.
#'   \item response: The response as a Surv object (specified by formula).
#'   \item time_points: The timepoints selected and used to train MTLR.
#'   \item C1: The regularization parameter used.
#'   \item Call: The original call to mtlr.
#'   \item Terms: The x-value terms used in mtlr. These are later used in \code{\link[MTLR]{predict.mtlr}}
#'   \item scale: The means and standard deviations of features when normalize = TRUE. These are used in \code{\link[MTLR]{predict.mtlr}}. Will be
#'   NULL if normalize = FALSE.
#'   \item xlevels: The levels of the features used. This is used again by \code{\link[MTLR]{predict.mtlr}}.
#' }
#' @seealso
#' \code{\link[MTLR]{predict.mtlr}} \code{\link[MTLR]{plot.mtlr}} \code{\link[MTLR]{mtlr_cv}} \code{\link[MTLR]{predict.mtlr}}
#' @export
mtlr <- function(formula,
                 data,
                 nintervals = NULL,
                 normalize = T,
                 C1 = 1,
                 train_biases = T,
                 threshold = 1e-05,
                 maxit = 5000,
                 lower = -20,
                 upper = 20){
  cl <- match.call() #Save a copy of the function call.

  #Data setup
  mf <- stats::model.frame(formula = formula, data)
  Terms <- attr(mf, "terms")
  xlevels <- stats::.getXlevels(Terms, mf)

  x <- stats::model.matrix(Terms, data = mf)
  x <- x[,-1,drop=FALSE] #Remove intercept term -- We will handle biases later.
  y <- stats::model.response(mf)

  if (!survival::is.Surv(y))
    stop("The response must be a Surv object.")
  if(attr(y, "type") != "right")
    stop("Currently only right censored data is supported.")

  time <- y[,1]
  delta <- y[,2] #Death indicator (censored = 0, death = 1)

  if(any(time < 0))
    stop("All event times must be positive.")
  if(C1 < 0)
    stop("C1 must be non-negative.")
  if(threshold <= 0)
    stop("The threshold must be positive.")

  if(is.null(nintervals))
    nintervals <- sqrt(ceiling(nrow(y))) #Default number of intervals is the sq. root of the number of observations.
  if(normalize){
    x <- scale(x)
    scale_centers <- attr(x, "scaled:center")
    scale_scale <- attr(x,"scaled:scale")
    scales = list(center= scale_centers, sd = scale_scale)
  }else{
    scales = NULL
  }

  ##############################################################
  #Prepare data for MTLR Rcpp functions.
  #The Rcpp functions are built to work on data that comes in with observations ordered by their event status
  #(censored observations first).
  #Here we order our data by the censor status.
  ord = order(delta)

  if(ncol(x)){
    x <- x[ord,1:ncol(x),drop=FALSE]
  }

  m <- nintervals + 1   #The number of time points to evaulate will be the number of time intervals + 1.
  quantiles <- seq(0,1,length.out = m+2)[-c(1,m+2)] #We will select time point based on the distribution of the times.
  time_points <- unname(stats::quantile(time, quantiles))


  time_points <- time_points[!duplicated(time_points)]#Small data (or common event times) we will have duplicated time points -- we remove them.

  #We make a matrix where each column is a vector of indicators if an observation is dead at each time point, e.g. (0,0,,...,0,1,1,...1).
  #(See 'Learning Patient-Specific Cancer Survival Distributions as a Sequence of Dependent Regressors' Page 3.)
  y_matrix = matrix(1 - Reduce(c,Map(function(ind) time[ord] > time_points[ind], seq_along(time_points))), ncol = length(time), byrow = T)

  #We first train the biases (by setting feature parameters to zero (dAsZero)). Then we train the feature values as if all patients
  #were uncensored (this creates "good" starting values since with censored patients the objective is non-convex). Then we finally
  #train the data with their true censor statuses.
  threshold_factor = threshold/.Machine$double.eps

  if(train_biases){
    zero_matrix <- matrix(0,ncol = ncol(x), nrow = nrow(x))   #We create a zero_matrix to train the biases

    bias_par <- stats::optim(par = rep(0,length(time_points)*(ncol(x) +1)),fn = mtlr_objVal,gr = mtlr_grad, yval = y_matrix,
                     featureVal = zero_matrix, C1=C1, delta = sort(delta),
                     method = "L-BFGS-B", lower = lower, upper = upper, control = c(maxit = maxit, factr = threshold_factor))
    if(bias_par$convergence == 52)
      stop(paste("Error occured while training MTLR. Optim Error: ", bias_par$message))
  }else{
    bias_par <- list(par = rep(0,length(time_points)*(ncol(x) +1)))
  }

  params_uncensored <- stats::optim(par = bias_par$par,fn = mtlr_objVal, gr = mtlr_grad, yval = y_matrix, featureVal = x, C1 = C1, delta = rep(1,nrow(x)),
                       method = "L-BFGS-B", lower = lower, upper = upper, control = c(maxit = maxit, factr = threshold_factor))
  if(params_uncensored$convergence == 52)
    stop(paste("Error occured while training MTLR. Optim Error: ", params_uncensored$message))

  final_params <- stats::optim(par = params_uncensored$par,fn = mtlr_objVal,gr = mtlr_grad, yval = y_matrix, featureVal = x, C1 = C1,delta = sort(delta),
                    method = "L-BFGS-B", lower = lower, upper = upper, control = c(maxit = maxit, factr = threshold_factor))
  if(final_params$convergence == 52)
    stop(paste("Error occured while training MTLR. Optim Error: ", final_params$message))

  weights <- matrix(final_params$par, ncol = ncol(x) + 1,byrow=FALSE)
  colnames(weights) <- c("Bias",colnames(x))
  rownames(weights) <- round(time_points,2)
  x <- x[order(ord),0:ncol(x),drop = FALSE]
  y_matrix <- y_matrix[,order(ord)]
  fit <- list(weight_matrix = weights,
              x = x,
              y = y_matrix,
              response = y,
              time_points = time_points,
              C1 = C1,
              Call = cl,
              Terms = Terms,
              scale = scales,
              xlevels = xlevels)
  class(fit) <- "mtlr"
  fit
}







