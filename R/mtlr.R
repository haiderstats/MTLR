#' @useDynLib MTLR
#' @importFrom Rcpp sourceCpp
NULL

#' Train a Multi-Task Logistic Regression (MTLR) Model
#'
#' Trains a MTLR model for survival prediction. Right, left, and interval censored data are all supported.
#'
#' @param formula a formula object with the response to the left of the "~" operator. The response must be a survival object returned
#' by the \code{\link[survival]{Surv}} function.
#' @param data a data.frame containing the features for survival prediction. These must be variables corresponding to the formula object.
#' @param time_points the time points for MTLR to create weights. If left as NULL, the time_points chosen will be based on equally spaced quantiles
#' of the survival times. In the case of interval censored data note that only the start time is considered and not the end time for selecting time points.
#' It is strongly recommended to specify time points if your data is heavily interval censored. If time_points is not NULL then nintervals is ignored.
#' @param nintervals Number of time intervals to use for MTLR. Note the number of time points will be nintervals + 1. If left as NULL
#' a default of sqrt(N) is used where N is the number of observations in the supplied dataset. This parameter is ignored if time_points is specified.
#' @param normalize if TRUE, variables will be normalized (mean 0, standard deviation of 1). This is STRONGLY suggested. If normalization
#' does not occur it is much more likely that MTLR will fail to converge. Additionally, if FALSE consider adjusting "lower" and "upper"
#' used for L-BFGS-B optimization.
#' @param C1 The L2 regularization parameter for MTLR. C1 can also be selected via \code{\link[MTLR]{mtlr_cv}}. See "Learning Patient-Specific Cancer Survival Distributions as a Sequence of Dependent
#'  Regressors" by Yu et al. (2011) for details.
#' @param train_biases if TRUE, biases will be trained before feature weights (and again trained while training feature weights). This
#' has shown to speed up total training time.
#' @param train_uncensored if TRUE, one round of training will occur assuming all event times are uncensored. This is done due to the non-convexity issue
#' that arises in the presence of censored data. However if ALL data is censored we recommend setting this option to FALSE as it has shown to give poor
#' results in this case.
#' @param seed_weights the initialization weights for the biases and the features. If left as NULL all weights are initialized to zero. If seed_weights are
#' specified then either nintervals or time_points must also be specified. The length of seed_weights should correspond to (number of features + 1)*(length of
#' time_points) = (number of features + 1)*(nintervals + 1).
#' @param threshold The threshold for the convergence tolerance (in the objective function) when training the feature weights.
#'  This threshold will be passed to \link[stats]{optim}.
#' @param maxit The maximum iterations to run for MTLR. This parameter will be passed to \link[stats]{optim}.
#' @param lower The lower bound for L-BFGS-B optimization. This parameter will be passed to \link[stats]{optim}.
#' @param upper The upper bound for L-BFGS-B optimization. This parameter will be passed to \link[stats]{optim}.
#' @details This function allows one to train an MTLR model given a dataset containing survival data. mtlr uses the Limited-Memory
#' Broyden–Fletcher–Goldfarb–Shanno (L-BFGS-B) approximation method to train feature weights. This training is outsourced to the internal
#' \link[stats]{optim} function in R. Currently only a few parameters (namely threshold, maxit,lower, upper) of optim are supported, more will
#' likely become available in the future.
#'
#' Weights are initialized to 0 prior to training. Under default settings, the bias weights
#' will be trained before considering feature weights. As Yu et al. (2011) specified, the introduction of censored observations creates a non-convex
#' loss function. To address this, weights are first trained assuming all event times are \emph{uncensored}. Once these starting weights have
#' been trained another round of training is performed using the true values of the event indicator (censored/uncensored). However, in the event of
#' all censored data this has shown to negatively effect the results. If all data is censored (either left, right, or interval2) we suggest setting
#' train_uncensored = FALSE.
#'
#' Yu et al. (2011) actually suggested two regularization parameters, C1 to control the size of the feature weights and C2 to control the smoothness.
#' In Ping Jin's masters thesis (Using Survival Prediction Techniques to Learn Consumer-Specific Reservation Price Distributions) he showed that C2
#' is not required for smoothness and C1 will suffice (Appendix A.2) so we do not support the C2 parameter in this implementation.
#'
#' If an error occurs from optim it is likely the weights are getting too large. Including fewer time points (or specifying better time points) in
#' addition to changing the lower/upper bounds of L-BFGS-B may resolve these issues. The most common failure has been that the objective value sees
#' infinite values due to extremely large feature weights.
#'
#'\strong{Censored data:} Right, left, and interval censored data are all supported both separately and mixed. The convention to input these types of
#'data follows the \link[survival]{Surv} object format.
#'Per the Surv documentation, "The [interval2] approach is to think of each observation as a time interval with (-infinity, t) for left censored,
#'(t, infinity) for right censored, (t,t) for exact and (t1, t2) for an interval. This is the approach used for type = interval2.
#'Infinite values can be represented either by actual infinity (Inf) or NA." See the examples below for an example of inputting this type of data.
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
#' @examples
#' #Access the Surv function and the leukemia/lung dataset.
#' library(survival)
#' simple_mod <- mtlr(Surv(time,status)~., data = leukemia)
#' simple_mod
#'
#' bigger_mod <- mtlr(Surv(time,status)~., data = lung)
#' bigger_mod
#'
#' #Note that observations with missing data were removed:
#' nrow(lung)
#' nrow(bigger_mod$x)
#'
#'
#' # Mixed censoring types
#' time1 = c(NA, 4, 7, 12, 10, 6, NA, 3) #NA for right censored
#' time2 = c(14, 4, 10, 12, NA, 9, 5, NA) #NA for left censored
#' #time1 == time2 indicates an exact death time. time2> time1 indicates interval censored.
#' set.seed(42)
#' dat = cbind.data.frame(time1, time2, importantfeature = rnorm(8))
#' formula = Surv(time1,time2,type = "interval2")~.
#' mixedmod = mtlr(formula, dat)
#'
#' @seealso
#' \code{\link[MTLR]{predict.mtlr}} \code{\link[MTLR]{mtlr_cv}} \code{\link[MTLR]{plot.mtlr}}  \code{\link[MTLR]{plotcurves}}
#' @export
mtlr <- function(formula,
                 data,
                 time_points = NULL,
                 nintervals = NULL,
                 normalize = T,
                 C1 = 1,
                 train_biases = T,
                 train_uncensored =T,
                 seed_weights = NULL,
                 threshold = 1e-05,
                 maxit = 5000,
                 lower = -15,
                 upper = 15){
  cl <- match.call() #Save a copy of the function call.
  if(any(dim(data)==0)){
    stop("Dimensions of the dataset must be non-zero.")
  }
  #Data setup
  mf <- stats::model.frame(formula = formula, data)
  Terms <- attr(mf, "terms")
  xlevels <- stats::.getXlevels(Terms, mf)

  x <- stats::model.matrix(Terms, data = mf)
  x <- x[,-1,drop=FALSE] #Remove intercept term -- We will handle biases later.
  y <- stats::model.response(mf)

  if (!survival::is.Surv(y))
    stop("The response must be a Surv object.")
  type = attr(y,"type")
  if(!type %in% c("right", "left","interval","interval2"))
    stop("Currently only right, left,interval, and interval2 censored data is supported. See the Details section of ?mtlr.")

  time <- y[,1]
  if(any(time < 0))
    stop("All event times must be non-negative.")
  if(type %in% c("interval","interval2")){
    time2 = y[,2]
    delta <- y[,3] #Death indicator (right censored = 0, death = 1, left censored = 2, interval censored = 3)
  }else{
    delta <- y[,2] #Death indicator (right censored = 0, death = 1, left censored = 2, interval censored = 3)
  }
  if(C1 < 0)
    stop("C1 must be non-negative.")
  if(threshold <= 0)
    stop("The threshold must be positive.")
  #Extra checks for when seed_weights is not NULL.
  if(!is.null(seed_weights)){
    if(is.null(time_points) & is.null(nintervals))
      stop("When using 'seed_weights' one of 'time_points' or 'nintervals' must also be specified")
    if(!is.null(time_points)){
      if(length(seed_weights) != (ncol(x)+1)*length(time_points))
        stop("The length of seed_weights must equal the product of the (number of features + 1) and the length of time_points.")
    }
    else if(!is.null(nintervals))
      if(length(seed_weights) != (ncol(x)+1)*(nintervals +1))
        stop("The length of seed_weights must equal the product of the (number of features + 1) and (nintervals + 1)")
  }
  if(is.null(nintervals) & is.null(time_points))
    nintervals <- ceiling(sqrt(nrow(y))) #Default number of intervals is the sq. root of the number of observations.
  if(normalize){
    x <- scale(x)
    scale_centers <- attr(x, "scaled:center")
    scale_scale <- attr(x,"scaled:scale")
    scales <- list(center= scale_centers, sd = scale_scale)
    if(any(scale_scale ==0)){
      zeroVar = names(scale_scale)[which(scale_scale ==0)]
      warning(paste("The feature(s): ", paste(zeroVar, collapse = " "), " have zero variance. These features have been set to zero to ensure a feature weight of zero."))
      x[,which(scale_scale ==0)] = 0
    }
  }else{
    scales <- NULL
  }

  ##############################################################
  #Prepare data for MTLR Rcpp functions.
  #The Rcpp functions are built to work on data that comes in with observations ordered by their event status
  #(censored observations first).
  #Here we order our data by the censor status. However, since left and interval censored are coded as 2,3
  #we make these negative first.
  delta = ifelse(delta == 1,1,-delta)
  ord <- order(delta)
  if(ncol(x)){
    x <- x[ord,1:ncol(x),drop=FALSE]
  }

  if(!is.null(nintervals) & is.null(time_points)){
    m <- nintervals + 1   #The number of time points to evaulate will be the number of time intervals + 1.
    quantiles <- seq(0,1,length.out = m+2)[-c(1,m+2)] #We will select time point based on the distribution of the times.
    time_points <- unname(stats::quantile(time, quantiles))
  } #Otherwise time_points is not null and we use the user inputted values
  if(!is.numeric(time_points))
    stop("time_points must be a numeric vector")



  time_points <- time_points[!duplicated(time_points)]#Small data (or common event times) we will have duplicated time points -- we remove them.

  #We make a matrix where each column is a vector of indicators if an observation is dead at each time point, e.g. (0,0,,...,0,1,1,...1).
  #(See 'Learning Patient-Specific Cancer Survival Distributions as a Sequence of Dependent Regressors' Page 3.)
  if(type =="right"){
    y_right <- matrix(as.numeric(Reduce(c,Map(function(ind) time[delta == 0] < c(time_points,Inf)[ind], seq_along(c(time_points,Inf))))),
                      nrow = length(time_points)+1, byrow = T)
    y_uncen <- matrix(as.numeric(Reduce(c,Map(function(ind) time[delta == 1] < c(time_points,Inf)[ind], seq_along(c(time_points,Inf))))),
                      nrow = length(time_points)+1, byrow = T)
    y_matrix <- cbind(y_right, y_uncen)
  } else if(type == "left"){
    y_left <- matrix(as.numeric(Reduce(c,Map(function(ind) time[delta == 0] >= c(0,time_points)[ind],
                                             seq_along(c(time_points,Inf))))), nrow = length(time_points) +1, byrow = T)
    y_uncen <- matrix(as.numeric(Reduce(c,Map(function(ind) time[delta == 1] < c(time_points,Inf)[ind], seq_along(c(time_points,Inf))))),
                      nrow = length(time_points)+1, byrow = T)
    y_matrix <- cbind(y_left,y_uncen)
  } else{
    y_int <- matrix(as.numeric(Reduce(c,Map(function(ind) time[delta == -3] <= c(time_points,Inf)[ind] & time2[delta == -3] > c(0,time_points)[ind],
                                            seq_along(c(time_points,Inf))))), nrow = length(time_points)+1, byrow = T)
    y_left <- matrix(as.numeric(Reduce(c,Map(function(ind) time[delta == -2] >= c(0,time_points)[ind],
                                             seq_along(c(time_points,Inf))))), nrow = length(time_points) +1, byrow = T)
    y_right <- matrix(as.numeric(Reduce(c,Map(function(ind) time[delta == 0] < c(time_points,Inf)[ind], seq_along(c(time_points,Inf))))),
                      nrow = length(time_points)+1, byrow = T)
    y_uncen <- matrix(as.numeric(Reduce(c,Map(function(ind) time[delta == 1] < c(time_points,Inf)[ind], seq_along(c(time_points,Inf))))),
                      nrow = length(time_points)+1, byrow = T)
    y_matrix <- cbind(y_int, y_left, y_right, y_uncen)
  }

  #We first train the biases (by setting feature parameters to zero (dAsZero)). Then we train the feature values as if all patients
  #were uncensored (this creates "good" starting values since with censored patients the objective is non-convex). Then we finally
  #train the data with their true censor statuses.
  threshold_factor <- threshold/.Machine$double.eps
  censInd = ifelse(delta <1,0,1)
  if(is.null(seed_weights)){
    if(train_biases){
      zero_matrix <- matrix(0,ncol = ncol(x), nrow = nrow(x))   #We create a zero_matrix to train the biases

      bias_par <- stats::optim(par = rep(0,length(time_points)*(ncol(x) +1)),fn = mtlr_objVal,gr = mtlr_grad, yval = y_matrix,
                               featureValue = zero_matrix, C1=C1, delta = rep(1,nrow(x)),
                               method = "L-BFGS-B", lower = lower, upper = upper, control = c(maxit = maxit, factr = threshold_factor))
      if(bias_par$convergence == 52)
        stop(paste("Error occured while training MTLR. Optim Error: ", bias_par$message))
    }else{
      bias_par <- list(par = rep(0,length(time_points)*(ncol(x) +1)))
    }
    if(train_uncensored){
      params_uncensored <- stats::optim(par = bias_par$par,fn = mtlr_objVal, gr = mtlr_grad, yval = y_matrix, featureValue = x, C1 = C1, delta = rep(1,nrow(x)),
                                        method = "L-BFGS-B", lower = lower, upper = upper, control = c(maxit = maxit, factr = threshold_factor))
      if(params_uncensored$convergence == 52)
      stop(paste("Error occured while training MTLR. Optim Error: ", params_uncensored$message))
    }else{
      params_uncensored = bias_par
    }

    final_params <- stats::optim(par = params_uncensored$par,fn = mtlr_objVal,gr = mtlr_grad, yval = y_matrix, featureValue = x, C1 = C1,delta = sort(censInd),
                                 method = "L-BFGS-B", lower = lower, upper = upper, control = c(maxit = maxit, factr = threshold_factor))
    if(final_params$convergence == 52)
      stop(paste("Error occured while training MTLR. Optim Error: ", final_params$message))
  } else{
    final_params <- stats::optim(par = seed_weights,fn = mtlr_objVal,gr = mtlr_grad, yval = y_matrix, featureValue = x, C1 = C1,delta = sort(censInd),
                                 method = "L-BFGS-B", lower = lower, upper = upper, control = c(maxit = maxit, factr = threshold_factor))
  }


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







