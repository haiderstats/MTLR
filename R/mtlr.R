#' @useDynLib MTLR
#' @importFrom Rcpp sourceCpp
NULL

#' Primary MTLR function. Run the mode
#'
#' @param formula Surv formula for survival data.
#' @param data Surv data.
#' @param nintervals Number of time intervals.
#' @param C1 The regularization parameter
#' @param normalize Boolean specifying if normalization should occur. Strongly suggested.
#' @param threshold The threshold for optim
#' @param maxit The maximum interations for optim
#' @param lower The lower L-BFGS-B optim
#' @param upper The upper.
#' @param train_biases Boolean specifying whether to run an interation training biases first.
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
  x <- x[,-1] #Remove intercept term -- We will handle biases later.
  y <- stats::model.response(mf)

  if (!survival::is.Surv(y))
    stop("The response must be a Surv object.")
  if(attr(y, "type") != "right")
    stop("Currently only right censored data is supported.")

  time <- y[,1]
  delta <- y[,2] #Death indicator (censored = 0, death = 1)

  if(any(time < 0))
    stop("All event times must be positive.")
  if(max(delta) == 0)
    stop("All events are censored.")
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
  #The Rcpp functions are built to work on data that comes in with observations ordered by their event status (censored observations first).
  #Here we order our data by the censor status.
  ord = order(delta)

  x <- x[ord,]

  m <- nintervals + 1   #The number of time points to evaulate will be the number of time intervals + 1.
  quantiles <- seq(0,1,length.out = m+2)[-c(1,m+2)] #We will select time point based on the distribution of the times.
  time_points <- unname(stats::quantile(time, quantiles))


  time_points <- time_points[!duplicated(time_points)]#Small data (or common event times) we will have duplicated time points -- we remove them.

  #We make a matrix where each column is a vector of indicators if an observation is dead at each time point, e.g. (0,0,,...,0,1,1,...1).
  #(See 'Learning Patient-Specific Cancer Survival Distributions as a Sequence of Dependent Regressors' Page 3.)
  y_matrix = matrix(1 - Reduce(c,Map(function(ind) time[ord] > time_points[ind], seq_along(time_points))), ncol = nrow(x), byrow = T)

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
  x = x[order(ord),]
  y_matrix = y_matrix[,order(ord)]
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







