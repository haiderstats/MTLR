



mltr <- function(formula,
                 data,
                 nintervals = NULL,
                 C1 = 1,
                 threshold = 1e-5,
                 maxit = 5000,
                 normalize = T){
  cl <- match.call() #Save a copy of the function call.

  #Data setup
  mf <- model.frame(formula = formula, data = data)
  x <- model.matrix(attr(mf, "terms"), data = mf)
  x <- x[,-1] #Remove intercept term -- We will handle biases later.
  y <- model.response(mf)

  if (!is.Surv(y))
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

  ##############################################################
  #Prepare data for MTLR Rcpp functions.
  #The Rcpp functions are built to work on data that comes in with observations ordered by their event status (censored observations first).
  #Here we order our data by the censor status.
  ord = order(delta)

  x = x[ord,]

  m = nintervals + 1   #The number of time points to evaulate will be the number of time intervals + 1.
  quantiles <- seq(0,1,length.out = m+2)[-c(1,m+2)] #We will select time point based on the distribution of the times.
  time_points <- unname(quantile(time, quantiles))


  time_points <- time_points[!duplicated(time_points)]#Small data (or common event times) we will have duplicated time points -- we remove them.

  zero_matrix = matrix(0,ncol = ncol(x), nrow = nrow(x))   #We create a zero_matrix to train the biases

  #We make a matrix where each column is a vector of indicators if an observation is dead at each time point, e.g. (0,0,,...,0,1,1,...1).
  #(See 'Learning Patient-Specific Cancer Survival Distributions as a Sequence of Dependent Regressors' Page 3.)
  y_matrix = matrix(1 - Reduce(c,Map(function(ind) time > time_points[ind], seq_along(time_points))), ncol = nrow(x), byrow = T)

  #We first train the biases (by setting feature parameters to zero (dAsZero)). Then we train the feature values as if all patients
  #were uncensored (this creates "good" starting values since with censored patients the objective is non-convex). Then we finally
  #train the data with their true censor statuses.
  bias_par = optim(par = rep(0,length(timePoints)*(ncol(d) +1)),fn = mtlr_objVal,gr = mtlr_grad, yval = yval,
                  featureVal = dAsZero,C1=C1, delta = sort(training$delta),
                  method = "L-BFGS-B", lower = -20,upper=20,control=c(maxit = 5000, factr = .45036e11))

  allParamsUnc = optim(par = biasPar$par,fn = mtlr_objVal,gr = mtlr_grad, yval = yval, featureVal = d,C1=C1,delta = rep(1,nrow(training)),
                       method = "L-BFGS-B", lower = -20,upper=20,control=c(maxit = 5000, factr = .45036e11))

  allParams = optim(par = allParamsUnc$par,fn = mtlr_objVal,gr = mtlr_grad, yval = yval, featureVal = d,C1=C1,delta = sort(training$delta),
                    method = "L-BFGS-B", lower = -20,upper=20,control=c(maxit = 5000, factr = .45036e11))


  return(0)
}







