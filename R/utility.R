#utility.R contains functions that are used by exported functions but are not exported themselves.


#Helper function for plot.mtlr -- Here we get the sum of the absolute values of the feature weights across time.
#This essentially measure how much "influence" each feature had on the survival of each observation.
get_param_influence <- function(mtlr_object){
  weights= mtlr_object$weight_matrix
  weights = weights[,-1, drop= FALSE] #Remove bias
  influence = apply(weights,2,function(x) sum(abs(x)))
  return(influence)
}


#Helper function for predict.mtlr

#We need some type of predict function for survival curves - here we build a spline to fit the survival model curve. This spline is
#the montotone spline using the hyman filtering of the cubic Hermite spline method,
#see https://en.wikipedia.org/wiki/Monotone_cubic_interpolation. Also see help(splinefun).
#Note that we make an alteration to the method because if the last two time points
#have the same probability (y value) then the spline is constant outside of the training data. We need this to be a decreasing function
#outside the training data so instead we take the linear fit of (0,1) and the last time point we have (p,t*) and then apply this linear
#function to all points outside of our fit.
predict_prob <- function(survival_curve,predicted_times, time_to_predict){
  spline = stats::splinefun(predicted_times, survival_curve, method = "hyman")
  maxTime = max(predicted_times)
  slope = (1-spline(maxTime))/(min(predicted_times) - max(predicted_times))
  predictedProbabilities = rep(0, length(time_to_predict))
  linearChange = which(time_to_predict > maxTime)
  if(length(linearChange) > 0){
    predictedProbabilities[linearChange] = pmax(1 + time_to_predict[linearChange]*slope,0)
    #If time_to_predict is less than predicted_times then we will enforce this probability to be a max of 1.
    #The spline fit does a linear line outside of the range of predicted_times.
    predictedProbabilities[-linearChange] = pmin(spline(time_to_predict[-linearChange]),1)
  }
  else{
    predictedProbabilities = pmin(spline(time_to_predict),1)
  }
  return(predictedProbabilities)
}


#We calculate the mean and median survival times assuming a monotone spline fit of the survival curve points.
predict_mean <- function(survival_curve, predicted_times){
  #If all the predicted probabilities are 1 the integral will be infinite. For this reason we slightly decrease the
  #last value.
  if(all(survival_curve==1)){
    return(Inf)
  }
  spline = stats::splinefun(predicted_times, survival_curve, method = "hyman")
  maxTime = max(predicted_times)
  slope = (1-spline(maxTime))/(min(predicted_times) - max(predicted_times))
  zeroProbabilitiyTime = min( predicted_times[which(survival_curve ==0)], maxTime + (0-spline(maxTime))/slope)
  splineWithLinear = function(time) ifelse(time < maxTime, spline(time),1 + time*slope)
  area = stats::integrate(splineWithLinear,0, zeroProbabilitiyTime,subdivisions = 1000,rel.tol = .0001)[[1]]
  return(area)
}


predict_median <- function(survival_curve, predicted_times){
  #If all the predicted probabilities are 1 the integral will be infinite.
  if(all(survival_curve==1)){
    return(Inf)
  }
  spline = stats::splinefun(predicted_times, survival_curve, method = "hyman")
  minProb = min(spline(predicted_times))
  maxProb = max(spline(predicted_times))
  if(maxProb < 0.5) return(min(predicted_times))
  if(minProb < 0.5){
    maximumSmallerThanMedian = predicted_times[min(which(survival_curve <.5))]
    minimumGreaterThanMedian = predicted_times[max(which(survival_curve >.5))]
    splineInv = stats::splinefun(spline(seq(minimumGreaterThanMedian, maximumSmallerThanMedian, length.out = 1000)),
                          seq(minimumGreaterThanMedian, maximumSmallerThanMedian, length.out = 1000))
    medianProbabilityTime = splineInv(0.5)
  } else{
    maxTime = max(predicted_times)
    slope = (1-spline(maxTime))/(min(predicted_times) - max(predicted_times))
    medianProbabilityTime = maxTime + (0.5-spline(maxTime))/slope
  }
  return(medianProbabilityTime)
}

foldme <- function(time, delta, nfolds,foldtype = c("fullstrat","censorstrat","random")){
  if(nfolds < 1)
    stop("Number of folds must be greater than 0.")
  type = match.arg(foldtype)
  foldIndex = switch(type,
                     fullstrat = {
                       Order= order(delta,time)
                       lapply(1:nfolds, function(x) Order[seq(x,length(time), by = nfolds)])
                     },
                     censorstrat = {
                       censored = sample(which(!delta))
                       uncensored = sample(which(!!delta))
                       Order= c(censored,uncensored)
                       lapply(1:nfolds, function(x) Order[seq(x,length(time), by = nfolds)])
                     },
                     random = {
                       Order = sample(length(time))
                       lapply(1:nfolds, function(x) Order[seq(x,length(time), by = nfolds)])
                     }
  )
  return(foldIndex)
}

log_loss <- function(object, newdata){
  #For the loss we need to compute losses differently for censored and uncensored patients.
  #For censored patients the loss will correspond the the (log) survival probability assigned by the model at the time of censoring.
  #For uncensored patients, we will consider the log of the probability assigned to the time interval when the patient died.
  #Then we take the negative of this loss and thus would like to minimize the loss.

  #Get the survival curves for all observations.
  curves = predict.mtlr(object,newdata, type = "survivalcurve")


  #We need to know which observations are censored in the new data.
  Terms <- object$Terms
  mf <- stats::model.frame(Terms, data=newdata,xlev=object$xlevels)
  response <- stats::model.response(mf)

  event_times <- response[,1]
  delta <- response[,2]

  censor_ind <- which(!delta)
  uncen_ind <- which(!!delta)

  #Initialize loss to 0.
  logloss <- 0

  #Censored patients
  if(length(censor_ind)){
    censor_prob <- sapply(censor_ind,
                          function(index) predict_prob(curves[,index+1],
                                                       curves[,1],
                                                       event_times[index]
                          ))
    #Currently we add a value of 0.00001 to all values since otherwise we might take the log of 0.
    #This is clearly not ideal but shouldn't make a large impact especially since this is only used for selecting C1.
    logloss <- logloss - sum(log(censor_prob + 1e-05))
  }

  #Uncensored patients
  if(length(uncen_ind)){
    uncen_curves <- curves[,uncen_ind+1,drop=FALSE]
    uncen_curves = rbind(uncen_curves, 0)
    #Get the probability of the event occuring in every interval.
    pmf_probs <- -apply(uncen_curves,2,diff)

    #This will find the index of which interval the event occured we then index pmf_probs the by observation index
    #and the event index to generate the vector of probabilities assigend to each observations respective event time.
    event_ind <- sapply(event_times[uncen_ind], function(x) findInterval(x, curves[,1]))
    index_col <- 1:sum(delta)
    index_mat <- matrix(c(event_ind,index_col),ncol = 2)
    probs <- pmf_probs[index_mat]
    logloss <- logloss  - sum(log(probs+1e-05))
  }
  return(logloss/nrow(newdata))
}





































