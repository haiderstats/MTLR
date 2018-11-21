#utility.R contains functions that are used by exported functions but are not exported themselves.


#Helper function for plot.mtlr -- Here we get the sum of the absolute values of the feature weights across time.
#This essentially measure how much "influence" each feature had on the survival of each observation.
get_param_influence <- function(mtlr_object){
  weights= mtlr_object$weight_matrix
  weights = weights[,-1] #Remove bias
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
predict_prob <- function(survivalCurve,predictedTimes, timeToPredict){
  spline = stats::splinefun(predictedTimes, survivalCurve, method = "hyman")
  maxTime = max(predictedTimes)
  slope = (1-spline(maxTime))/(0 - max(predictedTimes))
  predictedProbabilities = rep(0, length(timeToPredict))
  linearChange = which(timeToPredict > maxTime)
  if(length(linearChange) > 0){
    predictedProbabilities[linearChange] = pmax(1 + timeToPredict[linearChange]*slope,0)
    predictedProbabilities[-linearChange] = spline(timeToPredict[-linearChange])
  }
  else{
    predictedProbabilities = spline(timeToPredict)
  }
  return(predictedProbabilities)
}


#We calculate the mean and median survival times assuming a monotone spline fit of the survival curve points.
predict_mean <- function(survivalCurve, predictedTimes){
  #If all the predicted probabilities are 1 the integral will be infinite. For this reason we slightly decrease the
  #last value.
  if(all(survivalCurve==1)){
    return(Inf)
  }
  spline = stats::splinefun(predictedTimes, survivalCurve, method = "hyman")
  maxTime = max(predictedTimes)
  slope = (1-spline(maxTime))/(0 - max(predictedTimes))
  zeroProbabilitiyTime = min( predictedTimes[which(survivalCurve ==0)], maxTime + (0-spline(maxTime))/slope)
  splineWithLinear = function(time) ifelse(time < maxTime, spline(time),1 + time*slope)
  area = stats::integrate(splineWithLinear,0, zeroProbabilitiyTime,subdivisions = 1000,rel.tol = .0001)[[1]]
  return(area)
}

predict_median <- function(survivalCurve, predictedTimes){
  #If all the predicted probabilities are 1 the integral will be infinite.
  if(all(survivalCurve==1)){
    return(Inf)
  }
  spline = stats::splinefun(predictedTimes, survivalCurve, method = "hyman")
  minProb = min(spline(predictedTimes))
  if(minProb < .5){
    maximumSmallerThanMedian = predictedTimes[min(which(survivalCurve <.5))]
    minimumGreaterThanMedian = predictedTimes[max(which(survivalCurve >.5))]
    splineInv = stats::splinefun(spline(seq(minimumGreaterThanMedian, maximumSmallerThanMedian, length.out = 1000)),
                          seq(minimumGreaterThanMedian, maximumSmallerThanMedian, length.out = 1000))
    medianProbabilityTime = splineInv(0.5)
  }
  else{
    maxTime = max(predictedTimes)
    slope = (1-spline(maxTime))/(0 - max(predictedTimes))
    medianProbabilityTime = maxTime + (0.5-spline(maxTime))/slope
  }
  return(medianProbabilityTime)
}

