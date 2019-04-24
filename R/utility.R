#utility.R contains functions that are used by exported functions but are not exported themselves.


#Helper function for plot.mtlr -- Here we get the sum of the absolute values of the feature weights across time.
#This essentially measure how much "influence" each feature had on the survival of each observation.
get_param_influence <- function(mtlr_object){
  weights<- mtlr_object$weight_matrix
  weights <- weights[,-1, drop= FALSE] #Remove bias
  influence <- apply(weights,2,function(x) sum(abs(x)))
  return(influence)
}


#Helper function for predict.mtlr

#We need some type of predict function for survival curves - here we build a spline to fit the survival model curve. This spline is
#the montotone spline using the hyman filtering of the cubic Hermite spline method,
#see https://en.wikipedia.org/wiki/Monotone_cubic_interpolation. Also see help(splinefun).
#Note that for timepoints outside of the fit the spline is constant, i.e. all probabilities outside of the max predicted_times will be constant.
predict_prob <- function(survival_curve,predicted_times, time_to_predict){
  spline <- stats::splinefun(predicted_times, survival_curve, method = "hyman")
  maxTime <- max(predicted_times)
  predictedProbabilities <- rep(0, length(time_to_predict))
  maxProb <- which(time_to_predict > maxTime)
  if(length(maxProb) > 0){
    #We enforce that times outside of the max time are constant.
    predictedProbabilities[maxProb] <- spline(maxTime)
    #If time_to_predict is less than predicted_times then we will enforce this probability to be a max of 1.
    predictedProbabilities[-maxProb] <- pmin(spline(time_to_predict[-maxProb]),1)
  }
  else{
    predictedProbabilities <- pmin(spline(time_to_predict),1)
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
  spline <- stats::splinefun(predicted_times, survival_curve, method = "hyman")
  maxTime <- max(predicted_times)
  slope <- (1-spline(maxTime))/(min(predicted_times) - max(predicted_times))
  zeroProbabilitiyTime <- min( predicted_times[which(survival_curve ==0)], maxTime + (0-spline(maxTime))/slope)
  splineWithLinear <- function(time) ifelse(time < maxTime, spline(time),1 + time*slope)
  area <- stats::integrate(splineWithLinear,min(predicted_times), zeroProbabilitiyTime,subdivisions = 1000,rel.tol = .0001)[[1]]
  return(area)
}


predict_median <- function(survival_curve, predicted_times){
  #If all the predicted probabilities are 1 the median will be infinite.
  if(all(survival_curve==1)){
    return(Inf)
  }
  spline <- stats::splinefun(predicted_times, survival_curve, method = "hyman")
  minProb <- min(spline(predicted_times))
  maxProb <- max(spline(predicted_times))
  if(maxProb < 0.5) return(min(predicted_times))
  if(minProb < 0.5){
    maximumSmallerThanMedian <- predicted_times[min(which(survival_curve <.5))]
    minimumGreaterThanMedian <- predicted_times[max(which(survival_curve >.5))]
    splineInv <- stats::splinefun(spline(seq(minimumGreaterThanMedian, maximumSmallerThanMedian, length.out = 1000)),
                          seq(minimumGreaterThanMedian, maximumSmallerThanMedian, length.out = 1000))
    medianProbabilityTime <- splineInv(0.5)
  } else{
    maxTime <- max(predicted_times)
    slope <- (1-spline(maxTime))/(min(predicted_times) - max(predicted_times))
    medianProbabilityTime <- maxTime + (0.5-spline(maxTime))/slope
  }
  return(medianProbabilityTime)
}

loglik_loss <- function(object, newdata){
  #For the loss we need to compute losses differently for censored and uncensored patients.
  #For right censored patients the loss will correspond the the (log) survival probability assigned by the model at the time of censoring.
  #For left censored patients the loss will correspond the the (log) death probability (1- survival) assigned by the model at the time of censoring.
  #For interval censored patients the loss will correspond the the (log) probability assigned to the interval (probability of survival lower bound time -
  #probability of survival upper bound time).
  #For uncensored patients, we will consider the log of the probability assigned to the time interval when the patient died.
  #Then we take the negative of this loss and thus would like to minimize the loss.

  #Get the survival curves for all observations.
  curves <- predict.mtlr(object,newdata, type = "survivalcurve")


  #We need to know which observations are censored in the new data.
  Terms <- object$Terms
  mf <- stats::model.frame(Terms, data=newdata,xlev=object$xlevels)
  response <- stats::model.response(mf)

  type = attr(response,"type")
  event_times <- response[,1]

  if(type %in% c("interval","interval2")){
    event_times2 = response[,2]
    delta <- response[,3] #Death indicator (right censored = 0, death = 1, left censored = 2, interval censored = 3)
  }else{
    delta <- response[,2] #Death indicator (right censored = 0, death = 1, left censored = 2, interval censored = 3)
  }

  censor_ind <- which(delta!=1)
  uncen_ind <- which(delta==1)

  #Initialize loss to 0.
  logloss <- 0

  #Censored patients
  if(length(censor_ind)){
    if(type == "right" | type == "left"){
      censor_prob <- sapply(censor_ind,
                            function(index) predict_prob(curves[,index+1],
                                                         curves[,1],
                                                         event_times[index]
                            ))
      censor_curves = curves[,(censor_ind +1), drop=F]
      censor_event = event_times[censor_ind]
      censor_prob <- mapply(function(x,y) predict_prob(x, curves[,1],y),
                            censor_curves, censor_event)
      #Currently we add a value of 0.00001 to all values since otherwise we might take the log of 0.
      #This is clearly not ideal but shouldn't make a large impact especially since this is only used for selecting C1.
      if(type == "right"){
        #For right censor we have a loss of the probability of surviving to censor time, log(1) is optimal, log (0) is worst case.
        logloss <- logloss - sum(log(censor_prob + 1e-05))
      }else{
        #For left censor we have a loss of the probability of dying before censor time, log(1) is optimal, log (0) is worst case.
        logloss <- logloss - sum(log(1-censor_prob + 1e-05))
      }
    }else{
      left_right <- which(delta %in% c(0,2))
      interval <- which(delta %in% c(3))


      left_right_curves = curves[,(left_right +1), drop=F]
      left_right_event = event_times[left_right]
      left_right_prob <- mapply(function(x,y) predict_prob(x, curves[,1],y),
                            left_right_curves, left_right_event)

      interval_curves = curves[,(interval +1), drop=F]
      interval_eventL = event_times[interval]
      interval_eventU = event_times2[interval]

      interval_probL <- mapply(function(x,y) predict_prob(x, curves[,1],y),
                               interval_curves, interval_eventL)
      interval_probU <- mapply(function(x,y) predict_prob(x, curves[,1],y),
                               interval_curves, interval_eventU)

      left_right_prob <- ifelse(left_right %in% which(delta==0), left_right_prob, 1 - left_right_prob)
      if(length(interval_probL)){
        interval_prob <- interval_probL - interval_probU
        censor_prob <- c(left_right_prob, interval_prob)
      }else{
        censor_prob = left_right_prob
      }
      logloss <- logloss - sum(log(censor_prob + 1e-05))
    }

  }

  #Uncensored patients
  if(length(uncen_ind)){
    uncen_curves <- curves[,uncen_ind+1,drop=FALSE]
    uncen_curves <- rbind(uncen_curves, 0)
    #Get the probability of the event occuring in every interval.
    pmf_probs <- -apply(uncen_curves,2,diff)

    #This will find the index of which interval the event occured we then index pmf_probs the by observation index
    #and the event index to generate the vector of probabilities assigend to each observations respective event time.
    event_ind <- sapply(event_times[uncen_ind], function(x) findInterval(x, curves[,1]))
    index_col <- 1:length(uncen_ind)
    index_mat <- matrix(c(event_ind,index_col),ncol = 2)
    probs <- pmf_probs[index_mat]
    logloss <- logloss  - sum(log(probs+1e-05))
  }
  return(logloss/nrow(newdata))
}


concordance_loss <- function(object, newdata){
  #Get the survival curves for all observations.
  preds <- -1*predict.mtlr(object,newdata, type = "median")
  Terms <- object$Terms
  mf <- stats::model.frame(Terms, data=newdata,xlev=object$xlevels)
  response <- stats::model.response(mf)
  conc <- -1*unname(survival::survConcordance(response~preds)$concordance)
  return(conc)
}





































