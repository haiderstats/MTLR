
print.mtlr <- function(x, digits = max(3, getOption("digits")), ...){
  cat("\nCall: ", deparse(x$Call), "\n\n")
  cat("\nWeights:\n" )
  print(x$weight_matrix,digits = digits)
}

predict.mtlr <- function(object, newdata, type = c("response"), add_zero = T,...){

  if(missing(newdata)){
    newframe = object$x
  }
  if(!missing(newdata)){
    Terms <- object$Terms
    Terms <- stats::delete.response(Terms)
    newframe <- stats::model.matrix(Terms, data=newdata,
                             xlev=object$xlevels)
    newframe <- newframe[,-1] #Remove intercept term.
    if(!is.null(object$scale)){
      newframe <- scale(newframe, center = object$scale$center, scale = object$scale$sd)
    }
  }

  surv_probs <- mtlr_predict(c(object$weight_matrix), newframe)
  #Issue due to machine precision, we get survival probabilities of 1+e-16. So here we adjust for that.
  surv_probs[surv_probs > 1] <- 1
  if(add_zero){
      time_points <- c(0,object$time_points)
      surv_probs <- rbind(1,surv_probs)
  }else{
    time_points <- object$time_points
  }
  surv_curves <- cbind.data.frame(time = time_points,surv_probs)
  surv_curves
}




