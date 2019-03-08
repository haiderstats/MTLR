#' Create folds for cross-validation.
#'
#' Create the test folds for k-fold cross validation. These cross-validation types differ from typical stratified cross-validation as this function also
#' considers the range of event times in the data.
#' @param time a vector of event times.
#' @param delta a vector of indicators for uncensored/censored data. The type of censoring here is not considered so it is suggested this function not
#' be used for data with mixed censoring types. The specific indicator value does not matter as long as censored and uncensored observations have different
#' values for their indicator.
#' @param nfolds The number of folds to create.
#' @param foldtype type of cross validation folds. Full stratification, "fullstrat", sorts observations by their event time and their event indicators
#' and numbers them off into folds. This effectively give each fold approximately the same number of uncensored observations as well as keeps the range
#' of time points as equivalent as possible across folds. This type of cross-validation is completely deterministic.
#' Censored stratification, "censorstrat", will put approximately the same number of uncensored observations in each fold but not pay any attention to
#' event time. This is partially stochastic. The totally random cross-validation, "random", randomly assigns observations to folds without considering
#' event time nor event status.
#' @return a list of size nfolds where each list component contains the indices of the test data for each fold.
#' @seealso \code{\link[MTLR]{mtlr_cv}}
#' @export
create_folds <- function(time, delta, nfolds, foldtype = c("fullstrat","censorstrat","random")){
  if(nfolds < 1)
    stop("Number of folds must be greater than 0.")
  type <- match.arg(foldtype)
  foldIndex <- switch(type,
                      fullstrat = {
                        Order<- order(delta,time)
                        lapply(1:nfolds, function(x) Order[seq(x,length(time), by = nfolds)])
                      },
                      censorstrat = {
                        censored <- sample(which(!delta))
                        uncensored <- sample(which(!!delta))
                        Order<- c(censored,uncensored)
                        lapply(1:nfolds, function(x) Order[seq(x,length(time), by = nfolds)])
                      },
                      random = {
                        Order <- sample(length(time))
                        lapply(1:nfolds, function(x) Order[seq(x,length(time), by = nfolds)])
                      }
  )
  return(foldIndex)
}
