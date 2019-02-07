#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace std;

//Author: Humza Haider
//Email: hshaider@ualberta.ca

//Please see www.cs.cornell.edu/~cnyu/papers/nips11_survival.pdf for relevant information for objective value and gradient calculations.

/*mltr_objVal arguments:
 * (1) params:  The parameters for the MTLR model. These must be sorted with all the biases appearing first followed by the parameters
 * for the first feature at all time points, second feature for all time points, etc. E.g. with m time points and p features we have:
 * b1, b2,...,bm, theta1,1 ,....theta1,m, theta2,1, ... theta2,m, ..., theta_p,1,... theta p,m
 * (2) yval: A matrix with rows indicating the time points and columns being the patients status at that time point. E.g. for 6 time
 * intervals (5 time points)and for a patient who died in the third interval we would have 0,0,1,1,1,1 as a column. yval should be a matrix of all these
 * columns. (Should be N x m+1 where N is the number of patients and m is the number of time points). Intervals are (0,t1), (t1,t2),...,(tm-1, tm), (tm,Inf)
 * For censored individuals this will correspond to the censoring function (a value of 1 if the patients were censored at that time). For someone
 * right censored at between t3 and t4: (0,0,0,1,1,1)
 * left censored between t2 and t3: (1,1,0,0,0,0)
 * interval censored start between t2 and t3 ending between t4 and t5: (0,1,1,1,1,0)
 * died between t4 and t5: c(0,0,0,0,1,1)
 * (3) featureValue: These should be the matrix of feature values. A row indicates a single patient (corresponding to the first COLUMN of yval).
 * (4) C1: The regularzation parameter
 * (5) delta: A vector indicating is a patient is uncensored or censored (delta = 1 --> uncensored).
 */
//We will assume that the data is ordered such that all censored observations come first (type of censored observation does not mattter).
//We could abstract (fairly easily) away from this requirement, however, since this is the bottleneck of the code we try to outsource operations which
// would need to be repeated many times (i.e. ordering yval and feature val by delta). Additionally, we assume that delta =0 means any type of censoring
// not just right censoring (the code is written so all censored instances are handled the same but dependent on yval).

// [[Rcpp::export]]
double mtlr_objVal(arma::rowvec params, arma::mat yval, arma::mat featureValue, double C1,arma::vec delta) {
  //For the objective value see Equation 3 of the MTLR paper.

  //We passed params in as a NumericVector because these are easier to parse and turn into a matrix.
  double valToReturn = 0;
  int N = featureValue.n_rows;
  int m = yval.n_rows-1;

  arma::vec biases(m);
  for(int i =0; i < m; i++){
    biases(i) = params(i);
  }

  arma::mat thetas(m,featureValue.n_cols);
  int count = m;
  for(int j = 0; j < thetas.n_cols; j++) {
    for(int i =0; i < thetas.n_rows; i++){
      thetas(i,j) = params(count);
      count++;
    }
  }

  //score will correspond to the f_Theta function given by Yu. et al (2011)
  arma::vec score(m,arma::fill::ones);
  arma::vec cumScore(m+1, arma::fill::ones);

  //normalizingScore corresponds to f_Theta in the denominator of the pmf (Second equation on page 4 of Yu et al.)
  arma::vec normalizingScore(m+1, arma::fill::ones);

  arma::uvec oneIndex;

  //We will first sum the objective value for censored observations and then add in values for the uncensored observations.
  //By doing this we require that yval and featureVal have been sorted such that all censored observations come first.

  //Censored Piece
  double NCens = sum(1-delta);
  double censorScore;
  double biggestCensVal;
  //For the censored individuals (see Equation 4)... and take the log to get the log-likelihood.
  for(int i=0;i<NCens;i++){
    score =  thetas * featureValue.row(i).t() + biases;
    //Because we need to include the last score of 0 (for the vector of all 0s) we resize our score vector to be
    //1 unit bigger. Doing do automatically inserts a 0 at the end of our score vector.
    score.resize(m+1);
    cumScore  = reverse(cumsum(reverse(score)));
    // log trick: https://timvieira.github.io/blog/post/2014/02/11/exp-normalize-trick/
    double biggestVal = cumScore.max() > 0 ? cumScore.max() : 0;
    normalizingScore = cumScore - biggestVal;
    double normalizationTerm = biggestVal + log(sum(exp(normalizingScore)));

    //We catch if the vector is all 0's and if so we skip to assigning a score of 0.
    oneIndex = find(yval.col(i) == 1);

    if(oneIndex.size() ==0){
      censorScore = 0;
    }else{
      biggestCensVal  = cumScore.elem(oneIndex).max() > 0 ? cumScore.elem(oneIndex).max() : 0 ;
      cumScore -= biggestCensVal;
      //we multiply cumScore by ycol to 0 out the other large values in cumScore. Note while this will make exp(0) = 1,
      //it will again be zeroed out by multiplying by y.col(i) again.
      cumScore %= yval.col(i);
      censorScore = biggestCensVal + log(sum(yval.col(i) % exp(cumScore)));
    }

    valToReturn += censorScore - normalizationTerm;

  }
  //For the uncensored individuals.
  for(int i=NCens;i<N;i++){
    score =  thetas * featureValue.row(i).t() + biases;
    //Because we need to include the last score of 0 (for the vector of all 0s) we resize our score vector to be
    //1 unit bigger. Doing do automatically inserts a 0 at the end of our score vector.
    score.resize(m+1);
    cumScore  = reverse(cumsum(reverse(score)));
    double biggestVal  = cumScore.max();
    cumScore -= biggestVal;
    double normalizationTerm = biggestVal + log(sum(exp(cumScore)));
    double uncensoredScore = sum(yval.col(i).t() * score);
    valToReturn += uncensoredScore - normalizationTerm;
  }

  valToReturn = (C1/2)*accu(square(thetas)) - valToReturn/N;
  return valToReturn;
}



//See mtlr_objVal for argument definition.
//Here we will assume that the data is ordered such that all censored observations come first. We could abstract (fairly easily)
//away from this requirement but then we would be doing the same operation multiple times within the bottleneck of the code.
// [[Rcpp::export]]
arma::rowvec mtlr_grad(arma::rowvec params, arma::mat yval, arma::mat featureValue, double C1, arma::vec delta) {
  double N = featureValue.n_rows;
  double P = featureValue.n_cols;
  int m = yval.n_rows -1;

  arma::vec biases(m);
  for(int i =0; i < m; i++){
    biases(i) = params(i);
  }

  arma::mat thetas(m,featureValue.n_cols);
  int count = m;
  for(int j = 0; j < thetas.n_cols; j++) {
    for(int i =0; i < thetas.n_rows; i++){
      thetas(i,j) = params(count);
      count++;
    }
  }

  arma::vec BiasGradient = arma::zeros<arma::vec>(biases.n_elem);
  arma::mat ThetaGradient = arma::zeros<arma::mat>(thetas.n_rows, thetas.n_cols);
  arma::rowvec Xval(P, arma::fill::ones);

  //NormNumer is the numerator that results from the gradient of the normalization term in equation (3).
  //ThetaNorm is the gradient of the normalization term with respect to one of the feature weights.
  //partialThetaGrad is the gradient for a single observation. The sum of these gradients is the entire gradient for the feature weights (ThetaGradient).
  arma::mat  NormNumer(m,1,arma::fill::ones), ThetaNorm(m,P,arma::fill::ones), partialThetaGrad(m,P,arma::fill::ones);

  //partialBiasGrad is partialThetaGrad but for the biases.
  //score will correspond to the f_Theta function given by Yu. et al (2011) and cumScore will the the cumulative sum of scores.
  //BiasNorm is the gradient of the normalization term with respect to one of the bias weights.
  arma::vec partialBiasGrad(m,arma::fill::ones), score(m, arma::fill::ones), cumScore(m, arma::fill::ones), BiasNorm(m,arma::fill::ones);

  arma::uvec oneIndex;
  //NormDenom is the gradient of the denominator of the normalization term
  double NormDenom;

  //Censored Piece
  //partialBiasGrad is partialThetaGrad but for the biases.
  //BiasCensNumer is the gradient of the numerator of the normalization term with respect to one of the bias weights of a censored individual.
  //cumCensScore is the cumulative score for a censored observation.
  arma::vec BiasCensNumer(m,arma::fill::ones), cumCensScore(m, arma::fill::ones);

  //CensNormNumer is the numerator that results from the gradient of the normalization term in equation (4) for right censored. This generalizes to
  //other censored instances as well.
  //CensThetaNorm is the gradient of the normalization term with respect to one of the feature weights for a censored indiviual.
  arma::mat CensNormNumer(m,1,arma::fill::ones), CensThetaNorm(m,P,arma::fill::ones);

  //CensNormDenom is the gradient of the denominator of the normalization term for a censored observation
  double CensNormDenom;
  double NCens = sum(1-delta);
  double biggestCensVal;

  //Begin censored observation gradients.
  for(int i=0;i<NCens;i++){
    Xval = featureValue.row(i);
    score =  (thetas * Xval.t()) + biases;
    //Because we need to include the last score of 0 (for the vector of all 0s) we resize our score vector to be
    //1 unit bigger. Doing do automatically inserts a 0 at the end of our score vector.
    score.resize(m+1);
    cumScore  = reverse(cumsum(reverse(score)));

    oneIndex = find(yval.col(i) == 1);
    //Using the log trick, see https://timvieira.github.io/blog/post/2014/02/11/exp-normalize-trick/
    //We check for the sequence of all 0's first since this is an annoying edge case.
    if(oneIndex.size() ==0){
      cumCensScore = cumScore % yval.col(i);
      biggestCensVal = 0;
    }else{
      biggestCensVal  = cumScore.elem(oneIndex).max() > 0 ? cumScore.elem(oneIndex).max() : 0 ;
      cumCensScore = cumScore - biggestCensVal;
      //we multiply firstPieceVal by ycol to 0 out the other large values in cumScore. Note while this will make exp(0) = 1,
      //it will again be zeroed out by multiplying by y.col(i) again.
      cumCensScore %= yval.col(i);
    }

    CensNormNumer= cumsum(yval.col(i) % exp(cumCensScore));
    CensNormDenom = CensNormNumer.max();
    //The numerator had an extra row incase the observation was censored for the last time interval.
    //We needed that for calculating the denominator but need to remove it before calculating gradients
    //as it doesn't belong in the numerator.
    CensNormNumer.shed_row(m);

    BiasCensNumer = CensNormNumer/CensNormDenom;
    CensNormNumer *= Xval;
    CensThetaNorm = CensNormNumer/CensNormDenom;

    double biggestVal;
    // We need to shed the extra row we assigned to cumScore to handle censored data since the normalization term doesnt need this.
    cumScore.shed_row(m);
    biggestVal = cumScore.max() > 0 ? cumScore.max() : 0;
    cumScore -= biggestVal;

    NormNumer= cumsum(exp(cumScore));
    //Add an extra piece to account for the (0,0,0...,0,0) sequence.
    NormDenom = NormNumer[NormNumer.n_rows-1] +exp(-biggestVal);

    BiasNorm = NormNumer/NormDenom;
    NormNumer *= Xval;
    ThetaNorm = NormNumer/NormDenom;

    partialBiasGrad =  BiasCensNumer - BiasNorm;
    partialThetaGrad =  CensThetaNorm - ThetaNorm;
    BiasGradient -= partialBiasGrad/N;
    //Note we use += here since we later subtract from the regularization term.
    ThetaGradient += partialThetaGrad/N;
  }

  //yval has an extra row for censored data. Here we remove it before dealing with uncensored data.
  arma::mat yvalReduced = yval;
  yvalReduced.shed_row(m);
  //Begin uncensored observation gradients.
  for(int i=NCens;i<N;i++){
    Xval = featureValue.row(i);
    score =  (thetas * Xval.t()) + biases;

    cumScore  = reverse(cumsum(reverse(score)));
    double biggestVal;
    biggestVal = cumScore.max() > 0 ? cumScore.max() : 0;
    cumScore -= biggestVal;

    NormNumer= cumsum(exp(cumScore));
    NormDenom = NormNumer[NormNumer.n_rows-1] +exp(-biggestVal);
    BiasNorm = NormNumer/NormDenom;

    NormNumer *= Xval;
    ThetaNorm = NormNumer/NormDenom;
    partialBiasGrad =  sum(yvalReduced.col(i) - BiasNorm.each_col(), 1);
    partialThetaGrad =  yvalReduced.col(i) * Xval - ThetaNorm;
    BiasGradient -= partialBiasGrad/N;

    //Note we use += here since we later subtract from the regularization term.
    ThetaGradient += partialThetaGrad/N;
  }

  ThetaGradient = C1*thetas - ThetaGradient;
  arma::vec ThetaFlat = vectorise(ThetaGradient);
  return arma::join_cols<arma::mat>(BiasGradient, ThetaFlat).t();
}


//See mltr_objVal for argument definitions.
//This is used for predicting survival curves given feature weights.
// [[Rcpp::export]]
arma::mat mtlr_predict(arma::rowvec params,arma::mat featureValue){
  double N = featureValue.n_rows;
  double P = featureValue.n_cols;

  double m = params.size()/(P+1); // +1 to account for bias terms.
  arma::vec biases(m);
  for(int i =0; i < m; i++){
    biases(i) = params(i);
  }

  arma::mat thetas(m,featureValue.n_cols);
  int count = m;
  for(int j = 0; j < thetas.n_cols; j++) {
    for(int i =0; i < thetas.n_rows; i++){
      thetas(i,j) = params(count);
      count++;
    }
  }

  arma::mat resMat = thetas*featureValue.t();
  resMat.each_col() += biases;
  arma::mat B = exp(reverse(cumsum(reverse(resMat))));
  arma::mat Ones = arma::ones<arma::mat>(1,N);

  B.insert_rows(B.n_rows,Ones);
  arma::rowvec BDenoms = sum(B);
  arma::mat resultsMatrix = reverse(cumsum(reverse(B.each_row()/BDenoms)));
  resultsMatrix.shed_row(0);
  return resultsMatrix;
}













