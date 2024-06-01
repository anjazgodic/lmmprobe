// [[Rcpp::depends(RcppArmadillo)]]

#include <fstream>
#include <iostream>
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
List Col_sum(const arma::mat X) {
  
  arma::rowvec Xp = sum(X,0);

  return List::create(Named("Colsum") = Xp);
}







   
  