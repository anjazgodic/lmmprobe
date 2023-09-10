// [[Rcpp::depends(RcppArmadillo)]]

#include <fstream>
#include <iostream>
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
List Row_sum(const arma::mat X) {
  
  arma::vec Xp = sum(X,1);

  return List::create(Named("Rowsum") = Xp);
}







   
  