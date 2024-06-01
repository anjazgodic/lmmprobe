// [[Rcpp::depends(RcppArmadillo)]]

#include <fstream>
#include <iostream>
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
List MVM(const arma::mat X,const arma::vec v) {
  
  arma::mat Z = X.each_row()%arma::trans(v);
  
  return List::create(Named("Res") = Z);
}







   
  