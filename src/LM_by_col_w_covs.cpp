// [[Rcpp::depends(RcppArmadillo)]]

#include <fstream>
#include <iostream>
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
List LM_by_col_w_covs(const arma::vec y, const arma::mat X, const arma::colvec X_adj, double sigma2_lmm) {
  
  int n = X.n_rows, d = X.n_cols;
  
  arma::vec X_mean = mean(X,1);
  arma::vec Xp = X_mean; 
  
  arma::mat coef_mat(d,3);
  arma::mat se(d,3);
  arma::mat sig(d,1);
  arma::mat X1(n,1,arma::fill::ones);
  
  for(int col=0;col<d; ++col){
    
    arma::mat X2=arma::join_rows(X1,X.col(col), X_adj);
    int n = X2.n_rows, k = X2.n_cols;
    
    arma::colvec coef = arma::solve(X2, y); 
    arma::colvec resid = y - X2*coef; 
    
    double sig2 = arma::as_scalar(arma::trans(resid)*resid/(n-k));
    arma::colvec stderrest = 
      arma::sqrt(sigma2_lmm * arma::diagvec( arma::inv(arma::trans(X2)*X2)) );
    double s = std::pow(sig2,0.5);
    
    coef_mat.row(col) = arma::trans(coef);
    se.row(col) = arma::trans(stderrest);
    sig.row(col) = s; 
  }
  
  return List::create(Named("Coefficients") = coef_mat,Named("StdErr") = se,Named("sig")=sig);
}







