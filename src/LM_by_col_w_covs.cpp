// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <fstream>
#include <iostream>

using namespace Rcpp;

// [[Rcpp::export]]
List LM_by_col_w_covs(const arma::vec y, const arma::mat X,
                      const arma::colvec X_adj, double sigma2_lmm) {

  int n = X.n_rows, d = X.n_cols;
  int k = 3; // intercept, x_col, X_adj

  arma::mat coef_mat(d, 3);
  arma::mat se(d, 3);
  arma::mat sig(d, 1);

  // Precomputations
  double sum_y = arma::sum(y);
  double dot_y_y = arma::dot(y, y);

  double sum_X_adj = arma::sum(X_adj);
  double dot_X_adj_X_adj = arma::dot(X_adj, X_adj);
  double dot_X_adj_y = arma::dot(X_adj, y);

  // Predictors: [W1, X_col, X_adj]
  // bXXt (3x3):
  // [ n,          sum(X_col),    sum(X_adj) ]
  // [ sum(X_col), dot(X,X),      dot(X, X_adj) ]
  // [ sum(X_adj), dot(X, X_adj), dot(X_adj, X_adj) ]

  arma::mat bXXt(3, 3);
  bXXt(0, 0) = n;
  bXXt(0, 2) = sum_X_adj;
  bXXt(2, 0) = sum_X_adj;
  bXXt(2, 2) = dot_X_adj_X_adj;

  arma::vec bXty(3);
  bXty(0) = sum_y;
  bXty(2) = dot_X_adj_y;

  for (int col = 0; col < d; ++col) {
    if (col % 1000 == 0)
      Rcpp::checkUserInterrupt();

    const arma::subview_col<double> x_col = X.col(col);

    double sum_x = arma::sum(x_col);
    double dot_x_x = arma::dot(x_col, x_col);
    double dot_x_y = arma::dot(x_col, y);
    double dot_x_X_adj = arma::dot(x_col, X_adj);

    // Update matrix
    bXXt(0, 1) = sum_x;
    bXXt(1, 0) = sum_x;
    bXXt(1, 1) = dot_x_x;
    bXXt(1, 2) = dot_x_X_adj;
    bXXt(2, 1) = dot_x_X_adj;

    // Update vector
    bXty(1) = dot_x_y;

    if (!bXXt.is_sympd())
      continue;
    arma::mat bXXt_inv = arma::inv_sympd(bXXt);

    arma::vec coef = bXXt_inv * bXty;

    // RSS = y'y - beta'bXty
    double dot_beta_bXty = arma::dot(coef, bXty);
    double RSS = dot_y_y - dot_beta_bXty;
    if (RSS < 0)
      RSS = 0;

    double sig2 = RSS / (double)(n - k);

    // SE = sqrt(sigma2_lmm * diag(inv(X'X)))
    arma::vec stderrest = arma::sqrt(sigma2_lmm * bXXt_inv.diag());

    coef_mat.row(col) = coef.t();
    se.row(col) = stderrest.t();
    sig(col, 0) = std::sqrt(sig2);
  }

  return List::create(Named("Coefficients") = coef_mat, Named("StdErr") = se,
                      Named("sig") = sig);
}
