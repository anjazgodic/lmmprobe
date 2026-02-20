// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <fstream>
#include <iostream>

using namespace Rcpp;

// [[Rcpp::export]]
List LM_by_col(const arma::vec &y, const arma::mat &X, double sigma2_lmm) {

  int n = X.n_rows, d = X.n_cols;

  arma::mat coef_mat(d, 2);
  coef_mat.fill(NA_REAL);
  arma::mat se(d, 2);
  se.fill(NA_REAL);
  arma::mat sig(d, 1);
  sig.fill(NA_REAL);

  // Precompute y terms
  double sum_y = arma::sum(y);
  double mean_y = sum_y / (double)n;
  double dot_y_y = arma::dot(y, y);
  double Syy = dot_y_y - n * mean_y * mean_y;

  for (int col = 0; col < d; ++col) {
    if (col % 1000 == 0)
      Rcpp::checkUserInterrupt();

    const arma::subview_col<double> x = X.col(col);

    double sum_x = arma::sum(x);
    double mean_x = sum_x / (double)n;
    double dot_x_x = arma::dot(x, x);
    double Sxx = dot_x_x - n * mean_x * mean_x;

    if (Sxx < 1e-12)
      continue;

    double dot_x_y = arma::dot(x, y);
    double Sxy = dot_x_y - n * mean_x * mean_y;

    double beta1 = Sxy / Sxx;
    double beta0 = mean_y - beta1 * mean_x;

    double RSS = Syy - beta1 * Sxy;
    if (RSS < 0)
      RSS = 0.0;

    double sig2 = RSS / (double)(n - 2);

    double var_beta1 = sigma2_lmm / Sxx;
    double var_beta0 = sigma2_lmm * (1.0 / n + (mean_x * mean_x) / Sxx);

    coef_mat(col, 0) = beta0;
    coef_mat(col, 1) = beta1;
    se(col, 0) = std::sqrt(var_beta0);
    se(col, 1) = std::sqrt(var_beta1);
    sig(col, 0) = std::sqrt(sig2);
  }

  return List::create(Named("Coefficients") = coef_mat, Named("StdErr") = se,
                      Named("sig") = sig);
}
