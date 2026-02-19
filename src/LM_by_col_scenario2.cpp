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

  // Precompute centered Y
  double mean_y = arma::mean(y);
  arma::vec y_cen = y - mean_y;
  double Syy = arma::dot(y_cen, y_cen);

  // We iterate over columns. For each column x, we fit y ~ x.
  // Model: y = beta0 + beta1 * x + epsilon

  for (int col = 0; col < d; ++col) {

    // Access column without copy
    const arma::subview_col<double> x = X.col(col);

    double mean_x = arma::mean(x);
    arma::vec x_cen = x - mean_x;

    // Sxx = sum((xi - mean_x)^2)
    double Sxx = arma::dot(x_cen, x_cen);

    // Check for singularity (constant x)
    if (Sxx < 1e-12) {
      // Leaves as NA
      continue;
    }

    // Sxy = sum((xi - mean_x)*(yi - mean_y))
    double Sxy = arma::dot(x_cen, y_cen);

    double beta1 = Sxy / Sxx;
    double beta0 = mean_y - beta1 * mean_x;

    // Optimized RSS calculation:
    // RSS = Syy - beta1 * Sxy
    double RSS = Syy - beta1 * Sxy;
    if (RSS < 0)
      RSS = 0.0;

    // sigma^2 estimation (residual variance of the regression)
    double sig2 = RSS / (double)(n - 2);

    // Standard Errors using the passed sigma2_lmm scaling
    // Matches original logic: se = sqrt(sigma2_lmm * diag(inv(X'X)))

    double var_beta1 = sigma2_lmm / Sxx;
    double var_beta0 = sigma2_lmm * (1.0 / n + (mean_x * mean_x) / Sxx);

    coef_mat(col, 0) = beta0;
    coef_mat(col, 1) = beta1;
    se(col, 0) = std::sqrt(var_beta0);
    se(col, 1) = std::sqrt(var_beta1);
    sig(col, 0) = std::sqrt(sig2);

    if (col % 1000 == 0)
      Rcpp::checkUserInterrupt();
  }

  return List::create(Named("Coefficients") = coef_mat, Named("StdErr") = se,
                      Named("sig") = sig);
}
