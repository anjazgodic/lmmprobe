// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
arma::vec Fast_MVM(const arma::mat &X, const arma::vec &v) {
  int n = X.n_rows;
  int m = X.n_cols;
  arma::vec res(n, arma::fill::zeros);
  for (int j = 0; j < m; ++j) {
    double vj = v[j];
    if (!std::isfinite(vj))
      continue;
    for (int i = 0; i < n; ++i) {
      double xij = X(i, j);
      if (std::isfinite(xij)) {
        res[i] += xij * vj;
      }
    }
  }
  return res;
}
