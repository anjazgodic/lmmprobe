library(Rcpp)

# Define Clean versions of time measurement
time_call <- function(expression_to_run) {
    t <- system.time(expression_to_run)
    return(t["elapsed"])
}

# ---------------------------------------------------------
# PROOF OF CONCEPT: Uninitialized Memory in Armadillo
# ---------------------------------------------------------
sourceCpp(code = "
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat uninit_test(int n) {
  arma::mat m(n, 4);
  // Function returns uninitialized memory
  // In UNHIDEM, if is_sympd is false, this is what is returned
  return m;
}
")

cat("--------------------------------------------------------\n")
cat("TEST 1: Uninitialized Memory Check (Potential Correctness Bug)\n")
cat("--------------------------------------------------------\n")
m <- uninit_test(5)
print(m)
cat("If the above matrix contains random numbers (not zeros), it confirms the bug.\n\n")

# ---------------------------------------------------------
# BENCHMARK: LM_by_col optimization
# ---------------------------------------------------------
sourceCpp(code = '
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
List LM_by_col_old(const arma::vec y, const arma::mat X, double sigma2_lmm) {
  int n = X.n_rows, d = X.n_cols;
  arma::mat coef_mat(d,2);
  arma::mat se(d,2);
  arma::mat sig(d,1);
  arma::mat X1(n,1,arma::fill::ones);

  for(int col=0;col<d; ++col){
    arma::mat X2=arma::join_rows(X1,X.col(col));
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

// [[Rcpp::export]]
List LM_by_col_new(const arma::vec& y, const arma::mat& X, double sigma2_lmm) {

  int n = X.n_rows, d = X.n_cols;

  arma::mat coef_mat(d,2);
  coef_mat.fill(NA_REAL); // SAFETY: Initialize with NA
  arma::mat se(d,2);
  se.fill(NA_REAL);
  arma::mat sig(d,1);
  sig.fill(NA_REAL);

  // Precompute centered Y
  double mean_y = arma::mean(y);
  arma::vec y_cen = y - mean_y;
  double Syy = arma::dot(y_cen, y_cen);

  for(int col=0;col<d; ++col){

    const arma::subview_col<double> x = X.col(col);
    double mean_x = arma::mean(x);
    arma::vec x_cen = x - mean_x;

    double Sxx = arma::dot(x_cen, x_cen);

    if (Sxx < 1e-12) {
       continue;
    }

    double Sxy = arma::dot(x_cen, y_cen);
    double beta1 = Sxy / Sxx;
    double beta0 = mean_y - beta1 * mean_x;

    double RSS = Syy - beta1 * Sxy;
    if (RSS < 0) RSS = 0.0;

    double sig2 = RSS / (double)(n - 2);

    double var_beta1 = sigma2_lmm / Sxx;
    double var_beta0 = sigma2_lmm * (1.0/n + (mean_x*mean_x)/Sxx);

    coef_mat(col,0) = beta0;
    coef_mat(col,1) = beta1;
    se(col,0) = std::sqrt(var_beta0);
    se(col,1) = std::sqrt(var_beta1);
    sig(col,0) = std::sqrt(sig2);

    if (col % 1000 == 0) Rcpp::checkUserInterrupt();
  }

  return List::create(Named("Coefficients") = coef_mat,Named("StdErr") = se,Named("sig")=sig);
}

')

cat("--------------------------------------------------------\n")
cat("TEST 2: Benchmark LM optimization\n")
cat("--------------------------------------------------------\n")
set.seed(42)
n <- 2000
d <- 2000
cat(sprintf("Generating %dx%d matrix...\n", n, d))
X <- matrix(rnorm(n * d), n, d)
y <- rnorm(n)
sigma2_lmm <- 1.0

cat("Running Old...\n")
# Run once to warm up?
invisible(LM_by_col_old(y[1:100], X[1:100, 1:10], sigma2_lmm))

t_old <- system.time(LM_by_col_old(y, X, sigma2_lmm))["elapsed"]
cat(sprintf("Old Time: %.4f s\n", t_old))

cat("Running New...\n")
invisible(LM_by_col_new(y[1:100], X[1:100, 1:10], sigma2_lmm))

t_new <- system.time(LM_by_col_new(y, X, sigma2_lmm))["elapsed"]
cat(sprintf("New Time: %.4f s\n", t_new))

cat(sprintf("Speedup: %.2fx\n", t_old / t_new))

res_old <- LM_by_col_old(y, X, sigma2_lmm)
res_new <- LM_by_col_new(y, X, sigma2_lmm)
diff <- max(abs(res_old$Coefficients - res_new$Coefficients), na.rm = TRUE)

# ---------------------------------------------------------
# BENCHMARK: LM_by_col_w_covs optimization
# ---------------------------------------------------------
sourceCpp(code = '
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
List LM_by_col_w_covs_old(const arma::vec y, const arma::mat X, const arma::colvec X_adj, double sigma2_lmm) {
  int n = X.n_rows, d = X.n_cols;
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

// [[Rcpp::export]]
List LM_by_col_w_covs_new(const arma::vec& y, const arma::mat& X, const arma::colvec& X_adj, double sigma2_lmm) {
  int n = X.n_rows, d = X.n_cols;
  arma::mat coef_mat(d,3);
  arma::mat se(d,3);
  arma::mat sig(d,1);
  
  // Pre-allocate design matrix X2 (N x 3)
  // [1, X_col, X_adj]
  arma::mat X2(n, 3);
  X2.col(0).ones();
  X2.col(2) = X_adj;
  
  for(int col=0;col<d; ++col){
    // Copy column into X2
    X2.col(1) = X.col(col);
    
    int k = 3;
    
    // Check for singularity? solve() handles it (throws) or we can use solve(..., opts::fast)
    // To match old behavior but be faster, just call solve.
    // Ideally we catch exceptions, but old code didn't.
    
    arma::colvec coef = arma::solve(X2, y); 
    arma::colvec resid = y - X2*coef; 
    
    double sig2 = arma::dot(resid, resid) / (n - k);
    
    // Standard error
    // inv(X^T X)
    // Efficiently: inv_sympd(X2.t() * X2) is 3x3.
    // X2.t()*X2 calculation involved N*k^2 ops.
    // 3*3 dot products.
    arma::mat XtX = arma::trans(X2) * X2;
    arma::vec diag_inv = arma::diagvec(arma::inv(XtX)); // 3x3 inv is cheap
    
    arma::colvec stderrest = arma::sqrt(sigma2_lmm * diag_inv);
    double s = std::sqrt(sig2);
    
    coef_mat.row(col) = arma::trans(coef);
    se.row(col) = arma::trans(stderrest);
    sig.row(col) = s; 
    
    if (col % 1000 == 0) Rcpp::checkUserInterrupt();
  }
  return List::create(Named("Coefficients") = coef_mat,Named("StdErr") = se,Named("sig")=sig);
}
')

cat("--------------------------------------------------------\n")
cat("TEST 3: Benchmark LM w/ Covs optimization\n")
cat("--------------------------------------------------------\n")
X_adj <- rnorm(n)
cat("Running Old w/ Covs...\n")
t_old_c <- system.time(LM_by_col_w_covs_old(y, X, X_adj, sigma2_lmm))["elapsed"]
cat(sprintf("Old Time: %.4f s\n", t_old_c))

cat("Running New w/ Covs...\n")
t_new_c <- system.time(LM_by_col_w_covs_new(y, X, X_adj, sigma2_lmm))["elapsed"]
cat(sprintf("New Time: %.4f s\n", t_new_c))
cat(sprintf("Speedup: %.2fx\n", t_old_c/t_new_c))

