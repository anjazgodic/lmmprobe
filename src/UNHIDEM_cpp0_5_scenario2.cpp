// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <fstream>
#include <iostream>

using namespace Rcpp;

// [[Rcpp::export]]
List UNHIDEM_cpp0_5(const arma::vec y, const arma::mat Z, const arma::mat Vt,
                    const arma::uvec rowind1, const arma::uvec rowind2,
                    const arma::uvec colind1, const arma::uvec colind2,
                    const arma::colvec Wt, const arma::colvec t_Dt2,
                    const arma::colvec Vt1_2, const arma::colvec W_var,
                    const arma::colvec delta, const arma::colvec beta_vec,
                    double sigma2_lmm) {

  // Getting the dimensions and initializing outputs
  int n = Z.n_rows, d = Z.n_cols;
  arma::mat coef_mat(d, 4, arma::fill::zeros);
  arma::mat se(d, 4, arma::fill::ones);
  arma::mat sig(d, 1, arma::fill::ones);
  arma::mat like(d, 1, arma::fill::ones);

  // Pre-computations (Loop Invariant Terms)
  // W1 is a vector of ones.
  // We need to compute dot products for the 4x4 matrix bXXt and 4x1 vector
  // bXty efficiently.
  // bXt columns: [W1, Z_col, t_Dt1, t_Dt2]
  // t_Dt1 = Wt - alpha * Z_col (where alpha = delta * beta)

  double sum_y = arma::sum(y);
  double sum_Wt = arma::sum(Wt);
  double sum_t_Dt2 = arma::sum(t_Dt2); // dot(W1, t_Dt2)

  double dot_y_y = arma::dot(y, y);
  double dot_Wt_y = arma::dot(Wt, y);
  double dot_t_Dt2_y = arma::dot(t_Dt2, y);

  double dot_Wt_Wt = arma::dot(Wt, Wt);
  double dot_Wt_t_Dt2 = arma::dot(Wt, t_Dt2);
  double dot_t_Dt2_t_Dt2 = arma::dot(t_Dt2, t_Dt2);

  double sum_Vt1_2 = arma::sum(Vt1_2); // Used for E_b1W_sum
  // D2t terms precomputation
  // D2t_1N_1 = Vt(rowind1, colind1) + t_Dt1^2
  // D2t_2N_2 = Vt(rowind2, colind2) + t_Dt2^2
  // D2t_2N_2_sum = sum(Vt_part2) + dot(t_Dt2, t_Dt2)
  // Vt_part2 is Vt.submat(rowind2, colind2) which corresponds to [N..2N-1, 1]
  // Since rowind2/colind2 are constant, this sum is constant.
  // Note: Vt is passed as matrix. rowind2 are indices.
  // Assuming rowind2 are contiguous N..2N-1 and colind2 is scalar 1.
  arma::vec Vt_part1 = Vt.submat(rowind1, colind1); // For D2t_1N_1
  arma::vec Vt_part2 = Vt.submat(rowind2, colind2); // For D2t_2N_2

  double sum_Vt_part1 = arma::sum(Vt_part1); // Effectively sum(W_var)
  double sum_Vt_part2 = arma::sum(Vt_part2);

  double D2t_2N_2_sum_val = sum_Vt_part2 + dot_t_Dt2_t_Dt2;

  // Temporary vars for the loop
  arma::mat bXXt(4, 4);
  arma::vec bXty(4);
  arma::mat XXt_inv(4, 4);

  // Main Loop
  for (int col = 0; col < d; ++col) {
    if (col % 1000 == 0)
      Rcpp::checkUserInterrupt();

    double alpha = delta(col) * beta_vec(col);

    // Access Z column without copy
    const arma::subview_col<double> curr_Z = Z.col(col);

    // Compute Z-dependent dot products
    double sum_Z = arma::sum(curr_Z);
    double dot_Z_Z = arma::dot(curr_Z, curr_Z);
    double dot_Z_y = arma::dot(curr_Z, y);
    double dot_Z_Wt = arma::dot(curr_Z, Wt);
    double dot_Z_t_Dt2 = arma::dot(curr_Z, t_Dt2);

    // Construct bXXt manualy
    // [0] W1, [1] Z, [2] t_Dt1, [3] t_Dt2
    // t_Dt1 = Wt - alpha * Z

    // (0,0) sum(1) = n
    bXXt(0, 0) = n;
    // (0,1) sum(Z)
    bXXt(0, 1) = sum_Z;
    bXXt(1, 0) = sum_Z;
    // (0,2) sum(t_Dt1) = sum(Wt) - alpha*sum(Z)
    bXXt(0, 2) = sum_Wt - alpha * sum_Z;
    bXXt(2, 0) = bXXt(0, 2);
    // (0,3) sum(t_Dt2)
    bXXt(0, 3) = sum_t_Dt2;
    bXXt(3, 0) = sum_t_Dt2;

    // (1,1) dot(Z,Z)
    bXXt(1, 1) = dot_Z_Z;
    // (1,2) dot(Z, t_Dt1) = dot(Z, Wt) - alpha*dot(Z,Z)
    bXXt(1, 2) = dot_Z_Wt - alpha * dot_Z_Z;
    bXXt(2, 1) = bXXt(1, 2);
    // (1,3) dot(Z, t_Dt2)
    bXXt(1, 3) = dot_Z_t_Dt2;
    bXXt(3, 1) = bXXt(1, 3);

    // (2,2) dot(t_Dt1, t_Dt1)
    // = dot(Wt - aZ, Wt - aZ) = dot(Wt,Wt) - 2a dot(Wt,Z) + a^2 dot(Z,Z)
    double dot_t_Dt1_t_Dt1 =
        dot_Wt_Wt - 2 * alpha * dot_Z_Wt + alpha * alpha * dot_Z_Z;
    bXXt(2, 2) = dot_t_Dt1_t_Dt1;
    // (2,3) dot(t_Dt1, t_Dt2) = dot(Wt, t_Dt2) - alpha*dot(Z, t_Dt2)
    bXXt(2, 3) = dot_Wt_t_Dt2 - alpha * dot_Z_t_Dt2;
    bXXt(3, 2) = bXXt(2, 3);

    // (3,3) dot(t_Dt2, t_Dt2)
    bXXt(3, 3) = dot_t_Dt2_t_Dt2;

    // Construct bXXt2 (Corrected Version for Inverse)
    arma::mat bXXt2 = bXXt; // Copy 4x4
    // Modifications from original code:
    // bXXt2(2,2) = D2t_1N_1_sum(0);
    // D2t_1N_1_sum = sum(Vt_part1) + dot(t_Dt1, t_Dt1)
    bXXt2(2, 2) = sum_Vt_part1 + dot_t_Dt1_t_Dt1;

    // bXXt2(2,3) = E_b1W_sum(0);
    // E_b1W = Vt1_2 + t_Dt1*t_Dt2
    // Sum = sum(Vt1_2) + dot(t_Dt1, t_Dt2)
    double E_b1W_sum_val =
        sum_Vt1_2 + bXXt(2, 3); // bXXt(2,3) is dot(t_Dt1, t_Dt2)
    bXXt2(2, 3) = E_b1W_sum_val;
    bXXt2(3, 2) = E_b1W_sum_val;

    // bXXt2(3,3) = D2t_2N_2_sum(0);
    bXXt2(3, 3) = D2t_2N_2_sum_val;

    // Check singularity and invert
    // Since 4x4, we could unroll or use armadillo.
    // For safety and logic preservation, access inv_sympd.
    if (!bXXt2.is_sympd())
      continue;
    XXt_inv = arma::inv_sympd(bXXt2);

    // Compute bXty manually
    // [0] sum(y)
    bXty(0) = sum_y;
    // [1] dot(Z, y)
    bXty(1) = dot_Z_y;
    // [2] dot(t_Dt1, y) = dot(Wt, y) - alpha*dot(Z, y)
    bXty(2) = dot_Wt_y - alpha * dot_Z_y;
    // [3] dot(t_Dt2, y)
    bXty(3) = dot_t_Dt2_y;

    // t_beta = XXt_inv * bXty
    arma::vec t_beta = XXt_inv * bXty;

    // sigma^2 update
    // resid'resid = (y - X*beta)'(y - X*beta)
    // = y'y - 2*beta'X'y + beta'X'X*beta
    // Note: X'X here is the ORIGINAL X'X (bXXt), NOT bXXt2.
    // Calculations:
    // term1 = y'y
    // term2 = 2 * dot(t_beta, bXty)  (Since X'y is bXty)
    // term3 = t_beta' * bXXt * t_beta
    // Optimization: bXXt * t_beta can be computed 4x4 * 4x1 -> 4x1.

    // Explicitly:
    // RSS = y'y - 2*y'X*beta + beta'X'X*beta
    // RSS = y'y - 2*beta'*bXty + beta'*(bXXt*beta)

    arma::vec bXXt_beta = bXXt * t_beta;
    double beta_bXXt_beta = arma::dot(t_beta, bXXt_beta);
    double beta_bXty = arma::dot(t_beta, bXty);

    double RSS = dot_y_y - 2.0 * beta_bXty + beta_bXXt_beta;
    if (RSS < 0)
      RSS = 0; // Numerical safety

    int k = 4;
    double sigma2_update_val = RSS / (double)(n - k);
    double like_update_val = -n * std::log(sigma2_update_val) / 2.0;

    // SEs
    // Vbt = sigma2_lmm * diag( inv(bXXt2) * bXXt2 * inv(bXXt2) )
    // Note: XXt_inv is inv(bXXt2).
    // So XXt_inv * bXXt2 * XXt_inv = XXt_inv.
    // So Vbt is simply sigma2_lmm * diagvec(XXt_inv).

    arma::vec stderrest = arma::sqrt(sigma2_lmm * XXt_inv.diag());

    // Store results
    coef_mat.row(col) = t_beta.t();
    se.row(col) = stderrest.t();
    sig(col, 0) = sigma2_update_val;
    like(col, 0) = like_update_val;
  }

  // Test statistics
  arma::mat T_vals = coef_mat / se;

  return List::create(Named("Coefficients") = coef_mat, Named("StdErr") = se,
                      Named("Sigma") = sig, Named("T_statistics") = T_vals,
                      Named("Log_like") = like);
}
