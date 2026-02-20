// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <fstream>
#include <iostream>

using namespace Rcpp;

// [[Rcpp::export]]
List UNHIDEM_cpp0_5_w_covs(const arma::vec y, const arma::mat Z,
                           const arma::mat Vt, const arma::colvec X_adj,
                           const arma::uvec rowind1, const arma::uvec rowind2,
                           const arma::uvec rowind3, const arma::uvec colind1,
                           const arma::uvec colind2, const arma::uvec colind3,
                           const arma::colvec Wt, const arma::colvec t_Dt2,
                           const arma::colvec t_Dt3, const arma::colvec Vt1_2,
                           const arma::colvec Vt1_3, const arma::colvec Vt2_3,
                           const arma::colvec W_var, const arma::colvec delta,
                           const arma::colvec beta_vec, double sigma2_lmm) {

  // Getting the dimensions and initializing outputs
  int n = Z.n_rows, d = Z.n_cols;
  arma::mat coef_mat(d, 6, arma::fill::zeros);
  arma::mat se(d, 6, arma::fill::ones);
  arma::mat sig(d, 1, arma::fill::ones);
  arma::mat like(d, 1, arma::fill::ones);

  // Pre-computations (Loop Invariant Terms)
  // bXt columns: [W1, Z, t_Dt1, t_Dt2, t_Dt3, X_adj]
  // Indices: 0, 1, 2, 3, 4, 5

  double sum_y = arma::sum(y);
  double sum_Wt = arma::sum(Wt);
  double sum_t_Dt2 = arma::sum(t_Dt2);
  double sum_t_Dt3 = arma::sum(t_Dt3);
  double sum_X_adj = arma::sum(X_adj);

  double dot_y_y = arma::dot(y, y);
  double dot_Wt_y = arma::dot(Wt, y);
  double dot_t_Dt2_y = arma::dot(t_Dt2, y);
  double dot_t_Dt3_y = arma::dot(t_Dt3, y);
  double dot_X_adj_y = arma::dot(X_adj, y);

  double dot_Wt_Wt = arma::dot(Wt, Wt);
  double dot_Wt_t_Dt2 = arma::dot(Wt, t_Dt2);
  double dot_Wt_t_Dt3 = arma::dot(Wt, t_Dt3);
  double dot_Wt_X_adj = arma::dot(Wt, X_adj);

  double dot_t_Dt2_t_Dt2 = arma::dot(t_Dt2, t_Dt2);
  double dot_t_Dt2_t_Dt3 = arma::dot(t_Dt2, t_Dt3);
  double dot_t_Dt2_X_adj = arma::dot(t_Dt2, X_adj);

  double dot_t_Dt3_t_Dt3 = arma::dot(t_Dt3, t_Dt3);
  double dot_t_Dt3_X_adj = arma::dot(t_Dt3, X_adj);

  double dot_X_adj_X_adj = arma::dot(X_adj, X_adj);

  // Vt sums
  double sum_Vt1_2 = arma::sum(Vt1_2);
  double sum_Vt1_3 = arma::sum(Vt1_3);
  double sum_Vt2_3 = arma::sum(Vt2_3);

  // Vt diagonal parts
  // Vt.submat returns column vector usually if indices are contiguous
  arma::vec Vt_part1 = Vt.submat(rowind1, colind1);
  arma::vec Vt_part2 = Vt.submat(rowind2, colind2);
  arma::vec Vt_part3 = Vt.submat(rowind3, colind3);

  double sum_Vt_part1 = arma::sum(Vt_part1);
  double sum_Vt_part2 = arma::sum(Vt_part2);
  double sum_Vt_part3 = arma::sum(Vt_part3);

  // Constant "correction" terms for bXXt2
  // D2t_2N_2_sum = sum(Vt_part2) + dot(t_Dt2, t_Dt2)
  double D2t_2N_2_sum_val = sum_Vt_part2 + dot_t_Dt2_t_Dt2;
  // D2t_3N_3_sum = sum(Vt_part3) + dot(t_Dt3, t_Dt3)
  double D2t_3N_3_sum_val = sum_Vt_part3 + dot_t_Dt3_t_Dt3;
  // E_b1b2_sum = sum(Vt2_3) + dot(t_Dt2, t_Dt3)
  double E_b1b2_sum_val = sum_Vt2_3 + dot_t_Dt2_t_Dt3;

  // Use fixed size for matrices inside loop to avoid alloc?
  // 6x6 is small enough for stack or fixed.
  // Using arma::mat(6,6) still allocates.
  // But moving it OUT of loop helps? No, need to re-init.
  // Actually, we can reuse one matrix instance and reset elements.
  arma::mat bXXt(6, 6);
  arma::vec bXty(6);
  arma::mat XXt_inv(6, 6);

  for (int col = 0; col < d; ++col) {
    if (col % 1000 == 0)
      Rcpp::checkUserInterrupt();

    double alpha = delta(col) * beta_vec(col);
    const arma::subview_col<double> curr_Z = Z.col(col);

    // Z dots
    double sum_Z = arma::sum(curr_Z);
    double dot_Z_Z = arma::dot(curr_Z, curr_Z);
    double dot_Z_y = arma::dot(curr_Z, y);
    double dot_Z_Wt = arma::dot(curr_Z, Wt);
    double dot_Z_t_Dt2 = arma::dot(curr_Z, t_Dt2);
    double dot_Z_t_Dt3 = arma::dot(curr_Z, t_Dt3);
    double dot_Z_X_adj = arma::dot(curr_Z, X_adj);

    // Construct bXXt
    // Row 0 (W1)
    bXXt(0, 0) = n;
    bXXt(0, 1) = sum_Z;
    bXXt(1, 0) = sum_Z;
    bXXt(0, 2) = sum_Wt - alpha * sum_Z;
    bXXt(2, 0) = bXXt(0, 2);
    bXXt(0, 3) = sum_t_Dt2;
    bXXt(3, 0) = sum_t_Dt2;
    bXXt(0, 4) = sum_t_Dt3;
    bXXt(4, 0) = sum_t_Dt3;
    bXXt(0, 5) = sum_X_adj;
    bXXt(5, 0) = sum_X_adj;

    // Row 1 (Z)
    bXXt(1, 1) = dot_Z_Z;
    bXXt(1, 2) = dot_Z_Wt - alpha * dot_Z_Z;
    bXXt(2, 1) = bXXt(1, 2);
    bXXt(1, 3) = dot_Z_t_Dt2;
    bXXt(3, 1) = bXXt(1, 3);
    bXXt(1, 4) = dot_Z_t_Dt3;
    bXXt(4, 1) = bXXt(1, 4);
    bXXt(1, 5) = dot_Z_X_adj;
    bXXt(5, 1) = bXXt(1, 5);

    // Row 2 (t_Dt1 term)
    // dot(t_Dt1, t_Dt1)
    double dot_t_Dt1_t_Dt1 =
        dot_Wt_Wt - 2 * alpha * dot_Z_Wt + alpha * alpha * dot_Z_Z;
    bXXt(2, 2) = dot_t_Dt1_t_Dt1;
    // dot(t_Dt1, t_Dt2)
    double dot_t_Dt1_t_Dt2_val = dot_Wt_t_Dt2 - alpha * dot_Z_t_Dt2;
    bXXt(2, 3) = dot_t_Dt1_t_Dt2_val;
    bXXt(3, 2) = bXXt(2, 3);
    // dot(t_Dt1, t_Dt3)
    double dot_t_Dt1_t_Dt3_val = dot_Wt_t_Dt3 - alpha * dot_Z_t_Dt3;
    bXXt(2, 4) = dot_t_Dt1_t_Dt3_val;
    bXXt(4, 2) = bXXt(2, 4);
    // dot(t_Dt1, X_adj)
    bXXt(2, 5) = dot_Wt_X_adj - alpha * dot_Z_X_adj;
    bXXt(5, 2) = bXXt(2, 5);

    // Row 3 (t_Dt2)
    bXXt(3, 3) = dot_t_Dt2_t_Dt2;
    bXXt(3, 4) = dot_t_Dt2_t_Dt3;
    bXXt(4, 3) = bXXt(3, 4);
    bXXt(3, 5) = dot_t_Dt2_X_adj;
    bXXt(5, 3) = bXXt(3, 5);

    // Row 4 (t_Dt3)
    bXXt(4, 4) = dot_t_Dt3_t_Dt3;
    bXXt(4, 5) = dot_t_Dt3_X_adj;
    bXXt(5, 4) = bXXt(4, 5);

    // Row 5 (X_adj)
    bXXt(5, 5) = dot_X_adj_X_adj;

    // Use a copy for bXXt2 to preserve original logic?
    // bXXt2 is bXXt with some elements replaced.
    // And we need original bXXt for RSS calculation later?
    // "RSS = y'y - 2*beta'bXty + beta'bXXt*beta"
    // Yes, we need "original" non-corrected bXXt for RSS term 3.
    // So distinct matrix needed.
    arma::mat bXXt2 = bXXt;

    // Apply strict corrections
    // (2,2) -> D2t_1N_1_sum = sum(Vt_part1) + dot(t_Dt1, t_Dt1)
    bXXt2(2, 2) = sum_Vt_part1 + dot_t_Dt1_t_Dt1;

    // (2,3) -> E_b1W_sum = sum(Vt1_2) + dot(t_Dt1, t_Dt2)
    double E_b1W_sum_val = sum_Vt1_2 + dot_t_Dt1_t_Dt2_val;
    bXXt2(2, 3) = E_b1W_sum_val;
    bXXt2(3, 2) = E_b1W_sum_val;

    // (2,4) -> E_b2W_sum = sum(Vt1_3) + dot(t_Dt1, t_Dt3)
    double E_b2W_sum_val = sum_Vt1_3 + dot_t_Dt1_t_Dt3_val;
    bXXt2(2, 4) = E_b2W_sum_val;
    bXXt2(4, 2) = E_b2W_sum_val;

    // (3,3) -> D2t_2N_2_sum_val
    bXXt2(3, 3) = D2t_2N_2_sum_val;

    // (3,4) -> E_b1b2_sum_val
    bXXt2(3, 4) = E_b1b2_sum_val;
    bXXt2(4, 3) = E_b1b2_sum_val;

    // (4,4) -> D2t_3N_3_sum_val
    bXXt2(4, 4) = D2t_3N_3_sum_val;

    if (!bXXt2.is_sympd())
      continue;
    XXt_inv = arma::inv_sympd(bXXt2);

    // bXty
    bXty(0) = sum_y;
    bXty(1) = dot_Z_y;
    bXty(2) = dot_Wt_y - alpha * dot_Z_y;
    bXty(3) = dot_t_Dt2_y; // Corrected: was missing in original logic?
    bXty(4) = dot_t_Dt3_y;
    bXty(5) = dot_X_adj_y;

    arma::vec t_beta = XXt_inv * bXty;

    // RSS calculation using original bXXt
    arma::vec bXXt_beta = bXXt * t_beta;
    double beta_bXXt_beta = arma::dot(t_beta, bXXt_beta);
    double beta_bXty = arma::dot(t_beta, bXty);

    double RSS = dot_y_y - 2.0 * beta_bXty + beta_bXXt_beta;
    if (RSS < 0)
      RSS = 0;

    int k = 6;
    double sigma2_update_val = RSS / (double)(n - k);
    double like_update_val = -n * std::log(sigma2_update_val) / 2.0;

    arma::vec stderrest = arma::sqrt(sigma2_lmm * XXt_inv.diag());

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
