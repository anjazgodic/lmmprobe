// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <omp.h>

using namespace Rcpp;

// [[Rcpp::export]]
List Compute_Random_Effects(const List &Tt_invR, const List &Tt_invR_T,
                            const List &Y_split, const List &W_ast_split,
                            const List &intercept_split,
                            const Nullable<List> &X_split_opt,
                            double sigma2_lmm, const arma::mat &inv_b_var,
                            double coef_int, double coef_W, double coef_X,
                            bool has_X, int n_subj) {

  std::vector<arma::vec> b_vec(n_subj);

  std::vector<arma::mat> Tt_invR_vec(n_subj);
  std::vector<arma::mat> Tt_invR_T_vec(n_subj);
  std::vector<arma::vec> Y_vec(n_subj);
  std::vector<arma::vec> W_vec(n_subj);
  std::vector<arma::vec> int_vec(n_subj);
  std::vector<arma::vec> X_vec;

  if (has_X) {
    X_vec.resize(n_subj);
  }

  for (int i = 0; i < n_subj; ++i) {
    Tt_invR_vec[i] = as<arma::mat>(Tt_invR[i]);
    Tt_invR_T_vec[i] = as<arma::mat>(Tt_invR_T[i]);
    Y_vec[i] = as<arma::vec>(Y_split[i]);
    W_vec[i] = as<arma::vec>(W_ast_split[i]);
    int_vec[i] = as<arma::vec>(intercept_split[i]);
    if (has_X) {
      List X_l = as<List>(X_split_opt);
      X_vec[i] = as<arma::vec>(X_l[i]);
    }
  }

#pragma omp parallel for
  for (int i = 0; i < n_subj; ++i) {
    arma::vec resid;
    if (has_X) {
      resid = Y_vec[i] - (coef_int * int_vec[i]) - (coef_W * W_vec[i]) -
              (coef_X * X_vec[i]);
    } else {
      resid = Y_vec[i] - (coef_int * int_vec[i]) - (coef_W * W_vec[i]);
    }

    arma::vec part2 = Tt_invR_vec[i] * resid;
    arma::mat LHS = Tt_invR_T_vec[i] + sigma2_lmm * inv_b_var;
    b_vec[i] = arma::solve(LHS, part2);
  }

  List result(n_subj);
  for (int i = 0; i < n_subj; ++i) {
    result[i] = b_vec[i];
  }
  return result;
}

// [[Rcpp::export]]
List Compute_Posterior_Var(const List &Tt_invR_T, const List &b_vec,
                           double sigma2_lmm, const arma::mat &inv_b_var,
                           int n_subj) {

  std::vector<arma::mat> b2_vec(n_subj);
  std::vector<arma::mat> Tt_invR_T_vec(n_subj);
  std::vector<arma::vec> b_v(n_subj);

  for (int i = 0; i < n_subj; ++i) {
    Tt_invR_T_vec[i] = as<arma::mat>(Tt_invR_T[i]);
    b_v[i] = as<arma::vec>(b_vec[i]);
  }

#pragma omp parallel for
  for (int i = 0; i < n_subj; ++i) {
    arma::mat LHS = (1.0 / sigma2_lmm) * Tt_invR_T_vec[i] + inv_b_var;
    b2_vec[i] =
        arma::solve(LHS, arma::eye(inv_b_var.n_rows, inv_b_var.n_cols)) +
        (b_v[i] * b_v[i].t());
  }

  List result(n_subj);
  for (int i = 0; i < n_subj; ++i) {
    result[i] = b2_vec[i];
  }
  return result;
}
// [[Rcpp::export]]
List Compute_Expectations(const List &Tt_invR_T, double sigma2_lmm_value,
                          const arma::mat &inv_b_var_value, const List &Tt_invR,
                          const List &W_ast_var_split, const List &curr_T_split,
                          const List &curr_T_split_t, const List &W_ast_split,
                          const List &b_vec, int number_re, int n_subj) {

  std::vector<arma::mat> Tt_invR_T_vec(n_subj);
  std::vector<arma::mat> Tt_invR_vec(n_subj);
  std::vector<arma::mat> W_ast_var_vec(n_subj);
  std::vector<arma::mat> curr_T_vec(n_subj);
  std::vector<arma::mat> curr_T_t_vec(n_subj);
  std::vector<arma::mat> W_ast_vec(n_subj);
  std::vector<arma::mat> b_vec_vec(n_subj);

  for (int i = 0; i < n_subj; ++i) {
    Tt_invR_T_vec[i] = as<arma::mat>(Tt_invR_T[i]);
    Tt_invR_vec[i] = as<arma::mat>(Tt_invR[i]);
    W_ast_var_vec[i] = as<arma::mat>(W_ast_var_split[i]);
    curr_T_vec[i] = as<arma::mat>(curr_T_split[i]);
    curr_T_t_vec[i] = as<arma::mat>(curr_T_split_t[i]);
    W_ast_vec[i] = as<arma::mat>(W_ast_split[i]);
    b_vec_vec[i] = as<arma::mat>(b_vec[i]);
  }

  std::vector<arma::mat> Vt_bb_vec(n_subj);
  std::vector<arma::mat> Vt_Wb_vec(n_subj);

  std::vector<arma::vec> cov_Wb1_vec(n_subj);
  std::vector<arma::vec> cov_Wb2_vec(n_subj);
  std::vector<arma::vec> cov_b1b2_vec(n_subj);

  double sum_E_b1_sq = 0;
  double sum_E_b2_sq = 0;
  double sum_E_b1b2_c = 0;
  double sum_E_Wb1_c = 0;
  double sum_E_Wb2_c = 0;

#pragma omp parallel for reduction(+ : sum_E_b1_sq, sum_E_b2_sq, sum_E_b1b2_c, \
                                       sum_E_Wb1_c, sum_E_Wb2_c)
  for (int i = 0; i < n_subj; ++i) {
    arma::mat A = (1.0 / sigma2_lmm_value) * Tt_invR_T_vec[i] + inv_b_var_value;
    Vt_bb_vec[i] = arma::solve(A, arma::eye(number_re, number_re));

    Vt_Wb_vec[i] = Vt_bb_vec[i] * (Tt_invR_vec[i] * W_ast_var_vec[i]);

    arma::vec T_col0 = curr_T_vec[i].col(0);
    double Vt00 = Vt_bb_vec[i](0, 0);
    double b0 = b_vec_vec[i](0, 0);
    double VtWb00 = Vt_Wb_vec[i](0, 0);

    double term1_b1_sum = arma::dot(T_col0, T_col0) * Vt00;
    double term2_b1_sum = b0 * b0 * arma::dot(T_col0, T_col0);
    sum_E_b1_sq += (term1_b1_sum + term2_b1_sum);

    double term1_Wb1 = -VtWb00 * arma::sum(T_col0);
    double term2_Wb1 = b0 * arma::dot(W_ast_vec[i], T_col0);
    sum_E_Wb1_c += (term1_Wb1 + term2_Wb1);

    cov_Wb1_vec[i] = -T_col0 * VtWb00;

    if (number_re == 2) {
      double Vt11 = Vt_bb_vec[i](1, 1);
      double Vt01 = Vt_bb_vec[i](0, 1);
      arma::vec T_col1 = curr_T_vec[i].col(1);
      double b1 = b_vec_vec[i](1, 0);
      double VtWb10 = Vt_Wb_vec[i](1, 0);

      double term1_b2 = arma::dot(T_col1, T_col1) * Vt11;
      double term2_b2 = b1 * b1 * arma::dot(T_col1, T_col1);
      sum_E_b2_sq += (term1_b2 + term2_b2);

      double term1_Wb2 = -VtWb10 * arma::sum(T_col1);
      double term2_Wb2 = b1 * arma::dot(W_ast_vec[i], T_col1);
      sum_E_Wb2_c += (term1_Wb2 + term2_Wb2);

      double term1_b1b2 = Vt01 * arma::dot(T_col0, T_col1);
      double term2_b1b2 = (b0 * b1) * arma::dot(T_col0, T_col1);
      sum_E_b1b2_c += (term1_b1b2 + term2_b1b2);

      cov_Wb2_vec[i] = -T_col1 * VtWb10;
      cov_b1b2_vec[i] = (T_col0 * Vt01) % T_col1;
    }
  }

  List result;
  result["Vt_bb"] = Vt_bb_vec;
  result["Vt_Wb"] = Vt_Wb_vec;
  result["sum_E_b1_sq"] = sum_E_b1_sq;
  result["sum_E_Wb1_c"] = sum_E_Wb1_c;
  result["cov_Wb1"] = cov_Wb1_vec;

  if (number_re == 2) {
    result["sum_E_b2_sq"] = sum_E_b2_sq;
    result["sum_E_b1b2_c"] = sum_E_b1b2_c;
    result["sum_E_Wb2_c"] = sum_E_Wb2_c;
    result["cov_Wb2"] = cov_Wb2_vec;
    result["cov_b1b2"] = cov_b1b2_vec;
  }

  return result;
}

// [[Rcpp::export]]
NumericVector Compute_Trace(const List &Tt_invR_T, double sigma2_lmm_value,
                            const arma::mat &inv_b_var_value,
                            const List &et_invR_e, int number_re, int n_subj) {

  NumericVector result(n_subj);

  std::vector<arma::mat> Tt_invR_T_vec(n_subj);
  std::vector<double> et_invR_e_vec(n_subj);

  for (int i = 0; i < n_subj; ++i) {
    Tt_invR_T_vec[i] = as<arma::mat>(Tt_invR_T[i]);
    et_invR_e_vec[i] = as<double>(et_invR_e[i]);
  }

#pragma omp parallel for
  for (int i = 0; i < n_subj; ++i) {
    arma::mat A = (1.0 / sigma2_lmm_value) * Tt_invR_T_vec[i] + inv_b_var_value;
    arma::mat X_mat = arma::solve(A, Tt_invR_T_vec[i]);
    double tr = arma::trace(X_mat);
    result[i] = tr + et_invR_e_vec[i];
  }

  return result;
}

// [[Rcpp::export]]
arma::mat Compute_Full_Vt(const List &Vt_bb_vec, const List &Vt_Wb_vec,
                          const List &curr_T_split, const arma::vec &W_ast_var,
                          int number_re, int n_total) {

  int n_subj = Vt_bb_vec.size();
  std::vector<arma::mat> Vbb(n_subj);
  std::vector<arma::mat> VWb(n_subj);
  std::vector<arma::mat> T(n_subj);

  std::vector<int> subj_indices(n_subj + 1, 0);

  for (int i = 0; i < n_subj; ++i) {
    Vbb[i] = as<arma::mat>(Vt_bb_vec[i]);
    VWb[i] = as<arma::mat>(Vt_Wb_vec[i]);
    T[i] = as<arma::mat>(curr_T_split[i]);
    subj_indices[i + 1] = subj_indices[i] + T[i].n_rows;
  }

  // Pre-allocate output vectors
  arma::vec Vt1_1 = W_ast_var;
  arma::vec Vt2_2(n_total);
  arma::vec Vt1_2(n_total);

  arma::vec Vt3_3, Vt1_3, Vt2_3;
  if (number_re == 2) {
    Vt3_3.set_size(n_total);
    Vt1_3.set_size(n_total);
    Vt2_3.set_size(n_total);
  }

#pragma omp parallel for
  for (int i = 0; i < n_subj; ++i) {
    int start = subj_indices[i];
    int n_i = T[i].n_rows;
    if (n_i == 0)
      continue;

    arma::vec T0 = T[i].col(0);
    double Vbb00 = Vbb[i](0, 0);
    double VWb00 = VWb[i](0, 0);

    for (int k = 0; k < n_i; ++k) {
      double t0 = T0[k];
      Vt2_2[start + k] = t0 * t0 * Vbb00;
      Vt1_2[start + k] = -t0 * VWb00;
    }

    if (number_re == 2) {
      arma::vec T1 = T[i].col(1);
      double Vbb11 = Vbb[i](1, 1);
      double Vbb01 = Vbb[i](0, 1);
      double VWb10 = VWb[i](1, 0);

      for (int k = 0; k < n_i; ++k) {
        double t1 = T1[k];
        double t0 = T0[k];
        Vt3_3[start + k] = t1 * t1 * Vbb11;
        Vt1_3[start + k] = -t1 * VWb10;
        Vt2_3[start + k] = t0 * t1 * Vbb01;
      }
    }
  }

  int N = n_total;

  if (number_re == 1) {
    arma::mat Vt(2 * N, 2);
    // Col 1
    Vt.submat(0, 0, N - 1, 0) = Vt1_1;
    Vt.submat(N, 0, 2 * N - 1, 0) = Vt1_2;
    // Col 2
    Vt.submat(0, 1, N - 1, 1) = Vt1_2;
    Vt.submat(N, 1, 2 * N - 1, 1) = Vt2_2;
    return Vt;
  } else {
    arma::mat Vt(3 * N, 3);
    // Col 1
    Vt.submat(0, 0, N - 1, 0) = Vt1_1;
    Vt.submat(N, 0, 2 * N - 1, 0) = Vt1_2;
    Vt.submat(2 * N, 0, 3 * N - 1, 0) = Vt1_3;

    // Col 2
    Vt.submat(0, 1, N - 1, 1) = Vt1_2;
    Vt.submat(N, 1, 2 * N - 1, 1) = Vt2_2;
    Vt.submat(2 * N, 1, 3 * N - 1, 1) = Vt2_3;

    // Col 3
    Vt.submat(0, 2, N - 1, 2) = Vt1_3;
    Vt.submat(N, 2, 2 * N - 1, 2) = Vt2_3;
    Vt.submat(2 * N, 2, 3 * N - 1, 2) = Vt3_3;
    return Vt;
  }
}
