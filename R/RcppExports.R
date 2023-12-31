# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

Col_sum <- function(X) {
    .Call(`_lmmprobe_Col_sum`, X)
}

LM_by_col <- function(y, X, sigma2_lmm) {
    .Call(`_lmmprobe_LM_by_col`, y, X, sigma2_lmm)
}

LM_by_col_w_covs <- function(y, X, X_adj, sigma2_lmm) {
    .Call(`_lmmprobe_LM_by_col_w_covs`, y, X, X_adj, sigma2_lmm)
}

MVM <- function(X, v) {
    .Call(`_lmmprobe_MVM`, X, v)
}

Row_sum <- function(X) {
    .Call(`_lmmprobe_Row_sum`, X)
}

UNHIDEM_cpp0_5 <- function(y, Z, Vt, rowind1, rowind2, colind1, colind2, Wt, t_Dt2, Vt1_2, W_var, delta, beta_vec, sigma2_lmm) {
    .Call(`_lmmprobe_UNHIDEM_cpp0_5`, y, Z, Vt, rowind1, rowind2, colind1, colind2, Wt, t_Dt2, Vt1_2, W_var, delta, beta_vec, sigma2_lmm)
}

UNHIDEM_cpp0_5_w_covs <- function(y, Z, Vt, X_adj, rowind1, rowind2, rowind3, colind1, colind2, colind3, Wt, t_Dt2, t_Dt3, Vt1_2, Vt1_3, Vt2_3, W_var, delta, beta_vec, sigma2_lmm) {
    .Call(`_lmmprobe_UNHIDEM_cpp0_5_w_covs`, y, Z, Vt, X_adj, rowind1, rowind2, rowind3, colind1, colind2, colind3, Wt, t_Dt2, t_Dt3, Vt1_2, Vt1_3, Vt2_3, W_var, delta, beta_vec, sigma2_lmm)
}

