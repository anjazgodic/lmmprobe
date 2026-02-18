b2_vec_func <- function(ID, sigma2_lmm, Tt_invR_T, inv_b_var, b_vec, b_vec_int, b_vec_slope) {
  b2_vec <- future.apply::future_lapply(1:length(unique(ID)), FUN = b2_vec_fun, sigma2_lmm_value = sigma2_lmm, Tt_invR_T = Tt_invR_T, inv_b_var_value = inv_b_var, b_vec = b_vec)

  b_vec_all <- data.frame(cbind(b_vec_int, b_vec_slope))
  b_vec_out <- b_vec_all[rep(seq_len(nrow(b_vec_all)), table(ID)), , drop = FALSE]

  out <- list(b2_vec = b2_vec, b_vec_all = b_vec_all, b_vec_out = b_vec_out, b2_vec_out = NULL)
  return(out)
}
