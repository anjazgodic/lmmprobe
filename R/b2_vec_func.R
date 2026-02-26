b2_vec_func <- function(ID, sigma2_lmm, Tt_invR_T, inv_b_var, b_vec, b_vec_int, b_vec_slope) {
  n_subj <- length(unique(ID))

  b2_vec <- Compute_Posterior_Var(
    Tt_invR_T = Tt_invR_T,
    b_vec = b_vec,
    sigma2_lmm = sigma2_lmm,
    inv_b_var = inv_b_var,
    n_subj = n_subj
  )

  b_vec_all <- data.frame(ID = unique(ID))
  if (is.null(b_vec_slope)) {
    b_vec_all$b1 <- b_vec_int
  } else {
    b_vec_all$b1 <- b_vec_int
    b_vec_all$b2 <- b_vec_slope
  }

  # Expansion to full length N
  # b_vec_out should be N x K matrix
  # Original: b_vec_all[rep(seq_len(nrow(b_vec_all)), table(ID)), , drop = FALSE]
  # We need to preserve the ID order if ID is not sorted, but usually it is split by ID.
  # Let's stick to the original expansion logic but be careful.

  # Remove ID for the output matrix
  b_vec_all_no_id <- b_vec_all[, -1, drop = FALSE]
  b_vec_out <- b_vec_all_no_id[rep(seq_len(nrow(b_vec_all_no_id)), as.vector(table(ID))), , drop = FALSE]

  out <- list(b2_vec = b2_vec, b_vec_all = b_vec_all_no_id, b_vec_out = as.matrix(b_vec_out), b2_vec_out = NULL)
  return(out)
}
