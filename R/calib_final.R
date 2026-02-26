calib_final <- function(
  W_ast, ID, Y, W_ast_var, sigma2_lmm, Tt_invR_T, b_var, Tt_invR,
  curr_T_split, curr_T_split_t, b_vec, intercept, curr_T, b,
  number_re, X, Vt_bb_by_j_part = NULL, Vt_Wb_by_j_part = NULL,
  expectations = NULL
) {
  ##### Calibration model with current W and current b #####
  W_ast_split <- lapply(split(W_ast, ID), matrix, ncol = 1) # last iteration W from the E-Step
  n_subj <- length(unique(ID))

  # Initialize to NULL to avoid "object not found" errors
  sum_E_b1_squared <- sum_E_Wb1_c <- cov_Wb1_c <- NULL
  sum_E_b2_squared <- sum_E_b1b2_c <- sum_E_Wb2_c <- cov_b1b2_c <- cov_Wb2_c <- NULL

  if (!is.null(expectations)) {
    Vt_bb_by_j_part <- expectations$Vt_bb
    Vt_Wb_by_j_part <- expectations$Vt_Wb

    sum_E_b1_squared <- expectations$sum_E_b1_sq
    sum_E_Wb1_c <- expectations$sum_E_Wb1_c
    cov_Wb1_c <- unlist(expectations$cov_Wb1)

    if (number_re == 2) {
      sum_E_b2_squared <- expectations$sum_E_b2_sq
      sum_E_b1b2_c <- expectations$sum_E_b1b2_c
      sum_E_Wb2_c <- expectations$sum_E_Wb2_c
      cov_b1b2_c <- unlist(expectations$cov_b1b2)
      cov_Wb2_c <- unlist(expectations$cov_Wb2)
    }
  } else {
    if (is.null(Vt_bb_by_j_part) || is.null(Vt_Wb_by_j_part)) {
      W_ast_var_split <- lapply(split(W_ast_var, ID), matrix, ncol = 1) # last iteration W var from the E-Step
      inv_b_var <- solve(b_var)

      if (is.null(Vt_bb_by_j_part)) {
        Vt_bb_by_j_part <- Compute_Posterior_Var(
          Tt_invR_T = Tt_invR_T,
          b_vec = b_vec,
          sigma2_lmm = sigma2_lmm,
          inv_b_var = inv_b_var,
          n_subj = n_subj
        )
      }

      if (is.null(Vt_Wb_by_j_part)) {
        Vt_Wb_by_j_part <- future.apply::future_lapply(1:n_subj,
          FUN = Vt_Wb_by_j_part_fun,
          Tt_invR_T = Tt_invR_T, sigma2_lmm_value = sigma2_lmm,
          inv_b_var_value = inv_b_var, Tt_invR = Tt_invR,
          W_ast_var_split = W_ast_var_split
        )
      }
    }

    sum_E_b1_squared <- sum(unlist(future.apply::future_lapply(1:n_subj,
      FUN = E_b1_squared_fun,
      curr_T_split = curr_T_split, Vt_bb_by_j_part = Vt_bb_by_j_part,
      curr_T_split_t = curr_T_split_t, b_vec = b_vec
    )))

    sum_E_Wb1_c <- sum(unlist(future.apply::future_lapply(1:n_subj,
      FUN = E_Wb1_c_fun,
      curr_T_split = curr_T_split, Vt_Wb_by_j_part = Vt_Wb_by_j_part,
      W_ast_split = W_ast_split, b_vec = b_vec
    )))

    cov_Wb1_c <- unlist(future.apply::future_lapply(1:n_subj,
      FUN = cov_Wb1_c_fun,
      curr_T_split = curr_T_split, Vt_Wb_by_j_part = Vt_Wb_by_j_part
    ))

    if (number_re == 2) {
      sum_E_b2_squared <- sum(unlist(future.apply::future_lapply(1:n_subj,
        FUN = E_b2_squared_fun,
        curr_T_split = curr_T_split, Vt_bb_by_j_part = Vt_bb_by_j_part,
        curr_T_split_t = curr_T_split_t, b_vec = b_vec
      )))

      sum_E_b1b2_c <- sum(unlist(future.apply::future_lapply(1:n_subj,
        FUN = E_b1b2_c_fun,
        curr_T_split = curr_T_split, Vt_bb_by_j_part = Vt_bb_by_j_part,
        curr_T_split_t = curr_T_split_t, b_vec = b_vec
      )))

      sum_E_Wb2_c <- sum(unlist(future.apply::future_lapply(1:n_subj,
        FUN = E_Wb2_c_fun,
        curr_T_split = curr_T_split, Vt_Wb_by_j_part = Vt_Wb_by_j_part,
        W_ast_split = W_ast_split, b_vec = b_vec
      )))

      cov_b1b2_c <- unlist(future.apply::future_lapply(1:n_subj,
        FUN = cov_b1b2_c_fun,
        curr_T_split = curr_T_split, Vt_bb_by_j_part = Vt_bb_by_j_part,
        curr_T_split_t = curr_T_split_t
      ))

      cov_Wb2_c <- unlist(future.apply::future_lapply(1:n_subj,
        FUN = cov_Wb2_c_fun,
        curr_T_split = curr_T_split, Vt_Wb_by_j_part = Vt_Wb_by_j_part
      ))
    }
  }

  if (number_re == 1) {
    Xt_c <- as.matrix(cbind(intercept, W_ast, curr_T[, 1] * b[, 1]))
  } else {
    Xt_c <- as.matrix(cbind(intercept, W_ast, curr_T[, 1] * b[, 1], curr_T[, 2] * b[, 2], X[, "x_adj"]))
  }

  XXt_c <- crossprod(Xt_c)
  XXt_c[2, 2] <- sum(W_ast_var + (W_ast^2))
  XXt_c[2, 3] <- sum_E_Wb1_c
  XXt_c[3, 3] <- sum_E_b1_squared
  XXt_c[3, 2] <- sum_E_Wb1_c

  if (number_re == 2) {
    XXt_c[4, 4] <- sum_E_b2_squared
    XXt_c[2, 4] <- sum_E_Wb2_c
    XXt_c[4, 2] <- sum_E_Wb2_c
    XXt_c[3, 4] <- sum_E_b1b2_c
    XXt_c[4, 3] <- sum_E_b1b2_c
  }

  # Use Cholesky decomposition for better numerical stability and speed
  # Since XXt_c is a crossproduct, it is symmetric and positive-definite
  R_chol <- chol(XXt_c)
  XXt_c_inv <- chol2inv(R_chol)
  c_coefs <- backsolve(R_chol, forwardsolve(t(R_chol), crossprod(Xt_c, Y)))

  out <- list(
    c_coefs = c_coefs,
    Xt_c = Xt_c,
    XXt_c_inv = XXt_c_inv,
    Vt_bb_by_j_part = Vt_bb_by_j_part,
    cov_Wb1_c = cov_Wb1_c,
    cov_b1b2_c = cov_b1b2_c,
    cov_Wb2_c = cov_Wb2_c
  )

  return(out)
}
