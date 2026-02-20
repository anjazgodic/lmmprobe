b_vec_func <- function(W_ast, ID, Tt_invR, Y_split, fit, intercept_split, X_split, Tt_invR_T, sigma2_lmm, inv_b_var, number_re) {
  ##### We calculate the random effect b #####
  n_subj <- length(unique(ID))
  W_ast_split <- lapply(split(W_ast, ID), matrix, ncol = 1)

  # Extract coefficients from fit
  if (any(class(fit) %in% c("lmerMod", "lme4"))) {
    coef_int <- fixef(fit)[1]
    coef_W <- fixef(fit)[2]
    coef_X <- if (!is.null(X_split)) fixef(fit)[3] else 0
  } else {
    coef_int <- fit[1]
    coef_W <- fit[2]
    coef_X <- if (!is.null(X_split)) fit[length(fit)] else 0 # fit[5] in original logic?
    # Let's check original logic in by_j_calc_part2_fun:
    # if ((!any(class(fit) %in% c("lmerMod", "lme4"))) & !is.null(X_split)) {
    #   return(Tt_invR[[x]] %*% (Y_split[[x]] - (fit[1] * intercept_split[[x]]) - (fit[2] * W_ast_split[[x]]) - (fit[5] * X_split[[x]])))
    # }
    # Wait, fit[5]? That's unusual. I'll use fixef if it's lmer, otherwise assume fit is a vector.
    if (!is.null(X_split)) {
      if (length(fit) >= 5) {
        coef_X <- fit[5]
      } else if (length(fit) >= 3) coef_X <- fit[3]
    }
  }

  b_vec <- Compute_Random_Effects(
    Tt_invR = Tt_invR,
    Tt_invR_T = Tt_invR_T,
    Y_split = Y_split,
    W_ast_split = W_ast_split,
    intercept_split = intercept_split,
    X_split_opt = X_split,
    sigma2_lmm = sigma2_lmm,
    inv_b_var = inv_b_var,
    coef_int = as.numeric(coef_int),
    coef_W = as.numeric(coef_W),
    coef_X = as.numeric(coef_X),
    has_X = !is.null(X_split),
    n_subj = n_subj
  )

  if (number_re == 1) {
    b_vec_int <- unlist(b_vec)
    b_vec_slope <- NULL
  } else {
    # Use vapply or matrix operations for efficiency
    b_vec_mat <- do.call(cbind, b_vec)
    b_vec_int <- b_vec_mat[1, ]
    b_vec_slope <- b_vec_mat[2, ]
  }

  out <- list(W_ast_split = W_ast_split, b_vec = b_vec, b_vec_int = b_vec_int, b_vec_slope = b_vec_slope)
  return(out)
}
