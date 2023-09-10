by_j_calc_part2_fun <- function(x, Tt_invR, Y_split, fit, intercept_split, W_ast_split, X_split) {
  if(class(fit) %in% c("lmerMod", "lme4") & is.null(X_split)){
    return( Tt_invR[[x]] %*% (Y_split[[x]] - (fixef(fit)[1]*intercept_split[[x]]) - (fixef(fit)[2]*W_ast_split[[x]]) ) )
  } 
  if(class(fit) %in% c("lmerMod", "lme4") & !is.null(X_split)){
    return( Tt_invR[[x]] %*% (Y_split[[x]] - (fixef(fit)[1]*intercept_split[[x]]) - (fixef(fit)[2]*W_ast_split[[x]]) - (fixef(fit)[3]*X_split[[x]]) ) )
  }
  if((!class(fit) %in% c("lmerMod", "lme4")) & is.null(X_split)){
    return( Tt_invR[[x]] %*% (Y_split[[x]] - (fit[1]*intercept_split[[x]]) - (fit[2]*W_ast_split[[x]]) ) )
  } 
  if((!class(fit) %in% c("lmerMod", "lme4")) & !is.null(X_split)){
    return( Tt_invR[[x]] %*% (Y_split[[x]] - (fit[1]*intercept_split[[x]]) - (fit[2]*W_ast_split[[x]]) - (fit[5]*X_split[[x]]) ) ) 
  }
}

b_vec_fun <- function(x, Tt_invR_T, sigma2_lmm_value, b_var_value, by_j_calc_part2) {
  solve( Tt_invR_T[[x]] + sigma2_lmm_value*solve(b_var_value) ) %*% by_j_calc_part2[[x]]
}

b2_vec_fun <- function(x, sigma2_lmm_value, Tt_invR_T, b_var_value, b_vec) {
  solve((sigma2_lmm_value^(-1))*Tt_invR_T[[x]] +  solve(b_var_value)) + b_vec[[x]]%*%t(b_vec[[x]]) 
}

e_fun <- function(x, Y_split, fit, intercept_split, W_ast_split, curr_T_split, b_vec, X_split) {
  if(is.null(X_split)){
    (Y_split[[x]] - (fixef(fit)[1]*intercept_split[[x]]) - (fixef(fit)[2]*W_ast_split[[x]]) - (curr_T_split[[x]]%*%b_vec[[x]])) 
  } else{
    (Y_split[[x]] - (fixef(fit)[1]*intercept_split[[x]]) - (fixef(fit)[2]*W_ast_split[[x]]) - (fixef(fit)[3]*X_split[[x]]) - (curr_T_split[[x]]%*%b_vec[[x]])) 
  }
}

trace_eRe_fun <- function(x, Tt_invR_T, sigma2_lmm_value, b_var_value, et_invR_e) {
  sum(diag(Tt_invR_T[[x]]%*%solve((sigma2_lmm_value^(-1))*Tt_invR_T[[x]] + solve(b_var_value)))) +  et_invR_e[x]
}

Vt_bb_by_j_part_fun <- function(x, sigma2_lmm_value, Tt_invR_T, b_var_value) {
  solve( (sigma2_lmm_value^(-1))*Tt_invR_T[[x]] + solve(b_var_value) )
}

Vt_Wb_by_j_part_fun <- function(x, Tt_invR_T, sigma2_lmm_value, b_var_value, Tt_invR, W_ast_var_split) {
  (solve( Tt_invR_T[[x]] + sigma2_lmm_value*solve(b_var_value) )) %*% Tt_invR[[x]] %*% W_ast_var_split[[x]] 
}

Vt2_2_fun <- function(x, curr_T_split, Vt_bb_by_j_part, curr_T_split_t) {
  curr_T_split[[x]][,1]* Vt_bb_by_j_part[[x]][1,1] * curr_T_split_t[[x]][1,]
} #Var b1

Vt1_2_fun <- function(x, curr_T_split, Vt_Wb_by_j_part) {
  (-1)*curr_T_split[[x]][,1]* Vt_Wb_by_j_part[[x]][1,1]
} #Cov Wb1, Cov b1W

Vt1_3_fun <- function(x, curr_T_split, Vt_Wb_by_j_part) {
  (-1)*curr_T_split[[x]][,2]* Vt_Wb_by_j_part[[x]][2,1]
} #Cov Wb2, Cov b2W

Vt2_3_fun <- function(x, curr_T_split, Vt_bb_by_j_part, curr_T_split_t) {
  curr_T_split[[x]][,1]* Vt_bb_by_j_part[[x]][1,2] * curr_T_split_t[[x]][2,]
} #Cov b1b2

Vt3_2_fun <- function(x, curr_T_split, Vt_bb_by_j_part, curr_T_split_t) {
  curr_T_split[[x]][,2]* Vt_bb_by_j_part[[x]][2,1] * curr_T_split_t[[x]][1,]
} #Cov b2b1

Vt3_3_fun <- function(x, curr_T_split, Vt_bb_by_j_part, curr_T_split_t) {
  curr_T_split[[x]][,2]* Vt_bb_by_j_part[[x]][2,2] * curr_T_split_t[[x]][2,]
} #Var b2

E_b1_squared_fun <- function(x, curr_T_split, Vt_bb_by_j_part, curr_T_split_t, b_vec) {
  curr_T_split[[x]][,1]* Vt_bb_by_j_part[[x]][1,1] * curr_T_split_t[[x]][1,] +
    ((curr_T_split[[x]][,1]*b_vec[[x]][1,1])^2) 
} #Var b1 + Eb1*Eb1

E_b2_squared_fun <- function(x, curr_T_split, Vt_bb_by_j_part, curr_T_split_t, b_vec) {
  curr_T_split[[x]][,2]* Vt_bb_by_j_part[[x]][2,2] * curr_T_split_t[[x]][2,] +
    ((curr_T_split[[x]][,2]*b_vec[[x]][2,1])^2)
} #Var b2 + Eb2*Eb2

E_Wb1_c_fun <- function(x, curr_T_split, Vt_Wb_by_j_part, W_ast_split, b_vec) { 
  (-1)*curr_T_split[[x]][,1]* Vt_Wb_by_j_part[[x]][1,1] +
    ((W_ast_split[[x]])*(curr_T_split[[x]][,1]*b_vec[[x]][1,1]))
} #Cov W-b1 + EW*Eb1

E_Wb2_c_fun <- function(x, curr_T_split, Vt_Wb_by_j_part, W_ast_split, b_vec) { 
  (-1)*curr_T_split[[x]][,2]* Vt_Wb_by_j_part[[x]][2,1] + 
    ((W_ast_split[[x]])*(curr_T_split[[x]][,2]*b_vec[[x]][2,1]))
} #Cov W-b2 + EW*Eb2

E_b1b2_c_fun <- function(x, curr_T_split, Vt_bb_by_j_part, curr_T_split_t, b_vec) {
  curr_T_split[[x]][,1]* Vt_bb_by_j_part[[x]][1,2] * curr_T_split_t[[x]][2,] +
    ((curr_T_split[[x]][,1]*b_vec[[x]][1,1])*(curr_T_split[[x]][,2]*b_vec[[x]][2,1]))
} #Cov b1-b2 + Eb1*Eb2

e_calib_fun <- function(x, Y_split, c_coefs, intercept_split, W_ast_split, curr_T_split, b_vec, X_split) {
  if(is.null(X_split)){
    (Y_split[[x]] - (c_coefs[1]*intercept_split[[x]]) - (c_coefs[2]*W_ast_split[[x]]) - (c_coefs[3]*curr_T_split[[x]][,1]*b_vec[[x]][1,1]) ) 
  } else{
    (Y_split[[x]] - (c_coefs[1]*intercept_split[[x]]) - (c_coefs[2]*W_ast_split[[x]]) - (c_coefs[3]*curr_T_split[[x]][,1]*b_vec[[x]][1,1]) - (c_coefs[4]*curr_T_split[[x]][,2]*b_vec[[x]][2,1]) - (c_coefs[5]*X_split[[x]]) ) 
  }
}

cov_Wb1_c_fun <- function(x, curr_T_split, Vt_Wb_by_j_part) { 
  (-1)*curr_T_split[[x]][,1]* Vt_Wb_by_j_part[[x]][1,1]
} #Cov W-b1 

cov_Wb2_c_fun <- function(x, curr_T_split, Vt_Wb_by_j_part) { 
  (-1)*curr_T_split[[x]][,2]* Vt_Wb_by_j_part[[x]][2,1] 
} #Cov W-b2 

cov_b1b2_c_fun <- function(x, curr_T_split, Vt_bb_by_j_part, curr_T_split_t) {
  curr_T_split[[x]][,1]* Vt_bb_by_j_part[[x]][1,2] * curr_T_split_t[[x]][2,] 
} #Cov b1-b2

