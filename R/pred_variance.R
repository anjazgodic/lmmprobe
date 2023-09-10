pred_var <- function(W_ast, ID, W_ast_var, VCV, 
                  curr_T_split, curr_T_split_t, 
                  curr_T, curr_T2, b, 
                  number_re, c_coefs, 
                  sigma2_lmm, Tt_invR_T, b_var, 
                  Vt_bb_by_j_part, 
                  cov_Wb1_c, cov_b1b2_c, cov_Wb2_c){
  
  if(number_re == 1){
    V_w   <- W_ast_var+(W_ast^2)
    V_b   <- unlist(sfLapply(1:length(unique(ID)), fun=Vt2_2_fun, curr_T_split, Vt_bb_by_j_part, curr_T_split_t)) #Var b1
    V_wb  <- cov_Wb1_c
    S2_mu <- VCV[1,1]
    S2_alpha <- VCV[2,2]
    S2_gamma <- VCV[3,3]
    S2_mu_alpha <- VCV[1,2]
    S2_mu_gamma <- VCV[1,3]
    S2_alpha_gamma <- VCV[2,3]
    mu_hat <- c_coefs[1]
    alpha_hat <- c_coefs[2]
    gamma_hat <- c_coefs[3]
    
    var_Yhat_new_train <- S2_mu + S2_alpha*V_w + (alpha_hat^2)*V_w + (W_ast^2)*S2_alpha +
      curr_T2*S2_gamma*V_b + curr_T2*(gamma_hat^2)*V_b + 
      curr_T2*(b^2)*S2_gamma - 2*W_ast*S2_mu_alpha - 2*curr_T*b*S2_mu_gamma -
      2*curr_T*S2_alpha_gamma*V_wb - 2*curr_T*V_wb*alpha_hat*gamma_hat - 
      2*curr_T*S2_alpha_gamma*W_ast*b
  } else{
    V_w   <- W_ast_var+(W_ast^2)
    V_b1   <- unlist(sfLapply(1:length(unique(ID)), fun=Vt2_2_fun, curr_T_split, Vt_bb_by_j_part, curr_T_split_t)) #Var b1
    V_b2   <- unlist(sfLapply(1:length(unique(ID)), fun=Vt3_3_fun, curr_T_split, Vt_bb_by_j_part, curr_T_split_t)) #Var b2
    V_b1b2   <- cov_b1b2_c
    V_wb1  <- cov_Wb1_c
    V_wb2  <- cov_Wb2_c
    S2_mu <- VCV[1,1]
    S2_alpha <- VCV[2,2]
    S2_gamma <- VCV[3,3]
    S2_omega <- VCV[4,4]
    S2_X_adj <- VCV[5,5]
    S2_mu_alpha <- VCV[1,2]
    S2_mu_gamma <- VCV[1,3]
    S2_mu_omega <- VCV[1,4]
    S2_mu_X_adj <- VCV[1,5]
    S2_alpha_gamma <- VCV[2,3]
    S2_alpha_omega <- VCV[2,4]
    S2_alpha_X_adj <- VCV[2,5]
    S2_gamma_omega <- VCV[3,4]
    S2_gamma_X_adj <- VCV[3,5]
    S2_omega_X_adj <- VCV[4,5]
    mu_hat <- c_coefs[1]
    alpha_hat <- c_coefs[2]
    gamma_hat <- c_coefs[3]
    omega_hat <- c_coefs[4]
    
    var_Yhat_new_train <- 
      S2_mu + 
      S2_alpha*V_w + 
      (alpha_hat^2)*V_w + 
      (W_ast^2)*S2_alpha +
      curr_T2[,1]*S2_gamma*V_b1 + 
      curr_T2[,1]*(gamma_hat^2)*V_b1 + 
      curr_T2[,1]*(b[,1]^2)*S2_gamma +
      curr_T2[,2]*S2_omega*V_b2 +
      curr_T2[,2]*(omega_hat^2)*V_b2 + 
      curr_T2[,2]*(b[,2]^2)*S2_omega - 
      2*W_ast*S2_mu_alpha - 
      2*curr_T[,1]*b[,1]*S2_mu_gamma -
      2*curr_T[,2]*b[,2]*S2_mu_omega - 
      2*curr_T[,1]*S2_alpha_gamma*V_wb1 - 
      2*curr_T[,1]*V_wb1*alpha_hat*gamma_hat - 
      2*curr_T[,1]*S2_alpha_gamma*W_ast*b[,1] - 
      2*curr_T[,2]*S2_alpha_omega*V_wb2 - 
      2*curr_T[,2]*V_wb2*alpha_hat*omega_hat - 
      2*curr_T[,2]*S2_alpha_omega*W_ast*b[,2] - 
      2*curr_T[,1]*curr_T[,2]*S2_gamma_omega*V_b1b2 - 
      2*curr_T[,1]*curr_T[,2]*V_b1b2*gamma_hat*omega_hat - 
      2*curr_T[,1]*curr_T[,2]*S2_gamma_omega*b[,1]*b[,2] +
      curr_T2[,2]*S2_X_adj - 
      2*curr_T[,2]*S2_mu_X_adj - 
      2*curr_T[,2]*W_ast*S2_alpha_X_adj - 
      2*curr_T[,2]*curr_T[,1]*b[,1]*S2_gamma_X_adj - 
      2*curr_T[,2]*curr_T[,2]*b[,2]*S2_omega_X_adj
    
  }
  
  out <- list(var_Yhat_new_train=var_Yhat_new_train)
  
  return(out)
  
}