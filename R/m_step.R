m_step_it2plus <- function(M, N, W_ast_var, ID, number_re, curr_T, b, sigma2_lmm, Tt_invR_T, 
                b_var, Tt_invR, curr_T_split, curr_T_split_t, Y, Z, W_ast, delta, 
                beta_vec, X, mu_m, beta_m, beta_stderr, T_vals, alpha_m, gamma_m, 
                omega_m, beta_adj){
  
  ##### M-step #####
  t_Dt1 = t_Dt2 = t_Dt3 <- matrix(NA, ncol = M, nrow = N)
  Vt1_1 = Vt2_2 = Vt3_3 = Vt1_2 = Vt2_1 = Vt1_3 = Vt3_1 = Vt2_3 = Vt3_2 <-  matrix(NA, ncol = M, nrow = N)
  W_ast_var_split <- lapply(split(W_ast_var, ID), matrix, ncol=1) #this is V from my notes
  
  if(number_re == 1){
    t_Dt2 <- curr_T*b
    t_Dt3 <- NULL
  } else{
    t_Dt2 <- (curr_T[,1]*b[,1])
    t_Dt3 <- (curr_T[,2]*b[,2])  
  }
  
  Vt_bb_by_j_part <- sfLapply(1:length(unique(ID)), fun=Vt_bb_by_j_part_fun,
                              sigma2_lmm_value=sigma2_lmm, Tt_invR_T=Tt_invR_T, 
                              b_var_value=b_var)
  
  Vt_Wb_by_j_part <- sfLapply(1:length(unique(ID)), fun=Vt_Wb_by_j_part_fun, 
                              Tt_invR_T=Tt_invR_T, sigma2_lmm_value=sigma2_lmm,
                              b_var_value=b_var, Tt_invR=Tt_invR, 
                              W_ast_var_split=W_ast_var_split)
  
  Vt1_1 <- W_ast_var #Var W
  Vt2_2 <- unlist(sfLapply(1:length(unique(ID)), fun=Vt2_2_fun, curr_T_split, Vt_bb_by_j_part, curr_T_split_t)) #Var b1
  Vt1_2 = Vt2_1 <- unlist(sfLapply(1:length(unique(ID)), fun=Vt1_2_fun, curr_T_split, Vt_Wb_by_j_part)) #Cov Wb1, Cov b1W
  
  if(number_re == 1){
    Vt2_3 = Vt3_2 = Vt3_3 = Vt1_3 <- NULL
    Vt <- cbind(c(Vt1_1, Vt2_1), c(Vt1_2, Vt2_2), c(Vt1_3, Vt2_3))
  } else{
    Vt1_3 = Vt3_1 <- unlist(sfLapply(1:length(unique(ID)), fun=Vt1_3_fun, curr_T_split, Vt_Wb_by_j_part)) #Cov Wb2, Cov b2W
    Vt2_3 <- unlist(sfLapply(1:length(unique(ID)), fun=Vt2_3_fun, curr_T_split, Vt_bb_by_j_part, curr_T_split_t)) #Cov b1b2
    Vt3_2 <- unlist(sfLapply(1:length(unique(ID)), fun=Vt3_2_fun, curr_T_split, Vt_bb_by_j_part, curr_T_split_t)) #Cov b2b1
    Vt3_3 <- unlist(sfLapply(1:length(unique(ID)), fun=Vt3_3_fun, curr_T_split, Vt_bb_by_j_part, curr_T_split_t)) #Var b2
    Vt <- cbind(c(Vt1_1, Vt2_1, Vt3_1), c(Vt1_2, Vt2_2, Vt3_2), c(Vt1_3, Vt2_3, Vt3_3))
  }
  
  rowind1 <- (1:N)-1
  rowind2 <- ((N+1):(2*N))-1
  rowind3 <- (((2*N)+1):(3*N))-1
  colind1 <- 1-1
  colind2 <- 2-1
  colind3 <- 3-1
  
  if(number_re == 1){
    LR_cpp <- UNHIDEM_cpp0_5(Y, as.matrix(Z), as.matrix(Vt),
                             rowind1, rowind2, 
                             colind1, colind2,
                             W_ast, as.matrix(t_Dt2),  
                             as.matrix(Vt1_2),
                             W_ast_var, delta, 
                             beta_vec, sigma2_lmm)
  } else{
    LR_cpp <- UNHIDEM_cpp0_5_w_covs(Y, as.matrix(Z), as.matrix(Vt),
                                    as.matrix(X), 
                                    rowind1, rowind2, rowind3, 
                                    colind1, colind2, colind3,
                                    W_ast, t_Dt2, t_Dt3, 
                                    Vt1_2, Vt1_3, Vt2_3, 
                                    W_ast_var, delta, 
                                    beta_vec, sigma2_lmm)
  }
  
  mu_m <- LR_cpp$Coefficients[,1]
  beta_m <- LR_cpp$Coefficients[,2]
  beta_stderr <- LR_cpp$StdErr[,2]
  T_vals <- LR_cpp$T_statistics[,2]
  alpha_m <- LR_cpp$Coefficients[,3]
  gamma_m <- LR_cpp$Coefficients[,4]
  
  if(number_re == 1){
    omega_m <- NULL
    beta_adj <- NULL
  } else{
    omega_m <- LR_cpp$Coefficients[,5]
    omega_m[is.na(omega_m)|is.nan(omega_m)] <- 0
    beta_adj <- LR_cpp$Coefficients[,6]
    beta_adj[is.na(beta_adj)|is.nan(beta_adj)] <- 0
  }
  
  mu_m[is.na(mu_m)|is.nan(mu_m)] <- 0
  beta_m[is.na(beta_m)|is.nan(beta_m)] <- 0
  alpha_m[is.na(alpha_m)|is.nan(alpha_m)] <- 0
  gamma_m[is.na(gamma_m)|is.nan(gamma_m)] <- 0
  beta_stderr[is.na(beta_stderr)|is.nan(beta_stderr)] <- 1
  T_vals[is.na(T_vals)|is.nan(T_vals)] <- 0
  
  out = list(mu_m=mu_m, beta_m=beta_m, alpha_m=alpha_m, gamma_m=gamma_m, 
             omega_m=omega_m, beta_adj=beta_adj, beta_stderr=beta_stderr, 
             T_vals=T_vals)
  return(out)
  
}
