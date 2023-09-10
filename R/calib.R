calib <- function(W_ast, ID, W_ast_var, sigma2_lmm, Tt_invR_T, b_var, Tt_invR, 
                  curr_T_split, curr_T_split_t, b_vec, intercept, curr_T, b, 
                  number_re, X){
  
  ##### Calibration model with current W #####
  W_ast_split <- lapply(split(W_ast, ID), matrix, ncol=1) # current iteration W
  W_ast_var_split <- lapply(split(W_ast_var, ID), matrix, ncol=1) # current iteration W var
  
  Vt_bb_by_j_part <- sfLapply(1:length(unique(ID)), fun=Vt_bb_by_j_part_fun,
                              sigma2_lmm_value=sigma2_lmm, Tt_invR_T=Tt_invR_T, 
                              b_var_value=b_var)
  
  Vt_Wb_by_j_part <- sfLapply(1:length(unique(ID)), fun=Vt_Wb_by_j_part_fun, 
                              Tt_invR_T=Tt_invR_T, sigma2_lmm_value=sigma2_lmm,
                              b_var_value=b_var, Tt_invR=Tt_invR, 
                              W_ast_var_split=W_ast_var_split)
  
  E_b1_squared <- unlist(sfLapply(1:length(unique(ID)), fun=E_b1_squared_fun,
                                  curr_T_split=curr_T_split, Vt_bb_by_j_part=Vt_bb_by_j_part,
                                  curr_T_split_t=curr_T_split_t, b_vec=b_vec)) #Var b1 + Eb1*Eb1
  
  E_Wb1_c <- unlist(sfLapply(1:length(unique(ID)), fun=E_Wb1_c_fun, 
                             curr_T_split=curr_T_split, Vt_Wb_by_j_part=Vt_Wb_by_j_part, 
                             W_ast_split=W_ast_split, b_vec=b_vec)) #Cov W-b1 + EW*Eb1
  
  if(number_re == 1){
    E_b2_squared = E_b1b2_c = E_Wb2_c <- NULL
  } else{
    E_b2_squared <- unlist(sfLapply(1:length(unique(ID)), fun=E_b2_squared_fun,
                                    curr_T_split=curr_T_split, Vt_bb_by_j_part=Vt_bb_by_j_part,
                                    curr_T_split_t=curr_T_split_t, b_vec=b_vec)) #Var b2 + Eb2*Eb2
    
    E_b1b2_c <- unlist(sfLapply(1:length(unique(ID)), fun=E_b1b2_c_fun, 
                                curr_T_split=curr_T_split, Vt_bb_by_j_part=Vt_bb_by_j_part, 
                                curr_T_split_t=curr_T_split_t, b_vec=b_vec)) #Cov b1-b2 + Eb1*Eb2
    
    E_Wb2_c <- unlist(sfLapply(1:length(unique(ID)), fun=E_Wb2_c_fun, 
                               curr_T_split=curr_T_split, Vt_Wb_by_j_part=Vt_Wb_by_j_part, 
                               W_ast_split=W_ast_split, b_vec=b_vec)) #Cov W-b2 + EW*Eb2
  }
  
  if(number_re == 1){
    Xt_c <- as.matrix(cbind(intercept, W_ast, b[,1]))
  } else{
    Xt_c <- as.matrix(cbind(intercept, W_ast, curr_T[,1]*b[,1], curr_T[,2]*b[,2], X$x_adj))
  }
  
  XXt_c <- t(Xt_c)%*%Xt_c
  XXt_c[2,2] <- sum(W_ast_var+(W_ast^2))
  XXt_c[2,3] <- sum(E_Wb1_c)
  XXt_c[3,3] <- sum(E_b1_squared)
  XXt_c[3,2] <- sum(E_Wb1_c)
  
  if(number_re == 2){
    XXt_c[4,4] <- sum(E_b2_squared)
    XXt_c[2,4] <- sum(E_Wb2_c) 
    XXt_c[4,2] <- sum(E_Wb2_c) 
    XXt_c[3,4] <- sum(E_b1b2_c) 
    XXt_c[4,3] <- sum(E_b1b2_c) 
  }
  
  XXt_c_inv <- solve(XXt_c)
  c_coefs <- XXt_c_inv%*%(t(Xt_c)%*%Y)
  
  out <- list(c_coefs=c_coefs, Xt_c=Xt_c, XXt_c_inv=XXt_c_inv, 
              Vt_Wb_by_j_part=Vt_Wb_by_j_part)
  
  return(out)
  
}