b_vec_func <- function(W_ast, ID, Tt_invR, Y_split, fit, intercept_split, X_split, Tt_invR_T, sigma2_lmm, b_var, number_re){

  ##### We calculate the random effect b ##### 
  W_ast_split <- lapply(split(W_ast, ID), matrix, ncol=1)
  
  by_j_calc_part2 <- sfLapply(1:length(unique(ID)), fun=by_j_calc_part2_fun, Tt_invR=Tt_invR, 
                              Y_split=Y_split, fit=fit, intercept_split=intercept_split, 
                              W_ast_split=W_ast_split, X_split=X_split)
  
  b_vec <- sfLapply(1:length(unique(ID)), fun=b_vec_fun,  Tt_invR_T=Tt_invR_T, sigma2_lmm_value=sigma2_lmm, b_var_value=b_var, by_j_calc_part2=by_j_calc_part2)
  if(number_re == 1){
    b_vec_int <- unlist(b_vec)
    b_vec_slope <- NULL
  } else{
    b_vec_int <- unlist(b_vec)[seq(1,length(unlist(b_vec))-1, 2)]
    b_vec_slope <- unlist(b_vec)[seq(2,length(unlist(b_vec)), 2)]
  }
  
  out <- list(W_ast_split=W_ast_split, b_vec=b_vec, b_vec_int=b_vec_int, b_vec_slope=b_vec_slope)
  return(out)
}




