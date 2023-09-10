b2_vec_func <- function(ID, sigma2_lmm, Tt_invR_T, b_var, b_vec, b_vec_int, b_vec_slope){
  
  b2_vec <-  sfLapply(1:length(unique(ID)), fun=b2_vec_fun, sigma2_lmm, Tt_invR_T, b_var, b_vec)
  
  b_vec_all <- data.frame(cbind(b_vec_int, b_vec_slope))
  b_vec_out <- b_vec_all[rep(seq_len(nrow(b_vec_all)), each = unique(table(ID))), , drop=F] 
  
  b2_vec_collapse1 = lapply(b2_vec, function(x) { replicate(unique(table(ID)), x, simplify=FALSE) })
  b2_vec_collapse2 <- lapply(b2_vec_collapse1, function(x){ do.call(rbind, x) })
  b2_vec_out <- do.call(rbind, b2_vec_collapse2)
  rm(b2_vec_collapse1)
  rm(b2_vec_collapse2)
  
  out <- list(b2_vec=b2_vec, b_vec_all=b_vec_all, b_vec_out=b_vec_out, b2_vec_out=b2_vec_out)
  return(out)
  
}





