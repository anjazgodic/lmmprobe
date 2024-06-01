b2_vec_func <- function(ID, sigma2_lmm, Tt_invR_T, b_var, b_vec, b_vec_int, b_vec_slope){

  b2_vec <-  sfLapply(1:length(unique(ID)), fun=b2_vec_fun, sigma2_lmm, Tt_invR_T, b_var, b_vec)
  
  b_vec_all <- data.frame(cbind(b_vec_int, b_vec_slope))
  b_vec_out <- b_vec_all[rep(seq_len(nrow(b_vec_all)), table(ID)), , drop=F] 
  
  b2_vec_collapse1 = sapply(1:length(b2_vec), function(x) { replicate(table(ID)[x], b2_vec[[x]], simplify=FALSE) })
  if(any(class(b2_vec_collapse1) %in% c("list"))){
    b2_vec_collapse2 <- lapply(b2_vec_collapse1, function(x){ do.call(rbind, x) })
  } else{
    b2_vec_collapse2 = b2_vec_collapse1
  }
  b2_vec_out <- do.call(rbind, b2_vec_collapse2)
  rm(b2_vec_collapse1)
  rm(b2_vec_collapse2)
  
  out <- list(b2_vec=b2_vec, b_vec_all=b_vec_all, b_vec_out=b_vec_out, b2_vec_out=b2_vec_out)
  return(out)
  
}





