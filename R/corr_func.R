corr_func_simulation <- function(beta_t_new, beta_t, T_vals, B = 5, one.sided = FALSE) {
  
  ##### We clean up the T values first ##### 
  T_vals[is.na(T_vals)] <- 0
  T_vals[is.nan(T_vals)] <- 0
  if (!one.sided) {
    T_vals <- abs(T_vals)
  }
  xx <- data.frame(a = beta_t_new, b = beta_t, T_vals = T_vals)
  ##### We figure out how to break up the beta estimates into B groups ##### 
  brks <- c(min(T_vals), max(T_vals))
  if (B > 1) {
    min_prob <- 1/2 * I(one.sided) + 1/B * I(!one.sided)
    jumps <- 1/(B * (1 + 1 * I(one.sided)))
    brks <- with(xx, quantile(T_vals, probs = c(0, seq(min_prob, 1, 
                                                       jumps)), na.rm = TRUE))
  }
  if(length(brks) != length(unique(brks))) { 
    brks = brks + seq_along(brks) * .Machine$double.eps 
  } 
  ##### Labels the beta coefficients according to their group #####
  xx <- within(xx, group <- cut(T_vals, breaks = brks, labels = 1:(length(brks) - 1), 
                                include.lowest = TRUE))
  ##### Take the correlation between current and previous betas, within each of the B groups #####
  result <- by(xx[, 1:2], xx$group, function(x) {
    cor(x$a, x$b)
  })
  result.dataframe <- as.data.frame(as.matrix(result))
  cor_vec <- result.dataframe$V1[as.numeric(xx$group)]
  cor_vec[is.na(cor_vec)] <- 0
  return(cor_vec)
}