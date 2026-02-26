lfdr_T_GK_simulation <- function(T, pi0 = NULL, T_old=NULL, trunc = TRUE, monotone = TRUE,
                      df_val = NULL, one.sided = FALSE, ...) {
  
  if(is.null(T_old)){T_old <- T}
  lfdr_out <- rep(1, length(T) )
  rm_na <- !is.na(T)
  T <- T[rm_na]
  if (is.null(pi0)) {
    stop("pi0 must be given.")
  }
  min_val <- min(c(T,-5) )
  max_val <- max(c(T, 5) )
  n_vals <- (max_val - min_val)*1000 + 1
  T_seq <- seq(min_val,max_val,length.out = n_vals)
  myd <- density(T_old, kernel = "gaussian", bw = 1, n = n_vals,
                 from = min_val, to = max_val)
  mys <- smooth.spline(x = myd$x, y = myd$y)
  y <- predict(mys, T_seq)$y
  y[y <= 0] <- 1e-10
  pi0_trunc <- min(c(max(c(pi0, predict(mys, 0)$y/dt(0, df = df_val) )), 1) )
  lfdr <- pi0_trunc * dt(T_seq, df = df_val)/y
  
  if (trunc) {
    lfdr[lfdr > 1] <- 1
  }
  
  if(FALSE){
    plot(T_seq,y,type="l", ylim=c(0,1))
    lines(T_seq,pi0_trunc * dt(T_seq, df = df_val), col=2)
    lines(T_seq,1-lfdr,col=3)
  }
  
  
  if (monotone & !one.sided) {
    T_neg <- T_seq[T_seq < 0]
    lfdr_neg <- lfdr[T_seq < 0]
    T_pos <- T_seq[T_seq >= 0]
    lfdr_pos <- lfdr[T_seq >= 0]
    
    o <- order(T_neg, decreasing = FALSE)
    ro <- order(o)
    lfdr[T_seq < 0] <- cummax(lfdr_neg[o])[ro]
    
    o <- order(T_pos, decreasing = TRUE)
    ro <- order(o)
    lfdr[T_seq >= 0] <- cummax(lfdr_pos[o])[ro]
  }
  if (monotone & one.sided) {
    o <- order(T_seq, decreasing = TRUE)
    ro <- order(o)
    lfdr <- cummax(lfdr[o])[ro]
  }
  my_lfdr <- approxfun(x = T_seq, y = lfdr)
  lfdr_vals <- my_lfdr(T)
  
  lfdr_out[rm_na] <- lfdr_vals
  res <- list(lfdr = lfdr_out, f = my_lfdr, pi0_trunc = pi0_trunc)
  return(res)
}