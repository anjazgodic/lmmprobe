pi0_func_simulation <- function (p, lambda = seq(0.05, 0.95, 0.05), pi0.method = c("smoother", 
                                                                      "bootstrap"), smooth.df = 3, smooth.log.pi0 = FALSE, ...) 
{
  rm_na <- !is.na(p)
  p <- p[rm_na]
  pi0.method = match.arg(pi0.method)
  m <- length(p)
  lambda <- sort(lambda)
  ll <- length(lambda)
  if (min(p) < 0 || max(p) > 1) {
    stop("ERROR: p-values not in valid range [0, 1].")
  }
  else if (ll > 1 && ll < 4) {
    stop(sprintf(paste("ERROR:", paste("length(lambda)=", 
                                       ll, ".", sep = ""), "If length of lambda greater than 1,", 
                       "you need at least 4 values.")))
  }
  else if (min(lambda) < 0 || max(lambda) >= 1) {
    stop("ERROR: Lambda must be within [0, 1).")
  }
  if (ll == 1) {
    pi0 <- mean(p >= lambda)/(1 - lambda)
    pi0.lambda <- pi0
    pi0 <- min(pi0, 1)
    pi0Smooth <- NULL
  }
  else {
    ind <- length(lambda):1
    pi0 <- cumsum(tabulate(findInterval(p, vec = lambda))[ind])/(length(p) * 
                                                                   (1 - lambda[ind]))
    pi0 <- pi0[ind]
    pi0.lambda <- pi0
    if (pi0.method == "smoother") {
      if (smooth.log.pi0) {
        pi0 <- log(pi0)
        spi0 <- smooth.spline(lambda, pi0, df = smooth.df)
        pi0Smooth <- exp(predict(spi0, x = lambda)$y)
        pi0 <- min(pi0Smooth[ll], 1)
      }
      else {
        spi0 <- smooth.spline(lambda, pi0, df = smooth.df)
        pi0Smooth <- predict(spi0, x = lambda)$y
        pi0 <- min(pi0Smooth[ll], 1)
      }
    }
    else if (pi0.method == "bootstrap") {
      minpi0 <- quantile(pi0, prob = 0.1)
      W <- sapply(lambda, function(l) sum(p >= l))
      mse <- (W/(m^2 * (1 - lambda)^2)) * (1 - W/m) + (pi0 - 
                                                         minpi0)^2
      pi0 <- min(pi0[mse == min(mse)], 1)
      pi0Smooth <- NULL
    }
    else {
      stop("ERROR: pi0.method must be one of \"smoother\" or \"bootstrap\".")
    }
  }
  if (pi0 < 0) {
    stop("ERROR: The estimated pi0 <= 0. Check that you have valid p-values or use a different range of lambda.")
  }
  return(list(pi0 = pi0, pi0.lambda = pi0.lambda, lambda = lambda, 
              pi0.smooth = pi0Smooth))
}