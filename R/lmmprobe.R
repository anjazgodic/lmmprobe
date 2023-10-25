#' Sparse high-dimensional linear mixed modeling with PaRtitiOned empirical Bayes ECM (LMM-PROBE) algorithm. 
#' 
#' @name lmmprobe
#' @description Sparse high-dimensional linear mixed modeling with PaRtitiOned empirical Bayes ECM (LMM-PROBE) algorithm. Currently, the package offers functionality for two scenarios. Scenario 1: only a random intercept, each unit has the same number of observations; Scenario 2: a random intercept and a random slope, each unit has the same number of observations. We are actively expanding the package for more flexibility and scenarios. 
#' @param Y A matrix containing the outcome \code{Y}.
#' @param Z A matrix containing the sparse fixed-effect predictors on which to apply the lmmprobe algorithm.
#' @param V A matrix containing non-sparse predictors for the random effects. This matrix is currently only programmed for two scenarios. Scenario 1: only a random intercept, where V is a matrix with one column containing ID's and each unit has the same number of observations. Scenario 2: a random intercept and a random slope, where V is a matrix with two columns. The first column is ID and the second column is a continuous variable (e.g. time) for which a random slope is to be estimated. Each unit has the same number of observations. 
#' @param alpha Type I error; significance level.
#' @param ep Value against which to compare convergence criterion, we recommend 0.05.
#' @param B The number of groups to categorize estimated coefficients in to calculate correlation \eqn{\rho}. We recommend five.
#' @param adj A factor multiplying Silverman’s 'rule of thumb' in determining the bandwidth for density estimation, same as the 'adjust' argument of R's density function. Default is three. 
#' @param maxit Maximum number of iterations the algorithm will run for. Default is 10000.
#' @param cpus The number of CPUS user would like to use for parallel computations. Default is four. 
#' @param LR A learning rate parameter \code{r}. Using zero corresponds to the implementation described in Zgodic et al.  
#' @param C A learning rate parameter \code{c}. Using one corresponds to the implementation described in Zgodic et al. 
#' @return A list of the output of the lmmprobe function, including 
#' 
#' \code{beta_hat, beta_hat_var} MAP estimates of the posterior expectation (beta_hat) and variance (beta_hat_var) of the prior mean (\eqn{\beta}) of the regression coefficients assuming \eqn{\gamma=1}, 
#' 
#' \code{gamma} the posterior expectation of the latent \eqn{\gamma} variables, 
#' 
#' \code{preds} predictions of \eqn{Y},
#' 
#' \code{PI_lower, PI_upper} lower and upper prediction intervals for the predictions,
#' 
#' \code{sigma2_est} MAP estimate of the residual variance, 
#' 
#' \code{random_var} MAP estimate of the random effect(s) variance, 
#' 
#' \code{random_intercept} estimated random intercept terms,  
#' 
#' \code{random_slope} estimated random slope terms, if applicable. 
#' @examples
#' library(lmmprobe)
#' data(SLE)
#' ep <- 0.05
#' alpha <- 0.05
#' Y = SLE$Y
#' Z = SLE$Z
#' V = SLE$V
#' full_res <- lmmprobe(Y = Y, Z = Z, V = V, ep = ep, alpha = alpha)
#' @references Zgodic, A. and McLain, A. C. (2023). Sparse high-dimensional linear mixed modeling with a partitioned empirical Bayes ECM algorithm. arXiv preprint arXiv:2310.12285.
#' @export

utils::globalVariables(c(
  "by_j_calc_part2_fun", "b_vec_fun", 
  "b2_vec_fun", "e_fun", "trace_eRe_fun",
  "Vt2_2_fun", "Vt2_3_fun", "Vt3_2_fun", 
  "Vt3_3_fun", "Vt1_2_fun", "Vt1_3_fun",
  "Vt_bb_by_j_part_fun", "Vt_Wb_by_j_part_fun", 
  "E_b1_squared_fun", "E_b2_squared_fun", 
  "E_b1b2_c_fun", "E_Wb1_c_fun", "E_Wb2_c_fun",
  "e_calib_fun", "cov_Wb1_c_fun", 
  "cov_Wb2_c_fun", "cov_b1b2_c_fun",
  "Z", "Z_2", "beta_t", "delta", "W_ast", 
  "curr_T", "W_ast_var", "Tt_invR", 
  "Tt_invR_T", "ID", "curr_T_split",
  "sigma2_lmm", "b_var", "curr_T_split_t"
), add = F)

lmmprobe <- function(Y, Z, V, alpha = 0.05, maxit = 10000, 
                         ep = 0.05, cpus = 4,
                         B = 5, adj = 3, LR = 0, C = 1) {   
      
      ##### Setting initial values and initializing outputs #####
      one.sided = FALSE
      M <- dim(Z)[2]
      N <- dim(Z)[1]
      number_re = M_random <- dim(V)[2]
  
      intercept <- rep(1, N)
      
      if(number_re == 1) {
        colnames(V) <- c("id")
        X <- NULL
        RE_dat <- matrix(intercept, ncol = 1)
        colnames(RE_dat) <- c("intercept")
      } else{
        colnames(V) <- c("id", "x_adj")
        X <- V[,"x_adj", drop=F]
        RE_dat <- matrix(cbind(intercept, X[,"x_adj"]), ncol = 2)
        colnames(RE_dat) <- c("intercept", "x_adj")
      }
      
      ID <- V[,"id"]
      
      p <- 0
      if (!is.null(X)) {
        p <- dim(X)[2]
        beta_X <- rep(0, p)
      }
      
      delta <- beta_t <- beta_adj <- alpha_vec <- beta_var <- 
        gamma_vec <- omega_vec <- beta_tilde <- beta_tilde_var <- 
        mu_m <- beta_m <- alpha_m <- gamma_m <- omega_m <- 
        StdErr <- T_vals <- beta_stderr <- rep(0,M)
      b_var = b_var_value <- diag(1, number_re)
      W_ast_store = W_ast_var_store <- matrix(0, ncol = M, nrow = N) 
      b_vec_out = b2_vec_out <- matrix() 
      b_vec = b2_vec <- list()
      R <- diag(rep(1, unique(table(ID))))
      invR <- replicate(length(unique(ID)), solve(R), simplify=FALSE)
      R <- replicate(length(unique(ID)), R, simplify=FALSE)
      W_ast = e = trace_eRe <- rep(0, N)
      W_ast_var <- W_ast + 1
      count <- 0
      conv_check <- 0
      cor_vec <- rep(0, M)
      Z_2 <- Z * Z
      X_2 <- X * X
      RE_dat2 <- RE_dat * RE_dat
      Xt_conv1 = Xt_conv1_W <- 1
      Xt_conv1_store = Xt_conv1_store_W = loglik_current <- c()
      try2 <- 0
      sigma2_lmm <- var(Y)
      
      ##########################################
      ############## Iteration 1 ###############
      ##########################################
      count <- 1
      # start_time <- Sys.time()
      
      ##### Set up some objects ##### 
      df = N - 3 - p
      curr_T <- RE_dat
      curr_T2 <- RE_dat2
      curr_T_split <- lapply(split(curr_T, ID), as.matrix)
      curr_T2_split <- lapply(split(curr_T2, ID), as.matrix)
      intercept_split <- lapply(split(intercept, ID), matrix, ncol=1)
      curr_T_split_t <- lapply(curr_T_split, t)
      curr_T2_split_t <- lapply(curr_T2_split, t)
      Tt_invR <- (Map("%*%", curr_T_split_t, invR))
      Tt_invR_T <- Map("%*%", Map("%*%", curr_T_split_t, invR), curr_T_split)
      Y_split <- lapply(split(Y, ID), matrix, ncol=1)
      if(number_re == 2){
        X_split <- lapply(split(X$x_adj, ID), matrix, ncol=1)
      } else{
        X_split <- NULL
      }
      
      ##### These are some functions we use snowfall for parallel computations ##### 
      sfInit(cpus = cpus, type = "SOCK", parallel = TRUE, 
             socketHosts=NULL, restore=NULL,
             slaveOutfile=NULL, nostart=FALSE, useRscript=FALSE)
      # sfExport("by_j_calc_part2_fun", "b_vec_fun", 
      #          "b2_vec_fun", "e_fun", "trace_eRe_fun",
      #          "Vt2_2_fun", "Vt2_3_fun", "Vt3_2_fun", 
      #          "Vt3_3_fun", "Vt1_2_fun", "Vt1_3_fun",
      #          "Vt_bb_by_j_part_fun", "Vt_Wb_by_j_part_fun", 
      #          "E_b1_squared_fun", "E_b2_squared_fun", 
      #          "E_b1b2_c_fun", "E_Wb1_c_fun", "E_Wb2_c_fun",
      #          "e_calib_fun", "cov_Wb1_c_fun", 
      #          "cov_Wb2_c_fun", "cov_b1b2_c_fun",
      #          "Z", "Z_2", "beta_t", "delta", "W_ast", 
      #          "curr_T", "W_ast_var", "Tt_invR", 
      #          "Tt_invR_T", "ID", "curr_T_split",
      #          "sigma2_lmm", "b_var", "curr_T_split_t", 
      #          local = F)
      sfExportAll()
      sfLibrary(tidyr)
      sfLibrary(lme4)
      
      ##########################################
      ########## M-Step at Iteration 1########## 
      ##########################################
      
      ##### M-Step #####
      if(number_re == 1){
        LR_cpp <- LM_by_col(Y, as.matrix(Z), sigma2_lmm)
        b_var <- b_var_value  
        sigma2_lmm <- min(LR_cpp$sig)
      } else{
        LR_cpp <- LM_by_col_w_covs(Y, as.matrix(Z), as.matrix(X), sigma2_lmm)
        b_var <- b_var_value
        sigma2_lmm <- min(LR_cpp$sig)
      }
      
      mu_m <- LR_cpp$Coefficients[,1]
      beta_m <- LR_cpp$Coefficients[,2]
      StdErr <- LR_cpp$StdErr[,2]
      T_vals <- beta_m/StdErr
      alpha_vec <- rep(0, M)
      gamma_vec <- rep(0, M)
      
      if(number_re == 1){
        beta_adj <- NULL
        omega_vec <- NULL
      } else{
        beta_adj <- LR_cpp$Coefficients[,3]
        omega_vec <- rep(0, M)
      }
      
      ##### Saving/renaming some items from the M-step #####
      beta_t_new <- beta_m
      beta_var_new <- StdErr^2
      beta_adj_vec <- beta_adj
      
      ###########################################################
      ########## There is no mixing step at Iteration 1########## 
      ###########################################################
      ##### Mixing Step #####
      beta_t <- beta_t_new/2
      beta_var <- beta_var_new/2
      
      ##########################################
      ########## E-Step at Iteration 1########## 
      ##########################################
      
      ##### E-Step #####
      ##### Get and clean up T statistics ##### 
      T_vals <- beta_t/sqrt(beta_var)
      T_vals[beta_var<=0] <- 0
      T_vals[is.na(T_vals)] <- 0
      T_vals[is.nan(T_vals)] <- 0   
      p_vals <- pt(abs(T_vals),df=df,lower.tail = FALSE)*2
      
      ##### Estimate the probability of a null ##### 
      p_hat  <- (pi0_func_simulation(p_vals,lambda = 0.1)$pi0*M + 1)/(M+2) 
      ##### Estimate the lfdr ##### 
      delta <- 1-lfdr_T_GK_simulation(T = T_vals, pi0 = p_hat,T_old = NULL, trunc = TRUE, monotone = TRUE,adj = adj, df_val = df, one.sided = one.sided)$lfdr
      
      ##### We update W so that we can do the E-Step ##### 
      W_ast_store <- MVM(as.matrix(Z),delta*beta_t)$Res 
      W_ast_store[is.nan(W_ast_store) | is.na(W_ast_store)] <- 0
      W_ast <- c(Row_sum(as.matrix(W_ast_store))$Rowsum)
      
      ##### We update Var W so that we can do the E-Step ##### 
      W_ast_var <- NULL
      if(!is.null(Z_2)){
        W_ast_var_store <- MVM(as.matrix(Z_2),(beta_tilde^2)*(delta*(1-delta)))$Res 
        W_ast_var_store[is.nan(W_ast_var_store) | is.na(W_ast_var_store)] <- 0
        W_ast_var <- c(Row_sum(as.matrix(W_ast_var_store))$Rowsum)
      }
      
      ##### We calculate b with lmer to start ##### 
      if(number_re == 1){
        fit <- lmer( Y ~ W_ast + (1 | ID)) 
      } else{
        fit <- lmer( Y ~ W_ast + X$x_adj + (1 + RE_dat[,2] | ID))  
      }
      
      ##### We calculate the random effect b ##### 
      b_vec_core <- b_vec_func(W_ast, ID, Tt_invR, Y_split, fit, intercept_split, X_split, Tt_invR_T, sigma2_lmm, b_var, number_re)
      W_ast_split <- b_vec_core$W_ast_split
      b_vec <- b_vec_core$b_vec
      b_vec_int <- b_vec_core$b_vec_int
      b_vec_slope <- b_vec_core$b_vec_slope
      
      ##### We calculate b2 (for e, sigma2_lmm, b_var later) ##### 
      b2_vec_core <- b2_vec_func(ID, sigma2_lmm, Tt_invR_T, b_var, b_vec, b_vec_int, b_vec_slope)
      b2_vec <- b2_vec_core$b2_vec
      b_vec_all <- b2_vec_core$b_vec_all
      b = b_vec_out <- b2_vec_core$b_vec_out
      b2 = b2_vec_out <- b2_vec_core$b2_vec_out
      
      ##### We calculate error e and trace eRe #####
      e <- gather(data.frame(sfSapply(1:length(unique(ID)), fun=e_fun, Y_split, fit, intercept_split, W_ast_split, curr_T_split, b_vec, X_split)))[,2]
      e_split <- lapply(split(e, ID), matrix, ncol=1)
      e_split_t <- lapply(e_split, t)
      et_invR_e <- unlist(Map("%*%", Map("%*%", e_split_t, invR), e_split))
      trace_eRe <- sfSapply(1:length(unique(ID)), fun=trace_eRe_fun, Tt_invR_T, sigma2_lmm, b_var, et_invR_e)
      
      ##### Saving outputs from the E-step #####
      beta_tilde <- beta_t
      beta_tilde_var <- beta_var
      
      ##### Finishing up E-step #####
      if (sum(delta) == 0) {
        print("Delta is equal to zero")
        if(try2 == 0){
          cat("Warning loop completely recycled back to beta=0.\n 
                      Trying again with different starting values. \n")
          count <- 0
          beta_tilde <- rep(0.0001,M)
          beta_tilde_var <- rep(0,M)  
          delta <- rep(1,M)
          beta_adj_vec <- rep(0.0001,M)
          alpha_vec <- rep(0.0001,M)
          gamma_vec <- rep(0.0001,M)
          omega_vec <- rep(0.0001,M)
          cor_vec <- rep(0, M)  
          b_vec <- lapply(1:length(unique(ID)), matrix, data=rep(0, number_re), nrow=number_re, ncol=1)
          b_var <- b_var_value
          sigma2_lmm <- sigma2_lmm_value
          Xt_conv1 <- 1 
          try2 <- 1
          W_ast <- apply(Z,1,sum)*0.0001
          W_ast_var <- W_ast*0
        } else {
          try2 <- 2
          conv_check <- conv_check + 1
          W_ast <- apply(Z,1,sum)*0.0001
          W_ast_var <- W_ast*0
        }
      } else{
        ##### Updating W, W2, b, b2 #####
        W_ast_var <- NULL
        if(!is.null(Z_2)){
          W_ast_var_store <- MVM(as.matrix(Z_2),(beta_tilde^2)*(delta*(1-delta)))$Res 
          W_ast_var_store[is.nan(W_ast_var_store) | is.na(W_ast_var_store)] <- 0
          W_ast_var <- c(Row_sum(as.matrix(W_ast_var_store))$Rowsum)
        }
      } #end of big if (sum(delta) == 0)  statement
      
      ##### Checking convergence and iterations #####
      conv <- 0
      if ((conv_check == 0 & try2 < 2)|conv_check == 0) {
        conv <- 1
      }
      if (try2 == 2) {
        conv <- 2
        cat("Warning: loop completely recycled back to beta=delta=0 twice. Optimization failed.\n")
      }
      
      ##########################################
      ############ Iteration 2,... #############
      ##########################################
    
      ##### Starting the while loop #####
      while (count < maxit & conv_check < 1) {
        
        ##### Setting up some objects #####
        beta_vec = beta_tilde
        beta_var = beta_tilde_var
        T_old <- T_vals
        loglik_old <- loglik_current
        W_ast_old <- W_ast
        W_ast2_old <- W_ast_var
        count <- count + 1
        fact <- C*(count + 1)^(-(1+LR)) 
        
        ##########################################
        ######### M-Step at Iteration >1 #########
        ##########################################
        
        ##### M-step #####
        m_step_out <- m_step_it2plus(M, N, W_ast_var, ID, number_re, curr_T, b, sigma2_lmm, Tt_invR_T, 
                                     b_var, Tt_invR, curr_T_split, curr_T_split_t, Y, Z, W_ast, delta, 
                                     beta_vec, X, mu_m, beta_m, beta_stderr, T_vals, alpha_m, gamma_m, 
                                     omega_m, beta_adj)
        mu_m <- m_step_out$mu_m
        beta_m <- m_step_out$beta_m
        alpha_m <- m_step_out$alpha_m
        gamma_m <- m_step_out$gamma_m
        omega_m <- m_step_out$omega_m
        beta_adj <- m_step_out$beta_adj
        beta_stderr <- m_step_out$beta_stderr
        T_vals <- m_step_out$T_vals
        
        ##### This part estimates the (co)variance matrix b_var of the random effect (1x1 if only one random effect) #####
        if(number_re == 1){
          b_var <- mean(unlist(b2_vec))
        } else{
          b_var[1,1] <- mean(unlist(lapply(b2_vec, function(x) x[1,1])))
          b_var[2,2] <- mean(unlist(lapply(b2_vec, function(x) x[2,2])))
          b_var[1,2] <- mean(unlist(lapply(b2_vec, function(x) x[1,2])))
          b_var[2,1] <- mean(unlist(lapply(b2_vec, function(x) x[2,1])))
        }
  
        ##### This part estimates the error of the EM-LMM model because it is needed in calculations of b and b2 #####
        sigma2_lmm <- sum(trace_eRe)*(1/N) 
        
        ##### Saving some items from the M-step #####
        beta_t_new <- beta_m
        beta_var_new <- beta_stderr^2
        mu_vec <- mu_m
        alpha_vec <- alpha_m
        gamma_vec <- gamma_m
        omega_vec <- omega_m
        beta_adj_vec <- beta_adj
        
        ##########################################
        ##### Mixing Step at Iteration 2,...######
        ########################################## 
        
        ##### Mixing Step #####
        ##### Getting correlation rho #####
        if (count > 1) {
          cor_vec <- corr_func_simulation(beta_t_new, beta_vec, T_vals, B = B, one.sided = one.sided)
        }
        ##### Editing beta and variance(beta) using rho and the mixing step #####
        beta_t <- beta_t_new * fact + beta_vec * (1 - fact)
        beta_var <- ((beta_var_new^(-1) * fact) + (beta_var^(-1) * (1 - fact)))^(-1) 
        
        ##########################################
        ######## E-Step at Iteration 2,...########
        ##########################################  
        
        ##### E-step #####
        ##### Get and clean up T statistics ##### 
        T_vals <- beta_t/sqrt(beta_var)
        T_vals[beta_var<=0] <- 0
        T_vals[is.na(T_vals)] <- 0
        T_vals[is.nan(T_vals)] <- 0   
        M <- length(T_vals) 
        p_vals <- pt(abs(T_vals),df=df,lower.tail = FALSE)*2
        
        ##### Estimate the probability of a null ##### 
        p_hat  <- (pi0_func_simulation(p_vals,lambda = 0.1)$pi0*M + 1 )/(M+2)
        ##### Estimate the lfdr ##### 
        delta <- 1-lfdr_T_GK_simulation(T = T_vals, pi0 = p_hat,T_old = T_old, trunc = TRUE, monotone = TRUE,adj = adj, df_val = df, one.sided = one.sided)$lfdr
    
        ##### We update W so that we can do the E-Step ##### 
        W_ast_store <- MVM(as.matrix(Z),delta*beta_t)$Res 
        W_ast_store[is.nan(W_ast_store) | is.na(W_ast_store)] <- 0
        W_ast <- c(Row_sum(as.matrix(W_ast_store))$Rowsum)
        
        ##### We update WVar (for W2) so that we can do the E-Step #####
        W_ast_var <- NULL
        if(!is.null(Z_2)){
          W_ast_var_store <- MVM(as.matrix(Z_2),(beta_t^2)*(delta*(1-delta)))$Res 
          W_ast_var_store[is.nan(W_ast_var_store) | is.na(W_ast_var_store)] <- 0
          W_ast_var <- c(Row_sum(as.matrix(W_ast_var_store))$Rowsum)
        }
  
        if(all(W_ast==0)){
          print("all W are equal to zero")
          W_ast <- apply(Z,1,sum)*0.0001
          W_ast_var <- W_ast*0
        }
        
        ##### Calibration model with current W #####
        c_coefs <- calib(W_ast, ID, W_ast_var, sigma2_lmm, Tt_invR_T, b_var, Tt_invR, 
                         curr_T_split, curr_T_split_t, b_vec, intercept, curr_T, b, 
                         number_re, X)
        Vt_Wb_by_j_part <- c_coefs$Vt_Wb_by_j_part
        c_coefs <- c_coefs$c_coefs
        
        ##### We calculate error e and trace eRe #####
        e <- gather(data.frame(sfSapply(1:length(unique(ID)), fun=e_calib_fun,
                                        Y_split=Y_split, c_coefs=c_coefs,
                                        intercept_split=intercept_split,
                                        W_ast_split=W_ast_split, curr_T_split=curr_T_split,
                                        b_vec=b_vec, X_split=X_split)))[,2]
        e_split <- lapply(split(e, ID), matrix, ncol=1)
        e_split_t <- lapply(e_split, t)
        et_invR_e <- unlist(Map("%*%", Map("%*%", e_split_t, invR), e_split))
        
        cov_Wb1_c <- unlist(sfLapply(1:length(unique(ID)), fun=cov_Wb1_c_fun, curr_T_split, Vt_Wb_by_j_part)) #Cov W-b1
        if(number_re == 1){
          cov_Wb2_c <- NULL
        } else{
          cov_Wb2_c <- unlist(sfLapply(1:length(unique(ID)), fun=cov_Wb2_c_fun, curr_T_split, Vt_Wb_by_j_part)) #Cov W-b2
        }
        trace_eRe <- sfSapply(1:length(unique(ID)), fun=trace_eRe_fun, Tt_invR_T, sigma2_lmm, b_var, et_invR_e)
        sigma2_lmm <- sum(trace_eRe)*(1/N)
        
        ##### We calculate the random effect b ##### 
        b_vec_core <- b_vec_func(W_ast, ID, Tt_invR, Y_split, fit=c_coefs, intercept_split, X_split, Tt_invR_T, sigma2_lmm, b_var, number_re)
        W_ast_split <- b_vec_core$W_ast_split
        b_vec <- b_vec_core$b_vec
        b_vec_int <- b_vec_core$b_vec_int
        b_vec_slope <- b_vec_core$b_vec_slope
        
        ##### We calculate b2 (for e, sigma2_lmm, b_var later) ##### 
        b2_vec_core <- b2_vec_func(ID, sigma2_lmm, Tt_invR_T, b_var, b_vec, b_vec_int, b_vec_slope)
        b2_vec <- b2_vec_core$b2_vec
        b_vec_all <- b2_vec_core$b_vec_all
        b = b_vec_out <- b2_vec_core$b_vec_out
        b2 = b2_vec_out <- b2_vec_core$b2_vec_out
        
        ##### This part estimates the (co)variance matrix b_var of the random effect (1x1 if only one random effect) #####
        if(number_re == 1){
          b_var <- mean(unlist(b2_vec))
        } else{
          b_var[1,1] <- mean(unlist(lapply(b2_vec, function(x) x[1,1])))
          b_var[2,2] <- mean(unlist(lapply(b2_vec, function(x) x[2,2])))
          b_var[1,2] <- mean(unlist(lapply(b2_vec, function(x) x[1,2])))
          b_var[2,1] <- mean(unlist(lapply(b2_vec, function(x) x[2,1])))
        }
  
        ##### Calibration model with current W and current b ##### 
        c_coefs <- calib(W_ast, ID, W_ast_var, sigma2_lmm, Tt_invR_T, b_var, Tt_invR, 
                         curr_T_split, curr_T_split_t, b_vec, intercept, curr_T, b, 
                         number_re, X)
        Vt_Wb_by_j_part <- c_coefs$Vt_Wb_by_j_part
        c_coefs <- c_coefs$c_coefs
        
        ##### We calculate error e and trace eRe #####
        e <- gather(data.frame(sfSapply(1:length(unique(ID)), fun=e_calib_fun,
                                        Y_split=Y_split, c_coefs=c_coefs,
                                        intercept_split=intercept_split,
                                        W_ast_split=W_ast_split, curr_T_split=curr_T_split,
                                        b_vec=b_vec, X_split=X_split)))[,2]
        e_split <- lapply(split(e, ID), matrix, ncol=1)
        e_split_t <- lapply(e_split, t)
        et_invR_e <- unlist(Map("%*%", Map("%*%", e_split_t, invR), e_split))
        
        cov_Wb1_c <- unlist(sfLapply(1:length(unique(ID)), fun=cov_Wb1_c_fun, curr_T_split, Vt_Wb_by_j_part)) #Cov W-b1
        if(number_re == 1){
          cov_Wb2_c <- NULL
        } else{
          cov_Wb2_c <- unlist(sfLapply(1:length(unique(ID)), fun=cov_Wb2_c_fun, curr_T_split, Vt_Wb_by_j_part)) #Cov W-b2
        }
        trace_eRe <- sfSapply(1:length(unique(ID)), fun=trace_eRe_fun, Tt_invR_T, sigma2_lmm, b_var, et_invR_e)
        sigma2_lmm <- sum(trace_eRe)*(1/N)
        
        ##### Saving outputs from the E-step #####
        beta_tilde <- beta_t
        beta_tilde_var <- beta_var
        
        ##### Wrapping up E-Step #####
        if (sum(delta) == 0) {
          print("Delta is equal to zero")
          if(try2 == 0){
            cat("Warning loop completely recycled back to beta=0.\n 
                      Trying again with different starting values. \n")
            count <- 0
            beta_tilde <- rep(0.0001,M)
            beta_tilde_var <- rep(0,M)  
            delta <- rep(1,M)
            alpha_vec <- rep(0.0001,M)
            gamma_vec <- rep(0.0001,M)
            omega_vec <- rep(0.0001,M)
            beta_adj_vec <- rep(0.0001,M)
            cor_vec <- rep(0, M)  
            b_vec <- lapply(1:length(unique(ID)), matrix, data=rep(0, number_re), nrow=number_re, ncol=1)
            b_var <- b_var_value  #change this 
            sigma2_lmm <- sigma2_lmm_value #change this 
            Xt_conv1 <- 1 
            try2 <- 1
            W_ast <- apply(Z,1,sum)*0.0001
            W_ast_var <- W_ast*0
          } else {
            try2 <- 2
            conv_check <- conv_check + 1
            W_ast <- apply(Z,1,sum)*0.0001
            W_ast_var <- W_ast*0
          }
        } else{
          ##### Calculating convergence criterion #####
          if (count > 2) {
            if(any(W_ast2_old>0)){
              Xt_conv1_W <- pchisq(max((W_ast_old[W_ast2_old>0] - W_ast[W_ast2_old>0])^2 /
                                     W_ast2_old[W_ast2_old>0])/log(N), df = 1)
            } else{
              Xt_conv1_W <- tail(Xt_conv1_store_W, n=1)
            }
            
            Xt_conv1_store_W <- c(Xt_conv1_store_W, Xt_conv1_W)
            
            if ( Xt_conv1_W < ep ) {
              conv_check <- conv_check + 1
            }
          }
          if (try2 == 2) {
            conv <- 2
            cat("Warning: loop completely recycled back to beta=delta=0 twice. Optimization failed.\n")
          }
          
        } #end of big if (sum(delta) == 0)  statement
  
        if((count %% 50) == 0) { cat("Iteration=", count, "Convergence Crit=", round(Xt_conv1_W, digits = 3), "\n") }
        
      } #end of while
      
      ##########################################
      ### Calibration Step at Last Iteration ###
      ##########################################  
      
      ##### Final calibration step on training data ##### 
      cat("Converged", "\n")
      cat("Iteration=", count, "Convergence Crit=", round(Xt_conv1_W, digits = 3), "\n")
      df <- N - (sum(delta) + 3 + p)
      if(df <= 0){df <- NULL}
      
      ##### This part estimates the last b_var #####
      if(number_re == 1){
        b_var <- mean(unlist(b2_vec))
      } else{
        b_var[1,1] <- mean(unlist(lapply(b2_vec, function(x) x[1,1])))
        b_var[2,2] <- mean(unlist(lapply(b2_vec, function(x) x[2,2])))
        b_var[1,2] <- mean(unlist(lapply(b2_vec, function(x) x[1,2])))
        b_var[2,1] <- mean(unlist(lapply(b2_vec, function(x) x[2,1])))
      }
  
      ##### This part estimates last error of the EM-LMM model #####
      sigma2_lmm <- sum(trace_eRe)*(1/N)
      
      ##### Calibration model with current W and current b ##### 
      c_coefs <- calib_final(W_ast, ID, W_ast_var, sigma2_lmm, Tt_invR_T, b_var, Tt_invR, 
                       curr_T_split, curr_T_split_t, b_vec, intercept, curr_T, b, 
                       number_re, X)
      Xt_c <- c_coefs$Xt_c
      XXt_c_inv <- c_coefs$XXt_c_inv
      Vt_bb_by_j_part <- c_coefs$Vt_bb_by_j_part
      cov_Wb1_c <- c_coefs$cov_Wb1_c
      cov_b1b2_c <- c_coefs$cov_b1b2_c
      cov_Wb2_c <- c_coefs$cov_Wb2_c
      c_coefs <- c_coefs$c_coefs
      
      Y_pred <- Xt_c %*% c_coefs
      sigma2_c <- sum((Y - Y_pred)^2)/df
      VCV <- sigma2_c * t(XXt_c_inv) %*% (t(Xt_c) %*% Xt_c) %*% XXt_c_inv #sandwich style estimator
      
      ##### Variance of predictions for data ##### 
      prediction_var <- pred_var(W_ast, ID, W_ast_var, VCV, curr_T_split, curr_T_split_t, 
                           curr_T, curr_T2, b, number_re, c_coefs, 
                           sigma2_lmm, Tt_invR_T, b_var, Vt_bb_by_j_part, 
                           cov_Wb1_c, cov_b1b2_c, cov_Wb2_c)
      
      var_Yhat_new_train_sigma2_c <- (prediction_var$var_Yhat_new_train + sigma2_c)[,1]
      
      sfStop(nostop=FALSE)
            
      ##### Saving beta coefficients for each variable #####
      beta_hat <- (c_coefs[2])*beta_tilde
      beta_hat_var <- (c_coefs[2]^2)*beta_tilde_var
      
      results <- list(beta = beta_hat, 
                      beta_var = beta_hat_var, 
                      gamma = delta,
                      preds = Y_pred, 
                      PI_lower = c(Y_pred - var_Yhat_new_train_sigma2_c), 
                      PI_upper = Y_pred + var_Yhat_new_train_sigma2_c, 
                      random_var = b_var,
                      residual_var = sigma2_lmm, 
                      random_intercept = b_vec_int, 
                      random_slope = b_vec_slope)
      return(results)
}

