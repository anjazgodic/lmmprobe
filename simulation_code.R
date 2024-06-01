# Load libraries ####
rm(list=ls())
#Sys.setenv("PKG_CXFLAGS"="-std=c++11")
require(parallel)
library(MASS)
library(lme4)
library(tidyverse)
library(snowfall)
require(Rcpp)
require(RcppArmadillo)
# install.packages("RandomFieldsUtils_1.2.5.tar.gz", type = "source", dependencies = T)
# install.packages("RandomFields_3.3.14.tar.gz", type = "source", dependencies = T)
library(RandomFields)
RFoptions(spConform=FALSE)
library(lmmprobe)

##### Capture arguments passed from command line/bash script ##### 
# args <- commandArgs(trailingOnly = TRUE) #uncomment this if using this script from command line
args <- 1

##### Create simulation settings matrix ##### 
if(length(args)==1){
  settings1 <- expand.grid(
    n = c(250),
    M1 = c(0.05, 0.1),
    obs_per_sub_settings = c(4), 
    beta_settings = c(0.5, 0.75), 
    re_settings = c(1,2),
    b_var_number = c(1,2), 
    sigma2_lmm = c(100,150), 
    M = c(75^2),
    lat_sp = c(10),
    ep = 0.1,
    node_arg = 6,
    cv = c(1)
    )
  
  settings2 <- expand.grid(
    n = c(50),
    M1 = c(0.05, 0.1), 
    obs_per_sub_settings = c(3), 
    beta_settings = c(0.5, 0.75), 
    re_settings = c(1,2),
    b_var_number = c(1,2), 
    sigma2_lmm = c(100/10,150/10), 
    M = c(15^2),
    lat_sp = c(10), 
    ep = 0.1,
    node_arg = 6, 
    cv = c(1)
    )
  
  settings3 <- expand.grid(
    n = c(100), 
    M1 = c(0.05, 0.1), 
    obs_per_sub_settings = c(3), 
    beta_settings = c(0.5, 0.75), 
    re_settings = c(1,2),
    b_var_number = c(1,2), 
    sigma2_lmm = c(100/10,150/10), 
    M = c(25^2),
    lat_sp = c(10),
    ep = 0.05,
    node_arg = 6,
    cv = c(1)
    )
  
  settings <- rbind(settings1, settings2, settings3)
  rm(settings1, settings2, settings3)
  args <- as.numeric(settings[as.numeric(args),])
}

##### Saving current parameters ##### 
n = n_settings <- as.numeric(args[1])
M1_settings <- as.numeric(args[2])
obs_per_sub = obs_per_sub_settings <- as.numeric(args[3])
beta_value = beta_settings <- as.numeric(args[4])
number_re <- as.numeric(args[5])
b_var_number  <- as.numeric(args[6])
sigma2_lmm_value = sigma2_lmm_settings <- as.numeric(args[7])
M = M_settings <- as.numeric(args[8])

if(M_settings == 75^2){
  G_1_options <- c(50,100)
  G_2_options <- list(
    matrix(c(40, 0, 0, 25), nrow = 2, ncol = 2),
    matrix(c(60, 10, 10, 35), nrow = 2, ncol = 2)
    )
} else{
  G_1_options <- c(50/10,100/10)
  G_2_options <- list(
    matrix(c(40/10, 0, 0, 25/10), nrow = 2, ncol = 2), 
    matrix(c(60/10, 10/10, 10/10, 35/10), nrow = 2, ncol = 2)
    )
}

if(number_re == 1){
  b_var_value = b_var_settings <- G_1_options[b_var_number]
} else{
  b_var_value = b_var_settings <- G_2_options[[b_var_number]]
}

M1_value = M1_settings <- floor(M1_settings*M_settings)
lat_sp = lat_sp_settings <- as.numeric(args[9])
ep = as.numeric(args[10])
node_arg <- as.numeric(args[11])
cv_settings <- as.numeric(args[12])
  
##### Simulation set up and creating empty objects #####
reps <- 1:3
N <- n*obs_per_sub
maxit = 10000
alpha <- 0.05
one.sided = FALSE
B = 3
adj = 1
LR = 0
C = 1
lambda = 0.1
extra_obs <- obs_per_sub_settings
UseCores = node_arg-1

if (is.null(adj)) {
  adj <- 3/(1 + 1 * I(one.sided))
}
sigma_init = NULL

sigma2_simul <- 
  iter_simul <-  
  time_store <-  
  c()

mse_simul_test_xb <- 
  mse_simul_test <-  
  mad_simul_test_xb <-  
  mad_simul_test <- 
  c()

beta_true_simul <- 
  beta_hat_simul <- 
  beta_hat_var_simul <- 
  beta_hat_bias_simul <- 
  beta_hat_ecp_simul <- 
  p_k_hat_simul <- 
  matrix(NA, nrow = M-1, ncol = length(reps))
  
sum_delta <- 
  sum_delta_signal <- 
  c()

conf_int_ecp_simul_test <-  
  pred_int_ecp_simul_test <- 
  matrix(NA, nrow = N, ncol = length(reps))

conv_simul <-  
  conv_simul_W <- 
  matrix(NA, nrow = maxit+1, ncol = length(reps))

tr_xb_test <-  
  pred_xb_test <-  
  tr_test <-  
  pred_test <- 
  matrix(NA, nrow = N, ncol = length(reps))

random_effects_store <- 
var_pred_new_store <- 
  data.frame()

if(number_re == 1){
  sigma2_lmm_store <- 
    b_var_store <- 
    c()
  
  calibration_comparison_store <- c()
} else{
  sigma2_lmm_store <- 
    c()
  
  calibration_comparison_sigma2lmm_coefs_store <- 
    calibration_comparison_bvar_store <- 
    c()
  
  b_var_store <- 
    matrix(c(99,99,99,99), ncol = 2, nrow = 2)
}

##### Start simulation repetitions #####
for(simul_rep in reps){
      set.seed(20+1+simul_rep)
  
      start_time <- Sys.time()
      print(paste0("Simulation Number: ",simul_rep))
      
      ##### Generating training and testing data #####
      M <- M_settings
      
      ID <- rep(1:n, each = obs_per_sub)
      
      if(cv_settings == 1){
        ID_test <- rep(1:n, each = (obs_per_sub+extra_obs))
      } else{
        ID_test <- ID
      }
      
      tmp_ind <- rep(1:(obs_per_sub+extra_obs), n)
      
      ##### Generate X  #####
      x_grid <- seq(1, sqrt(M_settings), 1)
      
      print("Generating M predictors")
      if(cv_settings == 1){
        storage <- RFsimulate(
          model = RMstable(alpha = 2, scale = lat_sp, var = 1),
          x = x_grid,
          y = x_grid,
          grid = TRUE,
          n = n*(obs_per_sub+extra_obs)
          )
        
        x <- t(matrix(array(storage),M_settings,(n*(obs_per_sub+extra_obs))))
        
        x_test <- 
          storage_test <- 
          x[tmp_ind %in% c((obs_per_sub+1):(obs_per_sub+extra_obs)),]
        
        x <- x[tmp_ind %in% c(1:obs_per_sub),]
      
        } else{
        storage <- RFsimulate(
          model = RMstable(alpha = 2, scale = lat_sp, var = 1),
          x = x_grid,
          y = x_grid,
          grid = TRUE,
          n = N
          )
        
        x <- t(matrix(array(storage),M_settings,N))
        
        storage_test <- RFsimulate(
          model = RMstable(alpha = 2, scale = lat_sp, var = 1),
          x = x_grid,
          y = x_grid,
          grid = TRUE,
          n = N
          )
        
        x_test <- t(matrix(array(storage_test),M_settings,N))
      }

      if(cv_settings == 1){
          if(number_re == 2){
            x <- cbind(rep(1, nrow(x)), x)
            x_test <- cbind(rep(1, nrow(x_test)), x_test)
            x_adj <- matrix(rep(1:obs_per_sub, n)-1, ncol = 1)
            x_adj_test <- matrix(
              rep(c((obs_per_sub+1):(obs_per_sub+extra_obs)), n)-1,
              ncol = 1
              )
            x <- x[,-2]
            x_test <- x_test[,-2]
          } else{
            x[,1] <- rep(1, nrow(x))
            x_test[,1] <- rep(1, nrow(x_test))
            x_adj <- NULL
            x_adj_test <- NULL
          }
      } else{
          if(number_re == 2){
            x <- cbind(rep(1, nrow(x)), x)
            x_test <- cbind(rep(1, nrow(x_test)), x_test)
            x_adj <- matrix(rep(1:obs_per_sub, n)-1, ncol = 1) 
            x_adj_test <- matrix(rep(1:obs_per_sub, n)-1, ncol = 1) 
            x <- x[,-2]
            x_test <- x_test[,-2]
          } else{
            x[,1] <- rep(1, nrow(x))
            x_test[,1] <- rep(1, nrow(x_test))
            x_adj <- NULL
            x_adj_test <- NULL
          }
      }
      
      ##### Generate beta values #####
      j <- 1:M1_value 
      beta_true <- matrix(rep(0, M), ncol=1)
      beta_true[j, ] <- beta_value 
      beta_true_adj <- beta_value
      
      ##### Generate Random effects #####
      if(cv_settings == 1){
          re <- mvrnorm(n, mu=rep(0,number_re), Sigma=b_var_value)
          re <- matrix(rep(re, each = unique(table(ID_test))), ncol = ncol(re))
          re_test <- re[
            tmp_ind %in% c((obs_per_sub+1):(obs_per_sub+extra_obs)),
            ,
            drop=F
            ] 
          re <- re[tmp_ind %in% c(1:obs_per_sub),,drop=F]
      } else{
          re <- mvrnorm(n, mu=rep(0,number_re), Sigma=b_var_value)
          re <- matrix(rep(re, each = unique(table(ID))), ncol = ncol(re))
          re_test <- mvrnorm(n, mu=rep(0,number_re), Sigma=b_var_value)
          re_test <- matrix(
            rep(re_test, each = unique(table(ID))),
            ncol = ncol(re)
            )
      }
      
      ##### Generate Error #####
      if(cv_settings == 1){
        eps <- matrix(
          rnorm(
            (n*(obs_per_sub+extra_obs)),
            0, 
            sqrt(sigma2_lmm_value)
            ),
          ncol=1
          )
        
        eps_test <- eps[
          tmp_ind %in% c((obs_per_sub+1):(obs_per_sub+extra_obs)),
          ,
          drop=F
          ]
        
        eps <- eps[tmp_ind %in% c(1:obs_per_sub),,drop=F]
      
        } else{
        eps <- matrix(rnorm(N, 0, sqrt(sigma2_lmm_value)), ncol=1)
        eps_test <- matrix(rnorm(N, 0, sqrt(sigma2_lmm_value)), ncol=1)
      }
      
      ##### Generate Y and full dataset #####
      if(number_re == 1) {
        y <- x %*% beta_true + x[,1,drop=F]*(re) + eps
        y_test <- x_test %*% beta_true + 
          x_test[,1,drop=F]*(re_test) + 
          eps_test
      
        } else{
        y <- x %*% beta_true + 
          x_adj %*% beta_true_adj + 
          rowSums(cbind(x[,1],x_adj)*(re)) + 
          eps
        
        y_test <- x_test %*% beta_true + 
          x_adj_test %*% beta_true_adj + 
          rowSums(cbind(x_test[,1],x_adj_test)*(re_test)) + 
          eps_test
      }
      
      sim_data <- data.frame(
        cbind(
          ID,
          y,
          x_adj
          ,x
          )
        )
      
      if(cv_settings == 1){
        ID_test <- ID_test[
          tmp_ind %in% c((obs_per_sub+1):(obs_per_sub+extra_obs))
          ]
        }
      
      sim_data_test <- data.frame(
        cbind(
          ID_test,
          y_test,
          x_adj_test,
          x_test
          )
        )
      
      if(number_re == 1) {
        colnames(sim_data) <- 
          colnames(sim_data_test) <- 
          c("id","y",paste("x",1:length(beta_true),sep = ""))
      } else{
        colnames(sim_data) <- 
          colnames(sim_data_test) <- 
          c("id","y","x_adj",paste("x",1:length(beta_true),sep = ""))
      }
    
      ##### Setting initial values and initializing outputs #####
      ID <- sim_data[,"id"]
      Y <- sim_data[,"y"]
      Y_test = sim_data_test[,"y"]
      
      if(number_re == 1) {
        Z <- sim_data[,4:ncol(sim_data)]
        Z_test = sim_data_test[,4:ncol(sim_data_test)]
        X <- NULL
        X_test <- NULL
        RE_dat <- sim_data[,c("x1"), drop = F]
        RE_dat_test <- sim_data_test[,c("x1"), drop = F]
        V <- matrix(sim_data[,"id"], ncol=1)
        V_test <- matrix(sim_data_test[,"id"], ncol=1)
      } else{
        Z <- sim_data[,5:ncol(sim_data)]
        Z_test = sim_data_test[,5:ncol(sim_data_test)]
        X <- sim_data[,"x_adj", drop=F]
        X_test = sim_data_test[,"x_adj", drop=F]
        RE_dat <- sim_data[,c("x1", "x_adj"), drop = F]
        RE_dat_test <- sim_data_test[,c("x1", "x_adj"), drop = F]
        V <- as.matrix(RE_dat)
        V_test <- as.matrix(RE_dat_test)
      }
      
      ##### Call lmmprobe ####
      fit <- lmmprobe(
        Y = Y, 
        Z = Z, 
        V = V, 
        ID_data = ID,
        Y_test = Y_test, 
        Z_test = Z_test, 
        V_test = V_test, 
        ID_test,
        alpha = alpha,
        maxit = maxit,
        ep = ep,
        cpus = UseCores,
        B = B,
        adj = adj,
        LR = LR,
        C = C, 
        sigma_init = sigma_init
      )
      
      ##### Saving results in storage objects ##### 
      iter_simul[simul_rep] <- fit$count
      
      if(number_re == 1) {
        tr_xb_test[, simul_rep] <-
          as.numeric(
            beta_true[1] * fit$intercept + 
              as.matrix(Z_test) %*% beta_true[-1]
            )
        
        pred_xb_test[, simul_rep] <-
          (fit$c_coefs[1]) * fit$intercept +
          (as.matrix(Z_test) %*% (fit$beta * fit$gamma))
        
        tr_test[, simul_rep] <-
          as.numeric(
            beta_true[1] * fit$intercept + 
              as.matrix(Z_test) %*% beta_true[-1] + 
              re_test
            )
        
        if (cv_settings == 1) {
          pred_test[, simul_rep] <-
          (fit$c_coefs[1]) * fit$intercept +
          (as.matrix(Z_test) %*% (fit$beta * fit$gamma)) +
          (fit$c_coefs[3] * fit$b_vec_out[, 1])
        } else{
          pred_test[, simul_rep] <-
          (fit$c_coefs[1]) * fit$intercept +
          (as.matrix(Z_test) %*% (fit$beta * fit$gamma)) +
          (fit$c_coefs[3] * fit$b_vec_out_new[, 1])
        }
      } else{
        tr_xb_test[, simul_rep] <- as.numeric(
          beta_true[1] * fit$intercept +
            as.matrix(Z_test) %*% beta_true[-1] +
            as.matrix(X_test$x_adj) %*% beta_true_adj
        )
        
        pred_xb_test[, simul_rep] <- (
          (fit$c_coefs[1] * fit$intercept) +
            (
              as.matrix(Z_test) %*% (fit$beta * fit$gamma)
            ) +
            as.matrix(X_test$x_adj) %*% (fit$c_coefs[5])
        )
        
        tr_test[, simul_rep] <- as.numeric(
          beta_true[1] * fit$intercept +
            as.matrix(Z_test) %*% beta_true[-1] +
            as.matrix(X_test$x_adj) %*% beta_true_adj +
            RE_dat_test[, 1] * re_test[, 1] +
            RE_dat_test[, 2] * re_test[, 2]
        )
        
        if (cv_settings == 1) {
          pred_test[, simul_rep] <- (
                  (fit$c_coefs[1] * fit$intercept) +
                  (as.matrix(Z_test) %*% (fit$beta * fit$gamma)) +
                  (fit$c_coefs[3] * RE_dat_test[, 1] * fit$b_vec_out[, 1]) +
                  (fit$c_coefs[4] * RE_dat_test[, 2] * fit$b_vec_out[, 2]) +
                  as.matrix(X_test$x_adj) %*% (fit$c_coefs[5])
          )
        } else{
          pred_test[, simul_rep] <- (
                  (fit$c_coefs[1] * fit$intercept) +
                  (as.matrix(Z_test) %*% (fit$beta * fit$gamma)) +
                  (fit$c_coefs[3] * RE_dat_test[, 1] * fit$b_vec_out_new[, 1]) +
                  (fit$c_coefs[4] * RE_dat_test[, 2] * fit$b_vec_out_new[, 2]) +
                  as.matrix(X_test$x_adj) %*% (fit$c_coefs[5])
          )
        }
      }
      
      if (fit$try2 != 2) {
        mse_simul_test_xb[simul_rep] <-
          mean((tr_xb_test[, simul_rep] - pred_xb_test[, simul_rep]) ^ 2)
        
        mse_simul_test[simul_rep] <-
          mean((tr_test[, simul_rep] - pred_test[, simul_rep]) ^ 2)
        
        mad_simul_test_xb[simul_rep] <-
          mean(abs(tr_xb_test[, simul_rep] - pred_xb_test[, simul_rep]))
        
        mad_simul_test[simul_rep] <-
          mean(abs(tr_test[, simul_rep] - pred_test[, simul_rep]))
        
        beta_true_simul[, simul_rep] <- beta_true[-1]
        
        beta_hat_simul[, simul_rep] <- fit$beta
        
        p_k_hat_simul[, simul_rep] <- fit$gamma
        
        beta_hat_var_simul[, simul_rep] <- fit$beta_var
        
        beta_hat_bias_simul[, simul_rep] <-
          fit$beta - beta_true[-1]
        
        beta_hat_ecp_simul[, simul_rep] <-
          ifelse(
            beta_true[-1] >= fit$beta - 1.96 * sqrt(fit$beta_var) &
              beta_true[-1] <= fit$beta + 1.96 * sqrt(fit$beta_var),
            1,
            0
          )
        
        sum_delta[simul_rep] <- sum(fit$gamma) + 1
        
        sum_delta_signal[simul_rep] <-
          sum(fit$gamma[1:(M1_value - 1)]) + 1
        
        conf_int_ecp_simul_test[, simul_rep] <-
          ifelse((tr_test[, simul_rep] >=
                    (pred_test[, simul_rep] - 1.96 * sqrt(fit$var_Yhat_new))) &
                   (tr_test[, simul_rep] <=
                      (pred_test[, simul_rep] + 1.96 * sqrt(fit$var_Yhat_new))),
                 1,
                 0)
        
        pred_int_ecp_simul_test[, simul_rep] <-
          ifelse((tr_test[, simul_rep] >=
                    (
                      pred_test[, simul_rep] - 1.96 * sqrt(fit$var_Yhat_new_sigma2_c)
                    )) &
                   (tr_test[, simul_rep] <=
                      (
                        pred_test[, simul_rep] + 1.96 * sqrt(fit$var_Yhat_new_sigma2_c)
                      )),
                 1,
                 0)
        
        time_store[simul_rep] <- fit$time_unhidem
      } else {
        beta_true_simul[, simul_rep] <- beta_true[-1]
        
        mse_simul_test_xb[simul_rep] <-
          mse_simul_test[simul_rep] <-
          NA
        
        mad_simul_test_xb[simul_rep] <-
          mad_simul_test[simul_rep] <-
          NA
        
        beta_hat_simul[, simul_rep] <-
          beta_hat_simbeta_hat_var_simul[, simul_rep] <-
          beta_hat_bias_simul[, simul_rep] <-
          beta_hat_ecp_simul[, simul_rep] <-
          NA
        
        sum_delta[simul_rep] <- sum_delta_signal[simul_rep] <- NA
        
        conf_int_ecp_simul_test[, simul_rep] <-
          pred_int_ecp_simul_test[, simul_rep] <- NA
        
        time_store[simul_rep] <-
          time_unhidem
      }
      
      
      ##### Saving objects to compare convergence criteria #####
      conv_simul_W[1:length(fit$Xt_conv1_store_W), simul_rep] <-
        fit$Xt_conv1_store_W
      
      lower_ind <- 1
      upper_ind <-
        min(max(10000 - colSums(apply(
          conv_simul_W, 2, is.na
        ))) + 50, 10000)
      
      conv_simul_W_write <-
        cbind(lower_ind:upper_ind, conv_simul_W[lower_ind:upper_ind, ])
      
      ##### Preparing dataframes with results #####
      if (number_re == 1) {
        random_effects <-
          data.frame(
            cbind(
              rep(simul_rep, length(fit$b_vec)),
              fit$b_vec_int,
              re[!duplicated(re),]
              )
            )
        
        calibration_comparison <-
          data.frame(
            cbind(
              rep(simul_rep, length(fit$random_var)),
              fit$random_var,
              fit$residual_var,
              fit$c_coefs[1],
              fit$c_coefs[2],
              fit$c_coefs[3],
              fit$sigma2_c
            )
          )

        calibration_comparison_store <-
          data.frame(rbind(calibration_comparison_store, calibration_comparison))
        } else{
        
          random_effects <-
          data.frame(
            cbind(
              rep(simul_rep, length(fit$b_vec_int)),
              fit$b_vec_int,
              re[!duplicated(re), 1],
              fit$b_vec_slope,
              re[!duplicated(re), 2]
            )
          )
        calibration_comparison_sigma2lmm_coefs <-
          data.frame(
            cbind(
              rep(simul_rep, length(fit$residual_var)),
              fit$residual_var,
              fit$c_coefs[1],
              fit$c_coefs[2],
              fit$c_coefs[3],
              fit$c_coefs[4],
              fit$sigma2_c
            )
          )
        
          calibration_comparison_bvar <-
          data.frame(cbind(rep(simul_rep, nrow(fit$random_var)), 
                           fit$random_var))

        calibration_comparison_sigma2lmm_coefs_store <-
          data.frame(
            rbind(
              calibration_comparison_sigma2lmm_coefs_store,
              calibration_comparison_sigma2lmm_coefs
            )
          )

        calibration_comparison_bvar_store <-
          data.frame(rbind(
            calibration_comparison_bvar_store,
            calibration_comparison_bvar
          ))
      }
      
      random_effects_store <-
        rbind(random_effects_store, random_effects)
      
} #end of for loop for iterations


##### Preparing outputs #####
if (number_re == 1) {
  names(random_effects_store) <-
    c(
      "simul_rep",
      "b_vec_int",
      "re_int"
    )
  
  names(calibration_comparison_store) <-
    c(
      "simul_rep",
      "b_var",
      "sigma2_lmm",
      "mu_c",
      "alpha_c",
      "gamma_c",
      "sigma2_c"
    )
  } else{
  names(random_effects_store) <-
    c(
      "simul_rep",
      "b_vec_int",
      "re_int",
      "b_vec_slope",
      "re_slope"
    )
  
  names(calibration_comparison_sigma2lmm_coefs_store) <-
    c(
      "simul_rep",
      "sigma2_lmm",
      "mu_c",
      "alpha_c",
      "gamma_c",
      "omega_c",
      "sigma2_c"
    )

  names(calibration_comparison_bvar_store) <-
    c(
      "simul_rep",
      "b_var",
      "b_var"
      )
}

mse_mad <-
  data.frame(
    cbind(
      mse_simul_test_xb,
      mse_simul_test,
      mad_simul_test_xb,
      mad_simul_test,
      iter_simul,
      sum_delta,
      sum_delta_signal,
      time_store
    )
  )

betas_simul <-
  data.frame(
    cbind(
      beta_true_simul,
      beta_hat_simul, 
      p_k_hat_simul
    )
  )
names(betas_simul) <- c(
  paste(c("true_betas_simul"), 1:max(reps), sep = "_"),
  paste(c("lmmprobe_betas_simul"), 1:max(reps), sep = "_"), 
  paste(c("p_k_betas_simul"), 1:max(reps), sep = "_")
)

colnames(pred_int_ecp_simul_test) <- c(
  paste(c("PI_ECP_simul"), 1:max(reps), sep = "_")
)
