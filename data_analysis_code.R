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
library(groupdata2)
library(abind)
library(umap)
library(expm)
library(lmmprobe)

#### Reading in dataset ####
data(SLE)

#### Create folds for cross-validation ####
set.seed(222)
k = 5
n_obs <- nrow(real_data)
real_data$id <- factor(real_data$id)
real_data$ind <- 1:nrow(real_data)

df_folded <- fold(real_data, k = k, id_col = "id")

train_ind <- test_ind <- 
  test_ind_id <- train_ind_id  <- vector("list", k)

#### Create indexes for cross-validation type 1 ####
for(i in 1:k){
  test_ind[[i]] <- pull(df_folded[df_folded$.folds == i, ], "ind")
  test_ind_id[[i]] <- pull(df_folded[df_folded$.folds == i, ], "id")
  
  train_ind[[i]] <- pull(df_folded[df_folded$.folds != i, ], "ind")
  train_ind_id[[i]] <- pull(df_folded[df_folded$.folds != i, ], "id")
}
real_data$ind <- NULL

#### Create empty objects ####
resid_marg_xb <- resid_cond <- res_cov <- 
  matrix(NA, nrow = nrow(real_data), ncol = k)

tr_xb_train <- tr_train <- tr_test <- 
  matrix(NA, nrow = nrow(real_data), ncol = k)

pred_xb_train <- pred_train <- pred_test <- 
  matrix(NA, nrow = nrow(real_data), ncol = k)

mspe_test = mad_test <- c()
iter_store <- time_store <- c()
random_effects_store <-  c()
sigma2_lmm_store <- sigma2_c_store <- b_var_store <- c()
mu_store <- alpha_store <- gamma_store <- c()
variables_selected_store <- c()

#### Start cross-validated analyses ####
for(i in 1:k){
  
  cv_fold <- i
  print(paste0("This is CV fold ", i))
  
  #### Creating training and testing datasets ####
  data <- real_data[train_ind[[i]],]
  ID_data <- as.numeric(as.character(train_ind_id[[i]]))
  
  data_test <- real_data[test_ind[[i]],]
  ID_test <- as.numeric(as.character(test_ind_id[[i]]))
  
  Z <- data[,4:ncol(data)]
  Z_test = data_test[,4:ncol(data_test)]
  
  Y <- matrix(data[,"y"], ncol=1)
  Y_test <- matrix(data_test[,"y"], ncol=1)
  
  V <- matrix(data[,"id"], ncol=1)
  V_test <- matrix(data_test[,"id"], ncol=1)
  
  #### Train lmmprobe ####
  train_lmmprobe <- lmmprobe(
    Y = Y,Z = Z, V = V, ID_data = ID_data,
    Y_test = Y_test, Z_test = Z_test, V_test = V_test, ID_test = ID_test, 
    adj = 1, sigma_init = 1, ep = 0.1, B = 3
  )
  
  #### This is the calculation various quantities for a test set, W_ast_new, b_vec_new, etc. #####
  tr_test[1:nrow(data_test),cv_fold] <- as.matrix( 
    data_test[,"y",drop = F]
  )
  
  b_vec_out_test <- unique(train_lmmprobe$b_vec_out_new_train_test)
  
  bhat_test_test <- merge(
    cbind(ID_test, train_lmmprobe$RE_dat_test),
    b_vec_out_test,
    by.x = "ID_test", by.y="sort(ID_test)",
    all.x = T
  )
  
  pred_test[1:length(ID_test),cv_fold] <-
    (train_lmmprobe$c_coefs[1]*train_lmmprobe$intercept_test) + 
    (as.matrix(Z_test)%*%(train_lmmprobe$beta*train_lmmprobe$gamma)) + 
    (train_lmmprobe$c_coefs[3]*
       as.numeric(bhat_test_test[,"intercept"])*
       bhat_test_test[,"b_vec_int_new"]
     )
  
  #### End test calculations  ####
  
  #### Saving outputs LMM-PROBE ####
  iter_store[cv_fold] <- train_lmmprobe$count
  
  #### Saving true Y, Predictions, Residuals ####
  tr_train[1:nrow(data),cv_fold] <- as.matrix(data[,"y",drop = F])
  tr_xb_train[1:nrow(data),cv_fold] <- as.matrix(data[,"y",drop = F])
  
  pred_train[1:nrow(data),cv_fold] <-
    (train_lmmprobe$c_coefs[1]*train_lmmprobe$intercept) +
    (as.matrix(Z)%*%(train_lmmprobe$beta*train_lmmprobe$gamma)) +
    (train_lmmprobe$c_coefs[3]*
       train_lmmprobe$RE_dat[,1]*
       train_lmmprobe$b_vec_out[,1])
  
  pred_xb_train[1:nrow(data),cv_fold] <-
    (train_lmmprobe$c_coefs[1]*train_lmmprobe$intercept) +
    (as.matrix(Z)%*%(train_lmmprobe$beta*train_lmmprobe$gamma))
  
   V = sapply(
    simplify = F,
    1:length(unique(ID_data)),
    function(x) {
      solve(
        sqrtm(
          (train_lmmprobe$curr_T_split[[x]] %*% 
             train_lmmprobe$random_var %*% 
             train_lmmprobe$curr_T_split_t[[x]] + 
             train_lmmprobe$R[[x]] * train_lmmprobe$residual_var
           )
          )
        )
      }
    )
  
  sfStop()
  
  mspe_test[cv_fold] <- mean(
    (tr_test[,cv_fold]-pred_test[,cv_fold])^2,
    na.rm = T
    )
  
  mad_test[cv_fold] <- mean(
    abs(tr_test[, cv_fold] - pred_test[, cv_fold]),
    na.rm = T
    )
  
  time_store[cv_fold] <- train_lmmprobe$time_unhidem
  
  resid_marg_xb[, cv_fold] <-
    tr_xb_train[, cv_fold] - pred_xb_train[, cv_fold]
  
  resid_cond[, cv_fold] <-
    (tr_train[, cv_fold] - pred_train[, cv_fold])
  
  resid_split <- lapply(
    split(
      na.omit(resid_cond[, cv_fold]),
      ID_data
      ),
    as.matrix
    )
  
  I_RR <- sapply(
    1:length(unique(ID_data)),
    function(x) {
      train_lmmprobe$R[[x]] - 
        (V[[x]] %*% resid_split[[x]]) %*% t(V[[x]] %*% resid_split[[x]])
      }
    )
  
  res_cov[1:length(unique(ID_data)), cv_fold] <- sapply(
    1:length(unique(ID_data)),
    function(x) {
      sqrt(sum(I_RR[[x]] ^ 2))
      }
    )
  
  mu_store <- c(mu_store, train_lmmprobe$c_coefs[1])
  alpha_store <- c(alpha_store, train_lmmprobe$c_coefs[2])
  gamma_store <- c(gamma_store, train_lmmprobe$c_coefs[3])
  sigma2_c_store <- c(sigma2_c_store, train_lmmprobe$sigma2_c)
  sigma2_lmm_store <- c(sigma2_lmm_store, train_lmmprobe$residual_var)
  
  variables_selected_store <- c(
    variables_selected_store, 
    sum(train_lmmprobe$gamma > 0.5)
  )
  
  random_effects <- data.frame(
    cbind(
      rep(cv_fold, length(train_lmmprobe$b_vec_int)), 
      train_lmmprobe$b_vec_int, 
      train_lmmprobe$b_vec_slope
    )
  )
  random_effects_store <- rbind(random_effects_store, random_effects)
  
} # end of cross validation


#### Writing analysis outputs ####
results <- data.frame(
  cbind(
    fold_rep = 1:k,
    sigma2_lmm = sigma2_lmm_store,
    mu = mu_store,
    alpha = alpha_store,
    gamma = gamma_store,
    sigma2_c = sigma2_c_store,
    b_var = b_var_store,
    MSPE = mspe_test,
    MAD = mad_test,
    `Number of Variables Selected` = variables_selected_store,
    `Number of iterations` = iter_store,
    `Time by iteration` = time_store
  )
)

write.csv(
  results,
  file="data_analysis_results.csv", 
  row.names = F
)


#### Writing Random Effects ####
names(random_effects_store) <- c("fold_rep", "b_vec_int")
write.csv(
  random_effects_store,
  file="data_analysis_random_effects.csv",
  row.names = F
)

#### Writing Residuals ####
resids_analysis <- data.frame(
  resid_marg_xb,
  pred_xb_train,
  resid_cond,
  pred_train,
  res_cov
  )
names(resids_analysis) <- c(
  paste(c("res_marg_xb"), 1:k, sep="_"), 
  paste(c("pred_marg_xb"), 1:k, sep="_"), 
  paste(c("res_cond"), 1:k, sep="_"), 
  paste(c("pred_cond"), 1:k, sep="_"), 
  paste(c("res_cov"), 1:k, sep="_")
)
write.csv(
  resids_analysis,
  file="data_analysis_residuals.csv", 
  row.names = F
)



