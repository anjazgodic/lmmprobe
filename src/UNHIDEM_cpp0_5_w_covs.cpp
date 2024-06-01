// [[Rcpp::depends(RcppArmadillo)]]

#include <fstream>
#include <iostream>
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
List UNHIDEM_cpp0_5_w_covs(const arma::vec y, const arma::mat Z, const arma::mat Vt, const arma::colvec X_adj, const arma::uvec rowind1, const arma::uvec rowind2, const arma::uvec rowind3, const arma::uvec colind1, const arma::uvec colind2, const arma::uvec colind3, const arma::colvec Wt, const arma::colvec t_Dt2, const arma::colvec t_Dt3, const arma::colvec Vt1_2, const arma::colvec Vt1_3, const arma::colvec Vt2_3, const arma::colvec W_var, const arma::colvec delta, const arma::colvec beta_vec, double sigma2_lmm) {
  
  // Getting the dimensions and initializing outputs
  int n = Z.n_rows, d = Z.n_cols;
  arma::mat coef_mat(d,6);
  arma::mat se(d,6,arma::fill::ones);
  arma::mat sig(d,1,arma::fill::ones);
  arma::mat like(d,1,arma::fill::ones);
  arma::mat W1(n,1,arma::fill::ones);
  
  for(int col=0;col<d; ++col){
    // Getting Z_m and p_m
    arma::colvec curr_Z = Z.col(col); 
    //arma::colvec curr_Z2 = Z2.col(col); 
    arma::colvec curr_delta = delta.row(col); 
    
    // Taking hypothesis m off of the expected mean and variance.
    arma::colvec t_Dt1 = Wt - curr_Z*curr_delta(0)*beta_vec.row(col);
    //arma::colvec t_Wt2 = W_var - curr_Z2*curr_delta(0)*(beta_var.row(col)+arma::square(beta_vec.row(col))*(1-curr_delta(0)));
    
    //Calculating some covariances
    arma::colvec E_b1W = Vt1_2 + (t_Dt1%t_Dt2);
    arma::colvec E_b2W = Vt1_3 + (t_Dt1%t_Dt3);
    arma::colvec E_b1b2 = Vt2_3 + (t_Dt2%t_Dt3);
    arma::colvec E_b1W_sum = arma::trans(W1)*E_b1W;
    arma::colvec E_b2W_sum = arma::trans(W1)*E_b2W;
    arma::colvec E_b1b2_sum = arma::trans(W1)*E_b1b2;

    //Creating a big matrix Dt_t_Dt and D2t so that we can do MLE
    arma::colvec t_Dt1_t_Dt1 = t_Dt1%t_Dt1;
    arma::colvec t_Dt1_t_Dt2 = t_Dt1%t_Dt2;
    arma::colvec t_Dt1_t_Dt3 = t_Dt1%t_Dt3;
    arma::colvec t_Dt2_t_Dt1 = t_Dt2%t_Dt1;
    arma::colvec t_Dt2_t_Dt2 = t_Dt2%t_Dt2;
    arma::colvec t_Dt2_t_Dt3 = t_Dt2%t_Dt3;
    arma::colvec t_Dt3_t_Dt1 = t_Dt3%t_Dt1;
    arma::colvec t_Dt3_t_Dt2 = t_Dt3%t_Dt2;
    arma::colvec t_Dt3_t_Dt3 = t_Dt3%t_Dt3;

    arma::mat Dt_t_Dt_1  = arma::join_rows(t_Dt1_t_Dt1,t_Dt1_t_Dt2,t_Dt1_t_Dt3);
    arma::mat Dt_t_Dt_2  = arma::join_rows(t_Dt2_t_Dt1,t_Dt2_t_Dt2,t_Dt2_t_Dt3);
    arma::mat Dt_t_Dt_3  = arma::join_rows(t_Dt3_t_Dt1,t_Dt3_t_Dt2,t_Dt3_t_Dt3);
    arma::mat Dt_t_Dt  = arma::join_cols(Dt_t_Dt_1, Dt_t_Dt_2, Dt_t_Dt_3);
    
    arma::mat D2t = Vt + Dt_t_Dt;

    // Calculating some sums for X'X.
    //t_Wt2 = t_Wt2 + arma::square(t_Wt);
    arma::colvec D2t_1N_1 = D2t.submat(rowind1,colind1);
    arma::colvec D2t_2N_2 = D2t.submat(rowind2,colind2);
    arma::colvec D2t_3N_3 = D2t.submat(rowind3,colind3);
    arma::colvec D2t_1N_1_sum = arma::trans(W1)*D2t_1N_1;
    arma::colvec D2t_2N_2_sum = arma::trans(W1)*D2t_2N_2;
    arma::colvec D2t_3N_3_sum = arma::trans(W1)*D2t_3N_3;
    //arma::colvec X2_adj_sum = arma::trans(W1)*X2_adj;

    // Making the X, the (X'X) expectation matrix.
    arma::mat bXt_1 = arma::join_rows(W1,curr_Z);
    arma::mat bXt_2 = arma::join_rows(bXt_1,t_Dt1,t_Dt2,t_Dt3);
    arma::mat bXt  = arma::join_rows(bXt_2,X_adj);
    arma::mat bXXt = arma::trans(bXt)*bXt;
    arma::mat bXXt2= bXXt;
    bXXt2(2,2) = D2t_1N_1_sum(0);
    bXXt2(2,3) = E_b1W_sum(0);
    bXXt2(2,4) = E_b2W_sum(0);
    bXXt2(3,2) = E_b1W_sum(0);
    bXXt2(3,3) = D2t_2N_2_sum(0);
    bXXt2(3,4) = E_b1b2_sum(0);
    bXXt2(4,2) = E_b2W_sum(0);
    bXXt2(4,3) = E_b1b2_sum(0);
    bXXt2(4,4) = D2t_3N_3_sum(0);
    //bXXt2(5,5) = X2_adj_sum(0);

    // Checking if (X'X) is invertable and calculating it's inverse
    if(bXXt2.is_sympd()){
    arma::mat XXt_inv = arma::inv_sympd(bXXt2);
      
    int k = bXt.n_cols;
    // Updating beta, and sigma^2
    arma::colvec t_beta = XXt_inv*arma::trans(bXt)*y;
    arma::colvec resid = y - bXt*t_beta;
    arma::colvec sigma2_update = arma::trans(resid)*resid/(n-k);
    arma::colvec like_update = -n*log(sigma2_update)/2;
    
    // Estimating the covariance of beta and their SE's
    arma::mat Vbt = sigma2_lmm*arma::diagvec(arma::trans(XXt_inv)*(bXXt2)*XXt_inv);
    arma::colvec stderrest = arma::sqrt(Vbt);
    
    // Exporting results.
    coef_mat.row(col) = arma::trans(t_beta);
    se.row(col) = arma::trans(stderrest);
    sig(col,0) = sigma2_update(0); 
    like(col,0) = like_update(0); 
    }
  }
  
  // Test statistics
  arma::mat T_vals = coef_mat/se;
  
  return List::create(Named("Coefficients") = coef_mat,Named("StdErr") = se, 
                      Named("Sigma") = sig, Named("T_statistics") = T_vals, Named("Log_like") = like);
}








