#include "RcppArmadillo.h"
#include "SpMeta.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List theta_true_1_update(Rcpp::List theta_hat,
                               Rcpp::List Omega_inv,
                               int n,
                               arma::vec m,
                               Rcpp::List x, 
                               arma::vec beta_old,
                               arma::vec tau2_old,
                               Rcpp::List Q){
  
Rcpp::List theta_true(n);
for(int j = 0; j < n; ++j){
    
   arma::mat Q_mat = Rcpp::as<arma::mat>(Q[j]);
   arma::mat Omega_inv_mat = Rcpp::as<arma::mat>(Omega_inv[j]);
   arma::mat cov_theta_true = inv_sympd(Omega_inv_mat + 
                                        (1.00/tau2_old(j))*Q_mat);
   arma::mat mean_theta_true = cov_theta_true*(Omega_inv_mat*Rcpp::as<arma::vec>(theta_hat[j]) + Q_mat*Rcpp::as<arma::mat>(x[j])*beta_old/tau2_old(j));
   arma::mat ind_norms = arma::randn(1, m(j));
   arma::vec theta_true_vec = mean_theta_true + 
                              trans(ind_norms*arma::chol(cov_theta_true));
   theta_true[j] = theta_true_vec;   
    
   }

return(theta_true);

}



