#include "RcppArmadillo.h"
#include "SpMeta.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List rho_1_update(int n,
                        arma::vec m,
                        Rcpp::List x,
                        Rcpp::List W,
                        Rcpp::List theta_true,
                        arma::vec beta,
                        arma::vec tau2,
                        arma::vec rho_old,
                        Rcpp::List Q,
                        arma::vec Q_log_deter,
                        arma::vec metrop_sd_rho_trans,
                        arma::vec acctot_rho_trans){

arma::vec rho = rho_old;
for(int j = 0; j < n; ++j){
  
   double sign = 0.00;
   arma::vec temp_vec = Rcpp::as<arma::vec>(theta_true[j]) +
                        -Rcpp::as<arma::mat>(x[j])*beta;
   arma::mat W_mat = Rcpp::as<arma::mat>(W[j]);
   
   arma::mat Q_old = Rcpp::as<arma::mat>(Q[j]);
   double Q_log_deter_old = Q_log_deter(j);
   double rho_old = rho(j);
   double rho_transform_old = log(rho(j)/(1.00 - rho(j)));
   double denom = 0.50*Q_log_deter_old +
                  -0.50*dot(temp_vec, (Q_old*temp_vec))/tau2(j) +
                  rho_transform_old +
                  -2.00*log(1.00 + exp(rho_transform_old));
      
   double rho_transform = R::rnorm(rho_transform_old, 
                                   metrop_sd_rho_trans(j));
   rho(j) = 1.00/(1.00 + exp(-rho_transform));
   
   Q[j] = rho(j)*(arma::diagmat(arma::sum(W_mat, 1)) - W_mat) +
          (1.00 - rho(j))*arma::eye(m(j), m(j));
   log_det(Q_log_deter(j), sign, Rcpp::as<arma::mat>(Q[j]));
   
   double numer = 0.50*Q_log_deter(j) +
                  -0.50*dot(temp_vec, (Rcpp::as<arma::mat>(Q[j])*temp_vec))/tau2(j) +
                  rho_transform +
                  -2.00*log(1.00 + exp(rho_transform));  
      
   double ratio = exp(numer - denom);   
   double acc = 1;
   if(ratio < R::runif(0.00, 1.00)){
        
     Q[j] = Q_old;
     Q_log_deter(j) = Q_log_deter_old;
     rho(j) = rho_old;
     acc = 0;
      
     }
    
   acctot_rho_trans(j) = acctot_rho_trans(j) + 
                         acc;
      
   }

return Rcpp::List::create(Rcpp::Named("rho") = rho,
                          Rcpp::Named("acctot_rho_trans") = acctot_rho_trans,
                          Rcpp::Named("Q") = Q,
                          Rcpp::Named("Q_log_deter") = Q_log_deter);

}



