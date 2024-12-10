#include "RcppArmadillo.h"
#include "SpMeta.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec beta_1_update(int n,
                        int p_x,
                        Rcpp::List x, 
                        double sigma2_beta,
                        Rcpp::List theta_true,
                        arma::vec tau2_old,
                        Rcpp::List Q){
  
arma::mat cov_beta(p_x, p_x); cov_beta.fill(0.00);
arma::vec mean_beta(p_x); mean_beta.fill(0.00);
for(int j = 0; j < n; ++j){
  
   arma::mat Q_mat = Rcpp::as<arma::mat>(Q[j]);
   arma::mat x_mat = Rcpp::as<arma::mat>(x[j]);
   arma::mat x_mat_trans = trans(x_mat);
   cov_beta = cov_beta + 
              x_mat_trans*(Q_mat*x_mat)/tau2_old(j);
   mean_beta = mean_beta +
               x_mat_trans*(Q_mat*Rcpp::as<arma::vec>(theta_true[j]))/tau2_old(j);
  
   }
cov_beta = inv_sympd(cov_beta + 
                     (1.00/sigma2_beta)*eye(p_x, p_x));
mean_beta = cov_beta*mean_beta;

arma::mat ind_norms = arma::randn(1, p_x);
arma::vec beta = mean_beta + 
                 trans(ind_norms*arma::chol(cov_beta));

return(beta);

}
    

