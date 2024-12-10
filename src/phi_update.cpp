#include "RcppArmadillo.h"
#include "SpMeta.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List phi_update(int n,
                      arma::vec m,
                      Rcpp::List x, 
                      Rcpp::List theta_true,
                      arma::vec sigma2,
                      arma::vec beta,
                      arma::vec tau2_old,
                      Rcpp::List Q){
  
Rcpp::List phi(n);
for(int j = 0; j < n; ++j){
    
   arma::mat cov_phi = inv_sympd(eye(m(j), m(j))/sigma2(j) + 
                                 Rcpp::as<arma::mat>(Q[j])/tau2_old(j));
   arma::vec mean_phi = cov_phi*(Rcpp::as<arma::vec>(theta_true[j]) - Rcpp::as<arma::mat>(x[j])*beta)/sigma2(j);
    
   arma::mat ind_norms = arma::randn(1, m(j));
   arma::vec phi_temp = mean_phi + 
                        trans(ind_norms*arma::chol(cov_phi));
    
   phi[j] = phi_temp - 
            mean(phi_temp);
    
   }

return(phi);

}



