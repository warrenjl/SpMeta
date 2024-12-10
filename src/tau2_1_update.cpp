#include "RcppArmadillo.h"
#include "SpMeta.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec tau2_1_update(int n,
                        arma::vec m,
                        Rcpp::List x,
                        double a_tau2,
                        double b_tau2,
                        Rcpp::List theta_true,
                        arma::vec beta,
                        Rcpp::List Q){

arma::vec tau2(n); tau2.fill(0.00);
for(int j = 0; j < n; ++j){
    
   double a_tau2_update = m(j)/2.00 +
                          a_tau2;
   arma::vec temp_vec = Rcpp::as<arma::vec>(theta_true[j]) +
                        - Rcpp::as<arma::mat>(x[j])*beta;
   double b_tau2_update = 0.50*dot(temp_vec, (Rcpp::as<arma::mat>(Q[j])*temp_vec)) +
                          b_tau2;
   tau2(j) = 1.00/R::rgamma(a_tau2_update,
                            (1.00/b_tau2_update));
   
   }

return(tau2);

}



