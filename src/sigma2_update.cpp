#include "RcppArmadillo.h"
#include "SpMeta.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec sigma2_update(int n,
                        arma::vec m,
                        Rcpp::List x,
                        double a_sigma2,
                        double b_sigma2,
                        Rcpp::List theta_true,
                        arma::vec beta_old,
                        Rcpp::List phi_old){

arma::vec sigma2(n); sigma2.fill(0.00);
for(int j = 0; j < n; ++j){
    
   double a_sigma2_update = m(j)/2.00 +
                            a_sigma2;
   arma::vec temp_vec = Rcpp::as<arma::vec>(theta_true[j]) + 
                        -Rcpp::as<arma::mat>(x[j])*beta_old +
                        -Rcpp::as<arma::vec>(phi_old[j]);
   double b_sigma2_update = 0.50*dot(temp_vec, temp_vec) +
                            b_sigma2;
   sigma2(j) = 1.00/R::rgamma(a_sigma2_update,
                              (1.00/b_sigma2_update));
   
   }

return(sigma2);

}





