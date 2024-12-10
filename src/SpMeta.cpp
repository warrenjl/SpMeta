#include "RcppArmadillo.h"
#include "SpMeta.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List SpMeta(int mcmc_samples,
                  Rcpp::List theta_hat,
                  Rcpp::List se,
                  Rcpp::List x,
                  int model_indicator,
                  Rcpp::Nullable<Rcpp::List> neighbors = R_NilValue,
                  Rcpp::Nullable<arma::vec> metrop_var_rho_trans = R_NilValue,
                  Rcpp::Nullable<double> a_sigma2_prior = R_NilValue,
                  Rcpp::Nullable<double> b_sigma2_prior = R_NilValue,
                  Rcpp::Nullable<double> sigma2_beta_prior = R_NilValue,
                  Rcpp::Nullable<double> a_tau2_prior = R_NilValue,
                  Rcpp::Nullable<double> b_tau2_prior = R_NilValue,
                  Rcpp::Nullable<arma::vec> sigma2_init = R_NilValue,
                  Rcpp::Nullable<Rcpp::NumericVector> beta_init = R_NilValue,
                  Rcpp::Nullable<Rcpp::List> phi_init = R_NilValue,
                  Rcpp::Nullable<arma::vec> tau2_init = R_NilValue,
                  Rcpp::Nullable<arma::vec> rho_init = R_NilValue){
 
//Defining Parameters and Quantities of Interest
int p_x = Rcpp::as<arma::mat>(x[0]).n_cols;

int n = theta_hat.size();
arma::vec m(n); m.fill(0.00);
for(int j = 0; j < n; ++j){
   m(j) = Rcpp::as<arma::vec>(theta_hat[j]).size();
   }

Rcpp::List W(n);
if(neighbors.isNotNull()){ 
  W = Rcpp::as<Rcpp::List>(neighbors);
  }

Rcpp::List Omega_inv(n);
for(int j = 0; j < n; ++j){
  
   arma::mat temp_mat(m(j), m(j)); temp_mat.fill(0.00);
   for(int k = 0; k < m(j); ++ k){
      temp_mat(k,k) = (1.00/pow(Rcpp::as<arma::vec>(se[j])(k), 2)); 
      }
   Omega_inv[j] = temp_mat;
   
   }
  
Rcpp::List theta_true(mcmc_samples);
for(int j = 0; j < mcmc_samples; ++j){
  
   Rcpp::List temp_list(n); 
   theta_true[j] = temp_list;
   
   for(int k = 0; k < n; ++k){
  
      arma::vec temp_vec(m(k)); temp_vec.fill(0.00);
      Rcpp::as<Rcpp::List>(theta_true[j])[k] = temp_vec;
   
      }
   
   }
arma::mat sigma2(n, mcmc_samples); sigma2.fill(0.00);
arma::mat beta(p_x, mcmc_samples); beta.fill(0.00);
Rcpp::List phi(mcmc_samples);
for(int j = 0; j < mcmc_samples; ++j){
  
   Rcpp::List temp_list(n); 
   phi[j] = temp_list;
   
   for(int k = 0; k < n; ++k){
   
      arma::vec temp_vec(m(k)); temp_vec.fill(0.00);
      Rcpp::as<Rcpp::List>(phi[j])[k] = temp_vec;
  
      }
   
   } 
arma::mat tau2(n, mcmc_samples); tau2.fill(0.00);
arma::mat rho(n, mcmc_samples); rho.fill(0.00);

//Prior Information
double a_sigma2 = 0.01;
if(a_sigma2_prior.isNotNull()){
  a_sigma2 = Rcpp::as<double>(a_sigma2_prior);
  }

double b_sigma2 = 0.01;
if(b_sigma2_prior.isNotNull()){
  b_sigma2 = Rcpp::as<double>(b_sigma2_prior);
  }

double sigma2_beta = 10000.00;
if(sigma2_beta_prior.isNotNull()){
  sigma2_beta = Rcpp::as<double>(sigma2_beta_prior);
  }

double a_tau2 = 0.01;
if(a_tau2_prior.isNotNull()){
  a_tau2 = Rcpp::as<double>(a_tau2_prior);
  }
  
double b_tau2 = 0.01;
if(b_tau2_prior.isNotNull()){
  b_tau2 = Rcpp::as<double>(b_tau2_prior);
  }

//Initial Values
for(int j = 0; j < n; ++j){
   Rcpp::as<Rcpp::List>(theta_true[0])[j] = Rcpp::as<arma::vec>(theta_hat[j]);
   }

sigma2.col(0).fill(1.00);
if(sigma2_init.isNotNull()){
  sigma2.col(0) = Rcpp::as<arma::vec>(sigma2_init);
  }

beta.col(0).fill(0.00);
if(beta_init.isNotNull()){
  beta.col(0) = Rcpp::as<arma::vec>(beta_init);
  }

for(int j = 0; j < n; ++j){
  
   arma::vec temp_vec(m(j)); temp_vec.fill(0.00);
   Rcpp::as<Rcpp::List>(phi[0])[j] = temp_vec;
   if(phi_init.isNotNull()){
     Rcpp::as<Rcpp::List>(phi[0])[j] = Rcpp::as<arma::vec>(Rcpp::as<Rcpp::List>(phi_init)[j]);
     }
  
   }

tau2.col(0).fill(1.00);
if(tau2_init.isNotNull()){
  tau2.col(0) = Rcpp::as<arma::vec>(tau2_init);
  }

rho.col(0).fill(0.50);
if(rho_init.isNotNull()){
  rho.col(0) = Rcpp::as<arma::vec>(rho_init);
  }

Rcpp::List Q(n);
for(int j = 0; j < n; ++j){
   
   arma::mat temp_mat(m[j], m[j]); temp_mat.fill(0.00);
   Q[j] = temp_mat;
   
   }
arma::vec Q_log_deter(n); Q_log_deter.fill(0.00);
if(model_indicator > 0){
  
  double sign = 0.00;
  for(int j = 0; j < n; ++j){
    
     arma::mat W_mat = Rcpp::as<arma::mat>(W[j]);
     Q[j] = rho(j,0)*(arma::diagmat(arma::sum(W_mat, 1)) - W_mat) +
            (1.00 - rho(j,0))*arma::eye(m(j), m(j));
     log_det(Q_log_deter(j), sign, Rcpp::as<arma::mat>(Q[j]));
     
     }
  
  }

//Main Sampling Loop:  Nonspatial
if(model_indicator == 0){
  
  Rcpp::List phi_fixed(n);
  for(int j = 0; j < n; ++j){
     
     arma::vec temp_vec(m(j)); temp_vec.fill(0.00);
     phi_fixed[j] = temp_vec;
  
     }
  
  for(int j = 1; j < mcmc_samples; ++j){
    
     //theta_true
     theta_true[j] = theta_true_update(theta_hat,
                                       Omega_inv,
                                       n,
                                       m,
                                       x, 
                                       sigma2.col(j-1),
                                       beta.col(j-1),
                                       phi_fixed);
    
     //sigma2 Update
     sigma2.col(j) = sigma2_update(n,
                                   m,
                                   x,
                                   a_sigma2,
                                   b_sigma2,
                                   theta_true[j],
                                   beta.col(j-1),
                                   phi_fixed);
    
     //beta Update
     beta.col(j) = beta_update(n,
                               p_x,
                               x, 
                               sigma2_beta,
                               theta_true[j],
                               sigma2.col(j),
                               phi_fixed);
    
     //Progress
     if((j + 1) % 10 == 0){ 
       Rcpp::checkUserInterrupt();
       }
    
     if(((j + 1) % int(round(mcmc_samples*0.10)) == 0)){
      
       double completion = round(100*((j + 1)/(double)mcmc_samples));
       Rcpp::Rcout << "Progress: " << completion << "%" << std::endl;
      
       Rcpp::Rcout << "**************" << std::endl;
      
       }
    
     }
  
  return Rcpp::List::create(Rcpp::Named("theta_true") = theta_true,
                            Rcpp::Named("sigma2") = sigma2,
                            Rcpp::Named("beta") = beta);
  
  }

//Main Sampling Loop:  Spatial without Nugget
if(model_indicator == 1){
  
  //Metropolis Settings
  arma::vec acctot_rho_trans(n); acctot_rho_trans.fill(0);
  arma::vec metrop_sd_rho_trans(n); metrop_sd_rho_trans.fill(0.00);
  if(metrop_var_rho_trans.isNotNull()){
    metrop_sd_rho_trans = sqrt(Rcpp::as<arma::vec>(metrop_var_rho_trans));
    }
  
  for(int j = 1; j < mcmc_samples; ++j){
    
     //theta_true
     theta_true[j] = theta_true_1_update(theta_hat,
                                         Omega_inv,
                                         n,
                                         m,
                                         x, 
                                         beta.col(j-1),
                                         tau2.col(j-1),
                                         Q);
    
     //beta Update
     beta.col(j) = beta_1_update(n,
                                 p_x,
                                 x, 
                                 sigma2_beta,
                                 theta_true[j],
                                 tau2.col(j-1),
                                 Q);
     
     //tau2 Update
     tau2.col(j) = tau2_1_update(n,
                                 m,
                                 x,
                                 a_tau2,
                                 b_tau2,
                                 theta_true[j],
                                 beta.col(j),
                                 Q);
     
     //rho Update
     Rcpp::List rho_1_output = rho_1_update(n,
                                            m,
                                            x,
                                            W,
                                            theta_true[j],
                                            beta.col(j),
                                            tau2.col(j),
                                            rho.col(j-1),
                                            Q,
                                            Q_log_deter,
                                            metrop_sd_rho_trans,
                                            acctot_rho_trans);
     
     rho.col(j) = Rcpp::as<arma::vec>(rho_1_output[0]);
     acctot_rho_trans = Rcpp::as<arma::vec>(rho_1_output[1]);
     Q = Rcpp::as<Rcpp::List>(rho_1_output[2]);
     Q_log_deter = Rcpp::as<arma::vec>(rho_1_output[3]);
    
     //Progress
     if((j + 1) % 10 == 0){ 
       Rcpp::checkUserInterrupt();
       }
    
     if(((j + 1) % int(round(mcmc_samples*0.10)) == 0)){
      
       double completion = round(100*((j + 1)/(double)mcmc_samples));
       Rcpp::Rcout << "Progress: " << completion << "%" << std::endl;
      
       double accrate_rho_trans_min = round(100*(min(acctot_rho_trans)/(double)j));
       Rcpp::Rcout << "rho Acceptance (min): " << accrate_rho_trans_min << "%" << std::endl;
      
       double accrate_rho_trans_max = round(100*(max(acctot_rho_trans)/(double)j));
       Rcpp::Rcout << "rho Acceptance (max): " << accrate_rho_trans_max << "%" << std::endl;
      
       Rcpp::Rcout << "*************************" << std::endl;
      
       }
    
     }
  
  return Rcpp::List::create(Rcpp::Named("theta_true") = theta_true,
                            Rcpp::Named("beta") = beta,
                            Rcpp::Named("phi") = phi,
                            Rcpp::Named("tau2") = tau2,
                            Rcpp::Named("rho") = rho,
                            Rcpp::Named("acctot_rho_trans") = acctot_rho_trans);
  
  }

//Main Sampling Loop:  Spatial with Nugget
if(model_indicator == 2){
  
  //Metropolis Settings
  arma::vec acctot_rho_trans(n); acctot_rho_trans.fill(0);
  arma::vec metrop_sd_rho_trans(n); metrop_sd_rho_trans.fill(0.00);
  if(metrop_var_rho_trans.isNotNull()){
    metrop_sd_rho_trans = sqrt(Rcpp::as<arma::vec>(metrop_var_rho_trans));
    }
  
  for(int j = 1; j < mcmc_samples; ++j){
   
     //theta_true
     theta_true[j] = theta_true_update(theta_hat,
                                       Omega_inv,
                                       n,
                                       m,
                                       x, 
                                       sigma2.col(j-1),
                                       beta.col(j-1),
                                       phi[j-1]);
     
     //sigma2 Update
     sigma2.col(j) = sigma2_update(n,
                                   m,
                                   x,
                                   a_sigma2,
                                   b_sigma2,
                                   theta_true[j],
                                   beta.col(j-1),
                                   phi[j-1]);
      
     //beta Update
     beta.col(j) = beta_update(n,
                               p_x,
                               x, 
                               sigma2_beta,
                               theta_true[j],
                               sigma2.col(j),
                               phi[j-1]);
      
     //phi Update
     phi[j] = phi_update(n,
                         m,
                         x, 
                         theta_true[j],
                         sigma2.col(j),
                         beta.col(j),
                         tau2.col(j-1),
                         Q);
      
     //tau2 Update
     tau2.col(j) = tau2_update(n,
                               m,
                               a_tau2,
                               b_tau2,
                               phi[j],
                               Q);
    
     //rho Update
     Rcpp::List rho_output = rho_update(n,
                                        m,
                                        W,
                                        phi[j],
                                        tau2.col(j),
                                        rho.col(j-1),
                                        Q,
                                        Q_log_deter,
                                        metrop_sd_rho_trans,
                                        acctot_rho_trans);
      
     rho.col(j) = Rcpp::as<arma::vec>(rho_output[0]);
     acctot_rho_trans = Rcpp::as<arma::vec>(rho_output[1]);
     Q = Rcpp::as<Rcpp::List>(rho_output[2]);
     Q_log_deter = Rcpp::as<arma::vec>(rho_output[3]);
      
     //Progress
     if((j + 1) % 10 == 0){ 
       Rcpp::checkUserInterrupt();
       }
  
     if(((j + 1) % int(round(mcmc_samples*0.10)) == 0)){
    
       double completion = round(100*((j + 1)/(double)mcmc_samples));
       Rcpp::Rcout << "Progress: " << completion << "%" << std::endl;
     
       double accrate_rho_trans_min = round(100*(min(acctot_rho_trans)/(double)j));
       Rcpp::Rcout << "rho Acceptance (min): " << accrate_rho_trans_min << "%" << std::endl;
     
       double accrate_rho_trans_max = round(100*(max(acctot_rho_trans)/(double)j));
       Rcpp::Rcout << "rho Acceptance (max): " << accrate_rho_trans_max << "%" << std::endl;
     
       Rcpp::Rcout << "*************************" << std::endl;
    
       }
  
     }
                                  
  return Rcpp::List::create(Rcpp::Named("theta_true") = theta_true,
                            Rcpp::Named("sigma2") = sigma2,
                            Rcpp::Named("beta") = beta,
                            Rcpp::Named("phi") = phi,
                            Rcpp::Named("tau2") = tau2,
                            Rcpp::Named("rho") = rho,
                            Rcpp::Named("acctot_rho_trans") = acctot_rho_trans);
   
  }

return R_NilValue;

}

