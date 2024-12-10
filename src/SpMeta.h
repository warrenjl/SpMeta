#ifndef __SpMeta__
#define __SpMeta__

Rcpp::List theta_true_update(Rcpp::List theta_hat,
                             Rcpp::List Omega_inv,
                             int n,
                             arma::vec m,
                             Rcpp::List x, 
                             arma::vec sigma2_old,
                             arma::vec beta_old,
                             Rcpp::List phi_old);

Rcpp::List theta_true_1_update(Rcpp::List theta_hat,
                               Rcpp::List Omega_inv,
                               int n,
                               arma::vec m,
                               Rcpp::List x, 
                               arma::vec beta_old,
                               arma::vec tau2_old,
                               Rcpp::List Q);
  
arma::vec sigma2_update(int n,
                        arma::vec m,
                        Rcpp::List x,
                        double a_sigma2_prior,
                        double b_sigma2_prior,
                        Rcpp::List theta_true,
                        arma::vec beta_old,
                        Rcpp::List phi_old);
  
arma::vec beta_update(int n,
                      int p_x,
                      Rcpp::List x, 
                      double sigma2_beta_prior,
                      Rcpp::List theta_true,
                      arma::vec sigma2,
                      Rcpp::List phi_old);

arma::vec beta_1_update(int n,
                        int p_x,
                        Rcpp::List x, 
                        double sigma2_beta,
                        Rcpp::List theta_true,
                        arma::vec tau2_old,
                        Rcpp::List Q);

Rcpp::List phi_update(int n,
                      arma::vec m,
                      Rcpp::List x, 
                      Rcpp::List theta_true,
                      arma::vec sigma2,
                      arma::vec beta,
                      arma::vec tau2_old,
                      Rcpp::List Q);

arma::vec tau2_update(int n,
                      arma::vec m,
                      double a_tau2_prior,
                      double b_tau2_prior,
                      Rcpp::List phi,
                      Rcpp::List Q);

arma::vec tau2_1_update(int n,
                        arma::vec m,
                        Rcpp::List x,
                        double a_tau2,
                        double b_tau2,
                        Rcpp::List theta_true,
                        arma::vec beta,
                        Rcpp::List Q);

Rcpp::List rho_update(int n,
                      arma::vec m,
                      Rcpp::List W,
                      Rcpp::List phi,
                      arma::vec tau2,
                      arma::vec rho_old,
                      Rcpp::List Q,
                      arma::vec Q_log_deter,
                      arma::vec metrop_sd_rho_trans,
                      arma::vec acctot_rho_trans);

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
                        arma::vec acctot_rho_trans);

Rcpp::List SpMeta(int mcmc_samples,
                  Rcpp::List theta_hat,
                  Rcpp::List se,
                  Rcpp::List x,
                  int model_indicator,
                  Rcpp::Nullable<Rcpp::List> neighbors,
                  Rcpp::Nullable<arma::vec> metrop_var_rho_trans,
                  Rcpp::Nullable<double> a_sigma2_prior,
                  Rcpp::Nullable<double> b_sigma2_prior,
                  Rcpp::Nullable<double> sigma2_beta_prior,
                  Rcpp::Nullable<double> a_tau2_prior,
                  Rcpp::Nullable<double> b_tau2_prior,
                  Rcpp::Nullable<arma::vec> sigma2_init,
                  Rcpp::Nullable<Rcpp::NumericVector> beta_init,
                  Rcpp::Nullable<Rcpp::List> phi_init,
                  Rcpp::Nullable<arma::vec> tau2_init,
                  Rcpp::Nullable<arma::vec> rho_init);

#endif // __SpMeta__
