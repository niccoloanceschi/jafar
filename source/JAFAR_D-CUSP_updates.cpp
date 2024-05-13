#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::NumericMatrix update_loadings(int n, int p_m, int K, int K_Gm, 
                          const arma::mat& X_m, const arma::mat& facTfac,
                          const arma::mat& eta, const arma::mat& phi_m,
                          const arma::vec& mu_m, const arma::vec& s2_inv_m, 
                          const arma::vec& chi_m, const arma::vec& tau_m) {
  
  arma::mat new_load_m(p_m, K+K_Gm);
  
  arma::mat prior_prec_mj = diagmat(join_cols(chi_m,tau_m));
  arma::mat eta_phi_m_T   = (join_rows(eta,phi_m)).t(); 
  
  for (int j = 0; j < p_m; ++j) {

    arma::mat Q_load_mj = prior_prec_mj + s2_inv_m(j)*facTfac;
    arma::vec r_load_mj = s2_inv_m(j) * eta_phi_m_T * (X_m.col(j) - mu_m(j)*arma::ones(n));

    arma::mat L_load_mj  = trimatu(chol(Q_load_mj));
    arma::vec Lr_load_mj = solve(trimatl(L_load_mj.t()), r_load_mj);

    arma::vec mean_load_mj = solve(L_load_mj, Lr_load_mj);
    arma::vec std_load_mj  = solve(L_load_mj, arma::randn<arma::vec>(K+K_Gm));

    new_load_m.row(j) = (mean_load_mj + std_load_mj).t();
  }
  
  return Rcpp::wrap(new_load_m);
}

  

