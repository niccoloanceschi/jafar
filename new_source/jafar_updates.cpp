#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::NumericMatrix update_loadings(int n, int p_m, int Ktot, 
                          const arma::mat& X_m, const arma::mat& facTfac,
                          const arma::mat& fac, const arma::vec& mu_m,
                          const arma::vec& s2_inv_m, const arma::vec& prec_m) {
  
  arma::mat new_load_m(p_m, Ktot);
  
  arma::mat prior_prec_mj = diagmat(prec_m);
  
  for (int j = 0; j < p_m; ++j) {

    arma::mat Q_load_mj = prior_prec_mj + s2_inv_m(j)*facTfac;
    arma::vec r_load_mj = s2_inv_m(j) * fac.t() * (X_m.col(j) - mu_m(j));

    arma::mat L_load_mj  = trimatu(chol(Q_load_mj));
    arma::vec Lr_load_mj = solve(trimatl(L_load_mj.t()), r_load_mj);

    arma::vec mean_load_mj = solve(L_load_mj, Lr_load_mj);
    arma::vec std_load_mj  = solve(L_load_mj, arma::randn<arma::vec>(Ktot));

    new_load_m.row(j) = (mean_load_mj + std_load_mj).t();
  }
  
  return Rcpp::wrap(new_load_m);
}
