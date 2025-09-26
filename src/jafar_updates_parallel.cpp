#include <RcppArmadillo.h>
#include <RcppParallel.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]

using namespace Rcpp;
using namespace arma;
using namespace RcppParallel;

struct UpdateLoadingsWorker : public Worker {
  
  const int n;
  const int Ktot;
  const arma::mat& X_m;
  const arma::mat& facTfac;
  const arma::mat& fac;
  const arma::vec& mu_m;
  const arma::vec& s2_inv_m;
  const arma::mat prior_prec_mj;
  
  arma::mat& new_load_m;

  UpdateLoadingsWorker(int n_in, int Ktot_in,
                       const arma::mat& X_m_in, const arma::mat& facTfac_in,
                       const arma::mat& fac_in, const arma::vec& mu_m_in,
                       const arma::vec& s2_inv_m_in, const arma::vec& prec_m_in,
                       arma::mat& new_load_m_out) :
    n(n_in), Ktot(Ktot_in),
    X_m(X_m_in), facTfac(facTfac_in), fac(fac_in),
    mu_m(mu_m_in), s2_inv_m(s2_inv_m_in),
    prior_prec_mj(diagmat(prec_m_in)),
    new_load_m(new_load_m_out) {}

  void operator()(std::size_t begin, std::size_t end) {

    for (std::size_t j = begin; j < end; ++j) {
    
      arma::mat Q_load_mj = prior_prec_mj + s2_inv_m(j) * facTfac;
      arma::vec r_load_mj = s2_inv_m(j) * fac.t() * (X_m.col(j) - mu_m(j));
     
      arma::mat L_load_mj = trimatu(chol(Q_load_mj));
      arma::vec Lr_load_mj = solve(trimatl(L_load_mj.t()), r_load_mj);
      
      arma::vec mean_load_mj = solve(L_load_mj, Lr_load_mj);
      arma::vec std_load_mj = solve(L_load_mj, arma::randn<arma::vec>(Ktot));
     
      new_load_m.row(j) = (mean_load_mj + std_load_mj).t();
    }
  } 
};

// [[Rcpp::export]]
Rcpp::NumericMatrix update_loadings_parallel(int n, int p_m, int Ktot,
                                             const arma::mat& X_m, const arma::mat& facTfac,
                                             const arma::mat& fac, const arma::vec& mu_m,
                                             const arma::vec& s2_inv_m, const arma::vec& prec_m) {

  arma::mat new_load_m(p_m, Ktot);

  UpdateLoadingsWorker worker(n, Ktot, X_m, facTfac, fac, mu_m, s2_inv_m, prec_m, new_load_m);

  parallelFor(0, p_m, worker);

  return Rcpp::wrap(new_load_m);
}