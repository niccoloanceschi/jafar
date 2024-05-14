
#' Update of view-wise loadings within the Gibbs sampler for JAFAR under the D-CUSP prior
#'
#' @param n Number of observations (integer)
#' @param p_m Dimension of the considered view (integer)
#' @param K Number of shared factors
#' @param K_Gm Number of view-specific factors
#' @param X_m Observed predictors in the considered view (real matrix: n x p_m)
#' @param facTfac: t-crossp-roduct of overall factors (real matrix: (K+K_Gm) x (K+K_Gm))
#' @param eta Shared latent factors (real matrix: n x K)
#' @param phi_m View-specific latent factors (real matrix: n x K_Gm)
#' @param mu_m: Intercepts (real vector: p_m)
#' @param s2_inv_m: Idiosyncrstic components precision (real positive vector: p_m)
#' @param chi_m Shared-loadings prior variances (real positive vector: K)
#' @param tau_m View-specific loadings prior variances (real positive vector :K_Gm)
#' @return Updated loadings
#' 
update_loadings <- function(n, p_m, K, K_Gm, 
                            X_m, facTfac, eta,  phi_m,
                            mu_m, s2_inv_m, chi_m, tau_m){

  new_Loadings_m <- matrix(0,p_m,K+K_Gm)
  
  prior_prec_mj <- diag(c(chi_m,tau_m),K+K_Gm,K+K_Gm)
  eta_phi_m     <- cbind(eta,phi_m)
  
  for(j in c(1:p_m)){
    
    Q_Loadings_mj <- prior_prec_mj+s2_inv_m[j]*facTfac
    r_Loadings_mj <- s2_inv_m[j]*crossprod(eta_phi_m,X_m[,j]-rep(mu_m[j],n))
    
    L_Loadings_mj  <- chol(Q_Loadings_mj)
    Lr_Loadings_mj <- forwardsolve(t(L_Loadings_mj), r_Loadings_mj)
    
    mean_Loadings_mj <- backsolve(L_Loadings_mj, Lr_Loadings_mj)
    std_Loadings_mj  <- backsolve(L_Loadings_mj, rnorm(K+K_Gm))
    
    new_Loadings_m[j,] <- t(mean_Loadings_mj + std_Loadings_mj)
  }
  
  return(new_Loadings_m)
}