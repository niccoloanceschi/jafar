
# # INPUTS:
# #   n: integer
# #   p_m: integer
# #   K: integer
# #   K_Gm: integer
# #   X_m: real matrix (n x p_m)
# #   mu_m: real vector (p_m)
# #   s2_inv_m: real vector (p_m)
# #   eta: real matrix (n x K)
# #   phi_m: real matrix (n x K_Gm)
# #   chi_m: real vector (K)
# #   tau_m: real vector (K_Gm)
# #   facTfac: real matrix ((K+K_Gm) x (K+K_Gm))
# # EXTRA:
# #   noise_loadings: common random noise for numerical comparison: real vector (K+K_Gm)

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