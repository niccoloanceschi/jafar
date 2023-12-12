
JAFAR_elbo_vb <- function(t_iter,y,X_m,n,M,p_m,K,K_Gm,
                          which_cusp='orig',
                          a_sig, b_sig,                # response noise
                          var_spike_theta,             # spike value in response loadings variances
                          a_xi, b_xi,                  # mixture weight in response loadings variances
                          a_m, b_m,                    # idiosyncratic noises
                          prec0, prec0m,               # "local" means
                          var_spike, var_spike_loc,    # spike value in loadings variances
                          alpha, alpha_loc,            # beta dist stick breaking
                          m_Theta, v_Theta, lD_Theta,
                          m_eta, v_eta, lD_eta,
                          a_s2_inv, b_s2_inv, m_mu_y, v_mu_y,
                          a_xi_vb, b_xi_vb,
                          logP_D, logP_Dm, logP_Dm_leq_Diag, 
                          logP_Zm, logP_Zm_leq,
                          a_rho_m, b_rho_m, 
                          m_phi_m, v_phi_m, lD_phi_m,
                          m_mu_m, v_mu_m, 
                          m_Loadings_m, v_Loadings_m, lD_Loadings_m,
                          a_s2_inv_m, b_s2_inv_m,
                          a_nu_m, b_nu_m, 
                          a_chi_vb=NULL, b_chi_vb=NULL,
                          a_chi_m=NULL, b_chi_m=NULL,
                          a_tau_m=NULL, b_tau_m=NULL,
                          var_slab_theta=NULL,var_slab=NULL,var_slab_loc=NULL){
  
  elbo = 0
  
  facTfac_m <- list()
  etaTeta <- crossprod(m_eta)
  for(m in c(1:M)){
    facTfac_m[[m]] <- matrix(0,K+K_Gm[m],K+K_Gm[m])
    etaTphi <- crossprod(m_eta,m_phi_m[[m]])
    phiTphi <- crossprod(m_phi_m[[m]])
    
    facTfac_m[[m]][c(1:K),c(1:K)]   <- n*v_eta + etaTeta
    facTfac_m[[m]][-c(1:K),-c(1:K)] <- n*v_phi_m[[m]] + phiTphi
    facTfac_m[[m]][c(1:K),-c(1:K)]  <- etaTphi
    facTfac_m[[m]][-c(1:K),c(1:K)]  <- t(etaTphi)
  }
  
  # STEP 0: Data ---------------------------------------------------------
  
  ## 0.a Response ------------------------------------------------------------
  elbo <- elbo - as.numeric( 0.5*sum(y^2) - m_mu_y*sum(y) + 0.5*n*(v_mu_y+m_mu_y^2) -
      sum((y-rep(m_mu_y,n))*(m_eta%*%m_Theta)) + 0.5*sum(diag((n*v_eta+etaTeta)%*%(v_Theta+tcrossprod(m_Theta)))) ) * a_s2_inv/b_s2_inv
  elbo <- elbo + 0.5*n*(-log(b_s2_inv)+digamma(a_s2_inv))
  
  ## 0.b Omics ------------------------------------------------------------
  for(m in c(1:M)){
    for(j in 1:p_m[m]){
      Load_Load_T <- v_Loadings_m[[m]][j,,] + tcrossprod(m_Loadings_m[[m]][j,])
      
      elbo <- elbo - as.numeric( 0.5*sum(X_m[[m]][,j]^2) -
          m_mu_m[[m]][j]*sum(X_m[[m]][,j]) + 0.5*n*(v_mu_m[[m]][j]+m_mu_m[[m]][j]^2) - 
          (X_m[[m]][,j]-rep(m_mu_m[[m]][j],n))%*%(cbind(m_eta,m_phi_m[[m]])%*%m_Loadings_m[[m]][j,]) +
          0.5*sum(diag(facTfac_m[[m]]%*%Load_Load_T)) ) * a_s2_inv_m[[m]][j]/b_s2_inv_m[[m]][j]
      elbo <- elbo + 0.5*n*(-log(b_s2_inv_m[[m]][j])+digamma(a_s2_inv_m[[m]][j]))
    }
  }
  
  # STEP 1: Loadings ---------------------------------------------------------
  
  ## 1.a Theta ---------------------------------------------------------------
  
  logP_Dm_leq_Sum <- colSums(logP_Dm_leq_Diag)
  
  prob_chi <- exp(logP_D[2,])*(1-exp(logP_Dm_leq_Sum))
  
  if(which_cusp=='orig'){
    
    E_chi <- prob_chi*(a_chi_vb/b_chi_vb) + (1-prob_chi)*rep(1/var_spike_theta,K)
    
    elbo <- elbo + 0.5*sum(prob_chi*(-log(b_chi_vb)+digamma(a_chi_vb)) +
                             (1-prob_chi)*log(1/var_spike_theta))
  } else if (which_cusp=='mod'){
    E_chi <- prob_chi*rep(1/var_slab_theta,K) + (1-prob_chi)*rep(1/var_spike_theta,K)
    
    elbo <- elbo + 0.5*sum(prob_chi*log(1/var_slab_theta) + (1-prob_chi)*log(1/var_spike_theta))
  
  } else {
    stop("Choose a valid CUSP version ('orig' or 'mod')")
  }
  
  elbo <- elbo - 0.5*(m_Theta^2 + diag(v_Theta))%*%E_chi
  
  elbo <- elbo - 0.5*lD_Theta
    
  ## 1.b Lambda_m, Gamma_m ---------------------------------------------------
  
  for(m in c(1:M)){
    
    prob_chi_m <- (1-exp(logP_Dm_leq_Diag[m,]))*(1-exp(logP_D[1,]+colSums(logP_Dm_leq_Diag[-m,])))
    prob_tau_m <- (1-exp(diag(logP_Zm_leq[[m]])))
    
    if(which_cusp=='orig'){
      E_chi_m    <- prob_chi_m*a_chi_m[m,]/b_chi_m[m,] + (1-prob_chi_m)*rep(1/var_spike[m],K)
      E_tau_m    <- prob_tau_m*a_tau_m[[m]]/b_tau_m[[m]] + (1-prob_tau_m)*rep(1/var_spike_loc[m],K_Gm[m])
    
      elbo <- elbo + 0.5*p_m[m]*sum(prob_chi_m*(-log(b_chi_m[m,])+digamma(a_chi_m[m,])) +
                                    (1-prob_chi_m)*log(1/var_spike[m]))
      elbo <- elbo + 0.5*p_m[m]*sum(prob_tau_m*(-log(b_tau_m[[m]])+digamma(a_tau_m[[m]])) +
                                    (1-prob_tau_m)*log(1/var_spike_loc[m]))
    } else if (which_cusp=='mod'){
      E_chi_m    <- prob_chi_m*rep(1/var_slab[m],K) + (1-prob_chi_m)*rep(1/var_spike[m],K)
      E_tau_m    <- prob_tau_m*rep(1/var_slab_loc[m],K_Gm[m]) + (1-prob_tau_m)*rep(1/var_spike_loc[m],K_Gm[m])
    
      elbo <- elbo + 0.5*p_m[m]*sum(prob_chi_m*log(1/var_slab[m]) + (1-prob_chi_m)*log(1/var_spike[m]))
      elbo <- elbo + 0.5*p_m[m]*sum(prob_tau_m*log(1/var_slab_loc[m]) + (1-prob_tau_m)*log(1/var_spike_loc[m]))
    
    } else {
      stop("Choose a valid CUSP version ('orig' or 'mod')")
    }
      
    for(j in c(1:p_m[m])){
      
      elbo <- elbo - 0.5*(m_Loadings_m[[m]][j,]^2 + diag(v_Loadings_m[[m]][j,,]))%*%c(E_chi_m,E_tau_m)
      
      elbo <- elbo - 0.5*lD_Loadings_m[[m]][j]
    }
  }
  
  # STEP 2: Intercepts -------------------------------------------------------
  
  ## 2.a mu_y ----------------------------------------------------------------
  elbo <- elbo - 0.5*(m_mu_y^2 + v_mu_y)*prec0
  
  elbo <- elbo + 0.5*log(v_mu_y)
  
  ## 2.b mu_m ----------------------------------------------------------------
  for(m in 1:M){
    elbo <- elbo - 0.5*sum(m_mu_m[[m]]^2 + v_mu_m[[m]])*prec0m[m]
    
    elbo <- elbo + 0.5*sum(log(v_mu_m[[m]]))
  }
  
  # STEP 3: Idiosyncratic components -----------------------------------------
  
  ## 3.a s2_inv --------------------------------------------------------------
  elbo <- elbo + (a_sig-1)*(-log(b_s2_inv)+digamma(a_s2_inv)) - b_sig*a_s2_inv/b_s2_inv
  
  elbo <- elbo + a_s2_inv - log(b_s2_inv) + lgamma(a_s2_inv) + (1-a_s2_inv)*digamma(a_s2_inv)
  
  ## 3.b s2_inv_m ------------------------------------------------------------
  for(m in 1:M){
    elbo <- elbo + sum( (a_m[m]-1)*(-log(b_s2_inv_m[[m]])+digamma(a_s2_inv_m[[m]])) - b_m[m]*a_s2_inv_m[[m]]/b_s2_inv_m[[m]] )
  
    elbo <- elbo + sum( a_s2_inv_m[[m]] - log(b_s2_inv_m[[m]]) + lgamma(a_s2_inv_m[[m]]) + (1-a_s2_inv_m[[m]])*digamma(a_s2_inv_m[[m]]) )
    }
  
  # STEP 4: Factors ----------------------------------------------------------

  ## 4.a eta -----------------------------------------------------------------
  elbo <- elbo - 0.5*( sum(m_eta^2) - n*sum(diag(v_eta)) )
                      
  elbo <- elbo - 0.5*n*lD_eta
  
  ## 4.b phi_m ---------------------------------------------------------------
  for(m in 1:M){
    elbo <- elbo - 0.5*( sum(m_phi_m[[m]]^2) - n*sum(diag(v_phi_m[[m]])) )
    
    elbo <- elbo - 0.5*n*lD_phi_m[[m]]
  }
  
  # STEP 5: Latent indicatirs ------------------------------------------------
  
  ## 5.a Shared Component ----------------------------------------------------
  
  ### 5.a.1 delta ------------------------------------------------------------
  
  elbo <- elbo + sum(exp(logP_D[1,]))*digamma(a_xi_vb) + sum(exp(logP_D[2,]))*digamma(b_xi_vb)
  
  elbo <- elbo - sum(exp(logP_D)*logP_D)
  
  ### 5.a.2 delta_m ----------------------------------------------------------
  
  for(m in 1:M){
    
    a_digamma <- c(digamma(a_rho_m[m,]) - digamma(a_rho_m[m,]+b_rho_m[m,]),0)
    b_digamma <- c(0,digamma(b_rho_m[m,]) - digamma(a_rho_m[m,]+b_rho_m[m,]))
    
    elbo <- elbo + sum(exp(logP_Dm[m,,])*matrix(a_digamma+b_digamma,K,K))
    
    elbo <- elbo - sum(exp(logP_Dm[m,,])*logP_Dm[m,,])
  }
  
  ## 5.b Specific Components ---------------------------------------------------
  
  ### 5.b.1 z_m (latent indicators) ----
  for(m in 1:M){
    
    a_digamma <- c(digamma(a_nu_m[[m]]) - digamma(a_nu_m[[m]]+b_nu_m[[m]]),0)
    b_digamma <- c(0,cumsum(digamma(b_nu_m[[m]]) - digamma(a_nu_m[[m]]+b_nu_m[[m]])))
    
    elbo <- elbo + sum(exp(logP_Zm[[m]])*matrix(a_digamma+b_digamma,K_Gm[m],K_Gm[m]))
    
    elbo <- elbo - sum(exp(logP_Zm[[m]])*logP_Zm[[m]])
  }
    
  # STEP 6: Stick Breaking Elements ------------------------------------------
  
  ## 5.a xi (spike and slab mixture weight for Theta) ----
  
  elbo <- elbo + (a_xi-1)*(digamma(a_xi_vb)-digamma(a_xi_vb+b_xi_vb)) +
    (b_xi-1)*(digamma(b_xi_vb)-digamma(a_xi_vb+b_xi_vb))
  
  elbo <- elbo + lgamma(a_xi_vb) + lgamma(b_xi_vb) - lgamma(a_xi_vb+b_xi_vb) -
    (a_xi_vb-1)*digamma(a_xi_vb) - (b_xi_vb-1)*digamma(b_xi_vb) -
    (a_xi_vb+b_xi_vb-2)*digamma(a_xi_vb+b_xi_vb) 
    
  
  ## 6.b rho_m (stick breaking elements for Lambda_m) ----
  
  for(m in 1:M){
    elbo <- elbo + (alpha[m]-1)*sum(digamma(b_rho_m[m,])+digamma(a_rho_m[m,]+b_rho_m[m,]))
    
    elbo <- elbo + sum( lgamma(a_rho_m[m,]) + lgamma(b_rho_m[m,]) - lgamma(a_rho_m[m,]+b_rho_m[m,]) -
      (a_rho_m[m,]-1)*digamma(a_rho_m[m,]) - (b_rho_m[m,]-1)*digamma(b_rho_m[m,]) -
      (a_rho_m[m,]+b_rho_m[m,]-2)*digamma(a_rho_m[m,]+b_rho_m[m,]) )
  }
      
  ## 6.c nu_m (stick breaking elements for Gamma_m) ----
  
  for(m in 1:M){
    elbo <- elbo + (alpha_loc[m]-1)*sum(digamma(b_nu_m[[m]])+digamma(a_nu_m[[m]]+b_nu_m[[m]]))
    
    elbo <- elbo + sum( lgamma(a_nu_m[[m]]) + lgamma(b_nu_m[[m]]) - lgamma(a_nu_m[[m]]+b_nu_m[[m]]) -
      (a_nu_m[[m]]-1)*digamma(a_nu_m[[m]]) - (b_nu_m[[m]]-1)*digamma(b_nu_m[[m]]) -
      (a_nu_m[[m]]+b_nu_m[[m]]-2)*digamma(a_nu_m[[m]]+b_nu_m[[m]]) )
  }
    
  return(elbo)
  
}

