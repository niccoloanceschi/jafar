
gibbs_JAFAR_CUSP_init <- function(y, X_m, n, M, p_m,          # input data
                                  K0, K0_m,                   # initial number of factors
                                  a_sig, b_sig,               # response noise
                                  a_theta, b_theta,           # slab in response loadings variances
                                  var_spike_theta,            # spike value in response loadings variances
                                  a_xi, b_xi,                 # mixture weight in response loadings variances
                                  a_m, b_m,                   # idiosyncratic noises
                                  prec0, prec0m,              # "local" means
                                  var_spike, var_spike_vb,    # spike value in loadings variances
                                  a_chi, b_chi,               # slab in 'shared' loadings variances
                                  a_tau, b_tau,               # slab in 'local' loadings variances
                                  alpha, alpha_loc,           # beta dist stick breaking
                                  iterMax, rel_thr,
                                  seed) {
  
  # random initialization ----
  
  param_init <- within(list(), {
      
    K <- K0
    
    Theta   <- rep(NA,K)
    eta     <- matrix(rnorm(n*K),n,K)
    
    s2_inv  <- rgamma(1,shape=a_sig,rate=b_sig)
    mu_y    <- 0.
    
    chi      <- rep(1,K)
    delta    <- rep(1,K)
    xi       <- 0.5
    
    active_T <- c(1:K) # active elements in Theta
    
    K_T_eff <- K # NUMBER of active elements in Theta
    K_Lm_eff <- rep(K,M) # NUMBER of active columns in Lambda_m
    
    K_Gm     <- K0_m      # n. of retained columns in Gamma_m
    K_Gm_eff <- K0_m      # NUMBER of active columns in Gamma_m
    
    active_L <- active_G <- list() # active columns among retained ones
    
    Lambda_m <- Gamma_m <- phi_m <- s2_inv_m <- mu_m <- list()
    z_m <- nu_m <- w_m <- tau_m <- list()
    
    delta_m <- matrix(K,M,K) 
    
    rho_m   <- matrix(NA,M,K)
    xi_m    <- matrix(1/K,M,K) 
    chi_m   <- matrix(1,M,K)
    
    for(m in c(1:M)){
      active_L[[m]] <- c(1:K)
      active_G[[m]] <- c(1:K_Gm[m])
      
      Lambda_m[[m]] <- matrix(NA,p_m[m],K)
      Gamma_m[[m]]  <- matrix(NA,p_m[m],K_Gm[m])
      phi_m[[m]]    <- matrix(rnorm(K_Gm[m]*n),n,K_Gm[m])
      
      s2_inv_m[[m]] <- rgamma(p_m[m],shape=a_m[m],rate=b_m[m])
      mu_m[[m]]     <- rep(0,p_m[m])
      
      z_m[[m]]   <- rep(K_Gm[m],K_Gm[m])
      nu_m[[m]]  <- rep(NA,K_Gm[m])
      w_m[[m]]   <- rep(1/K_Gm[m],K_Gm[m])
      tau_m[[m]] <- rep(1,K_Gm[m])
    }  
  })
  
  return(param_init)
}

