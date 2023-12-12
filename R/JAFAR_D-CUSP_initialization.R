
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
                                  seed,
                                  vb_init, vb_rank_only) {
  
  # vb initialization ----
  
  if(vb_init){
    
    print("Inialization via VB")
    ris_VB <- vb_JAFAR_CUSP(y, X_m, n, M, p_m,                 # input data
                            iterMax=iterMax,                   # optimization hyper-params
                            rel_thr=rel_thr, seed=seed,        # optimization hyper-params
                            K=K0, K_Gm=K0_m,                   # latent dim bounds
                            a_sig=a_sig, b_sig=b_sig,          # response noise
                            a_theta=a_theta, b_theta=b_theta,  # slab in response loadings variances
                            var_spike_theta=var_spike_theta,   # spike value in response loadings variances
                            a_xi=a_xi, b_xi=b_xi,              # mixture weight in response loadings variances
                            a_m=a_m, b_m=b_m,                  # idiosyncratic noises
                            prec0=prec0, prec0m=prec0m,        # "local" means
                            var_spike=var_spike_vb,            # spike value in loadings variances
                            a_chi=a_chi, b_chi=b_chi,          # slab in 'shared' loadings variances
                            a_tau=a_tau, b_tau=b_tau,          # slab in 'local' loadings variances
                            alpha=alpha, alpha_loc=alpha_loc)  # beta dist stick breaking
    print("VB completed")
    
    act_J <- ris_VB$P_J_active>0
    K_J <- sum(act_J)
    
    if(K_J > 1){
      
      logP_Dm_leq_Diag <- t(apply(ris_VB$logP_Dm_leq, 1, diag))
      logP_Zm_leq_Diag <- lapply(ris_VB$logP_Zm_leq, diag)
      
      act_T_J <- (exp(ris_VB$logP_D[2,])>0 & act_J)
      act_L_J <- lapply(1:M, function(m) ((1-exp(logP_Dm_leq_Diag[m,])>0) & act_J))
      act_L   <- lapply(1:M, function(m) ((1-exp(logP_Dm_leq_Diag[m,])>0) & !act_J))
      act_G   <- lapply(1:M, function(m) ((1-exp(logP_Zm_leq_Diag[[m]])>0)))
      
      if(!vb_rank_only){
        
        param_init <- within(list(), {
        
          K    <- K_J+1
          K_Gm <- sapply(act_G,sum)+sapply(act_L,sum)+rep(1,M)
        
          Theta   <- c(ris_VB$m_Theta[act_J],sqrt(var_spike_theta)*rnorm(1))
          eta     <- cbind(ris_VB$m_eta[,act_J],rnorm(n))
          
          s2_inv  <- ris_VB$b_s2_inv/(ris_VB$a_s2_inv-1)
          mu_y    <- ris_VB$m_mu_y
          
          chi <- rep(1/var_spike_theta,K_J+1)
          
          delta   <- rep(0,K_J+1)
          
          xi      <- ris_VB$a_xi_vb / (ris_VB$a_xi_vb+ris_VB$b_xi_vb)
          
          s2_inv_m <- lapply(1:M, function(m) ris_VB$b_s2_inv_m[[m]]/(ris_VB$a_s2_inv_m[[m]]-1))
          mu_m     <- ris_VB$m_mu_m
          
          K_T_eff <- sum(act_T_J) # NUMBER of active elements in Theta
          K_Lm_eff <- rep(K_J,M) # NUMBER of active columns in Lambda_m
          
          K_Gm_eff <- K0_m-1      # NUMBER of active columns in Gamma_m
          
          active_T <- c(1:K_J)[act_T_J[act_J]]
          active_L <- lapply(1:M, function(m) c(1:K_J)[act_L_J[[m]][act_J]])
          active_G <- lapply(1:M, function(m) c(1:K_Gm_eff[m]))
          
        })
        
        param_init$chi[c(act_T_J[act_J],F)] <- ris_VB$a_chi_vb/ris_VB$b_chi_vb[act_T_J]
        param_init$delta[c(act_T_J[act_J],F)] <- 1
        
        Lambda_m <- lapply(1:M, function(m) matrix(rnorm(p_m[m]),p_m[m],1))
        Gamma_m  <- lapply(1:M, function(m) matrix(rnorm(p_m[m]),p_m[m],1))
        phi_m    <- lapply(1:M, function(m) matrix(rnorm(n),n,1)) 
        
        delta_m <- matrix(1,M,K_J+1)
        rho_m   <- matrix(1,M,K_J+1)
        xi_m    <- matrix(NA,M,K_J+1) 
        chi_m   <- matrix(1/var_spike,M,K_J+1)
        
        z_m   <- as.list(rep(c(1),M))
        nu_m  <- as.list(rep(c(1),M))
        w_m   <- as.list(rep(c(1),M))
        tau_m <- as.list(rep(c(1/var_spike),M))
        
        for(m in 1:M){
          
          delta_m[m,1:K_J] <- round(colSums(exp(ris_VB$logP_Dm[m,,])*matrix(1:K0,K0,K0)))[act_J]
          
          rho_m[m,1:K_J]   <- ris_VB$a_rho_m[m,act_J[-K0]] / (ris_VB$a_rho_m[m,act_J[-K0]] + ris_VB$b_rho_m[m,act_J[-K0]])
          
          xi_m[m,]         <- rho_m[m,] * c(1,cumprod(1-rho_m[m,1:K_J]))
          
          chi_m[m,c(act_L_J[[m]][act_J],F)] <- ris_VB$a_chi_m[m,act_L_J[[m]]]/ris_VB$a_chi_m[m,act_L_J[[m]]]
          
          Lambda_m[[m]]  <- cbind(ris_VB$m_Loadings_m[[m]][,c(act_J,rep(F,K0_m[m])),drop=F],Lambda_m[[m]])
          
          if(sum(act_G[[m]])>0){
            Gamma_m[[m]] <- cbind(ris_VB$m_Loadings_m[[m]][,c(rep(F,K0),act_G[[m]]),drop=F],Gamma_m[[m]])
            phi_m[[m]]   <- cbind(ris_VB$m_phi_m[[m]][,act_G[[m]],drop=F],phi_m[[m]])
            
            z_m[[m]]   <- c(round(colSums(exp(ris_VB$logP_Zm[[m]])*matrix(1:K0_m[m],K0_m[m],K0_m[m]))[act_G[[m]]]), z_m[[m]])
            
            nu_m[[m]]  <- c(ris_VB$a_nu_m[[m]][act_G[[m]][-K0_m[m]]] / (ris_VB$a_nu_m[[m]][act_G[[m]][-K0_m[m]]] + ris_VB$b_nu_m[[m]][act_G[[m]][-K0_m[m]]]), nu_m[[m]])
            tau_m[[m]] <- c(ris_VB$a_tau_m[[m]][act_G[[m]]] / ris_VB$b_tau_m[[m]][act_G[[m]]], tau_m[[m]])
          }
          
          if(sum(act_L[[m]])>0){
            Gamma_m[[m]] <- cbind(ris_VB$m_Loadings_m[[m]][,c(act_L[[m]],rep(F,K0_m[m])),drop=F],Gamma_m[[m]])
            phi_m[[m]]   <- cbind(ris_VB$m_eta[,act_L[[m]],drop=F],phi_m[[m]])
            
            z_m[[m]]   <- z_m[[m]] + sum(act_L[[m]]) 
            z_m[[m]]   <- c(round(colSums(exp(ris_VB$logP_Dm[m,,])*matrix(1:K0,K0,K0))[act_L[[m]]]), z_m[[m]])
            
            nu_m[[m]]  <- c(ris_VB$a_rho_m[m,act_L[[m]][-K0]] / (ris_VB$a_rho_m[m,act_L[[m]][-K0]] + ris_VB$b_rho_m[m,act_L[[m]][-K0]]), nu_m[[m]])
            tau_m[[m]] <- c(ris_VB$a_chi_m[m,act_L[[m]]] / ris_VB$b_chi_m[m,act_L[[m]]], tau_m[[m]])
          }
          
          if(length(nu_m[[m]])>1){w_m[[m]] <- nu_m[[m]]*c(1,cumprod(1-nu_m[[m]][-length(nu_m[[m]])]))}
          
        }
        
        param_init <- c(param_init,list(Lambda_m=Lambda_m,Gamma_m=Gamma_m,phi_m=phi_m,
                                        delta_m=delta_m,rho_m=rho_m,xi_m=xi_m,chi_m=chi_m,
                                        z_m=z_m,nu_m=nu_m,w_m=w_m,tau_m=tau_m))
        
        return(param_init)
      }
      
      K0    <- K_J+1
      K0_m <- sapply(act_G,sum)+sapply(act_L,sum)+rep(1,M)
    }
  }
  
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

