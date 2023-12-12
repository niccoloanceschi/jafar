
is.scalar <- function(x) is.atomic(x) && length(x) == 1L

colSumsStable <- function(x){if(is.null(dim(drop(x)))){x}else{apply(x,2,prod)}}

get_pred_coeff_JAFAR <- function(M,K,K_m,p_m,Theta,s2_inv_y,Lambda_m,Gamma_m,s2_inv_m,rescale_pred=FALSE){
  beta_m <- list()
  C_inv <- diag(1,K,K)
  
  for(m in 1:M){
    
    # TODO: check behavior ----- #
    mar_std_m = rep(1,p_m[m])
    if(rescale_pred){
      mar_std_m = sqrt(1/s2_inv_m[[m]] + rowSums(Lambda_m[[m]]^2) + rowSums(Gamma_m[[m]]^2))
    }
    Ga_m = Gamma_m[[m]]/mar_std_m
    s2_m = s2_inv_m[[m]]*(mar_std_m^2)
    La_m = Lambda_m[[m]]/mar_std_m
    # TODO: check behavior ----- #
    
    s2_Ga_m = s2_m*Ga_m
    GaT_s2_La_m = crossprod(s2_Ga_m,La_m)
    
    D_m_inv = chol(diag(1.,K_m[m],K_m[m])+crossprod(Ga_m,s2_Ga_m))
    D_GaT_s2_La_m = backsolve(D_m_inv,forwardsolve(t(D_m_inv),GaT_s2_La_m))
    
    C_inv = C_inv + crossprod(La_m,s2_m*La_m) - crossprod(GaT_s2_La_m,D_GaT_s2_La_m)
    
    beta_m[[m]] = t(s2_m*La_m - s2_Ga_m%*%D_GaT_s2_La_m)
  }
  
  C_chol = chol(C_inv)
  Theta_C = backsolve(C_chol,forwardsolve(t(C_chol),Theta))
  
  var_y = 1/s2_inv_y + sum(Theta*Theta_C)
  for(m in 1:M){beta_m[[m]] <- Theta_C%*%beta_m[[m]]}
  
  return(list(pred_coeff_m=beta_m, pred_var=var_y))
}

vb_JAFAR_CUSP <- function(y, X_m, n, M, p_m,                      # input data
                           iterMax=1e3, iter_print=10,            # optimization hyper-params
                           rel_thr=1e-3, seed=123,               # optimization hyper-params
                           K=NULL, K_Gm=NULL,                     # latent dim bounds
                           a_sig=NULL, b_sig=NULL,                # response noise
                           a_theta=NULL, b_theta=NULL,            # slab in response loadings variances
                           var_spike_theta=NULL,                  # spike value in response loadings variances
                           a_xi=NULL, b_xi=NULL,                  # mixture weight in response loadings variances
                           a_m=NULL, b_m=NULL,                    # idiosyncratic noises
                           prec0=NULL, prec0m=NULL,               # "local" means
                           var_spike=NULL, var_spike_loc=NULL,    # spike value in loadings variances
                           a_chi=NULL, b_chi=NULL,                # slab in 'shared' loadings variances
                           a_tau=NULL, b_tau=NULL,                # slab in 'local' loadings variances
                           alpha=NULL, alpha_loc=NULL             # beta dist stick breaking
                          ){
  
  set.seed(seed)
  
  # ------ Hyper-parameters Check ------ #
  
  # bounds on number of latent factors
  if(is.null(K)){K=min(floor(log(p_m)*3))}
  if(is.null(K_Gm)){K_Gm=floor(log(p_m)*3)} else if(is.scalar(K_Gm)){K_Gm=rep(K_Gm,M)}
  
  # hyperparams response noise
  if(is.null(a_sig)){a_sig=1}
  if(is.null(b_sig)){b_sig=0.3}
  
  # hyperparams idiosyncratic noises across modalities
  if(is.null(a_m)){a_m=rep(1,M)} else if(is.scalar(a_m)){a_m=rep(a_m,M)}
  if(is.null(b_m)){b_m=rep(0.3,M)} else if(is.scalar(b_m)){b_m=rep(b_m,M)}
  
  # hyperparams "local" mean for responses & across modalities
  if(is.null(prec0)){prec0=1/0.5}
  if(is.null(prec0m)){prec0m=rep(1/0.5,M)} else if(is.scalar(prec0m)){prec0m=rep(prec0m,M)}
  
  # hyperparams slab in response loadings variances
  if(is.null(a_theta)){a_theta=3}
  if(is.null(b_theta)){b_theta=0.1}
  
  # hyperparams spike in response loadings variances
  if(is.null(var_spike_theta)){var_spike_theta=0.005}
  
  # hyperparams mixutre weight in response loadings variances
  if(is.null(a_xi)){a_xi=5}
  if(is.null(b_xi)){b_xi=5}
  
  # hyperparam spike value in loadings variances
  if(is.null(var_spike)){var_spike=rep(0.005,M)} else if(is.scalar(var_spike)){var_spike=rep(var_spike,M)}
  if(is.null(var_spike_loc)){var_spike_loc=rep(0.005,M)} else if(is.scalar(var_spike_loc)){var_spike_loc=rep(var_spike_loc,M)}
  
  # hyperparams slab in 'shared' loadings variances
  if(is.null(a_chi)){a_chi=rep(1,M)} else if(is.scalar(a_chi)){a_chi=rep(a_chi,M)}
  if(is.null(b_chi)){b_chi=rep(0.2,M)} else if(is.scalar(b_chi)){b_chi=rep(b_chi,M)}
  
  # hyperparams slab in 'local' loadings variances
  if(is.null(a_tau)){a_tau=rep(1,M)} else if(is.scalar(a_tau)){a_tau=rep(a_tau,M)}
  if(is.null(b_tau)){b_tau=rep(0.2,M)} else if(is.scalar(b_tau)){b_tau=rep(b_tau,M)}
  
  # hyperparam beta dist stick breaking
  if(is.null(alpha)){alpha=rep(5,M)} else if(is.scalar(alpha)){alpha=rep(alpha,M)}
  if(is.null(alpha_loc)){alpha_loc=rep(5,M)} else if(is.scalar(alpha_loc)){alpha_loc=rep(alpha_loc,M)}
  
  # ------ Working Variables Initialization ------ #
  
  m_Theta   <- rep(NA,K)
  v_Theta   <- matrix(NA,K,K)
  
  m_eta     <- matrix(rnorm(n*K),n,K)
  v_eta     <- diag(1,K,K)
  
  a_s2_inv  <- a_sig+0.5*n   
  b_s2_inv  <- b_sig+0.5*n
  
  m_mu_y    <- 0.
  v_mu_y    <- 0.05
  
  a_chi_vb  <- a_theta + 0.5
  b_chi_vb  <- b_theta + rep(0.5*var_spike_theta,K)
  
  a_xi_vb   <- a_xi+0.5*K
  b_xi_vb   <- a_xi+0.5*K
  
  # logP_D[m,l,h] = log(q[\delta = l])
  logP_D      <- matrix(log(0.5),2,K)                     
  
  # logP_Dm_leq[m,l,h] = log(q[\delta_{m h} \leq l]) = log(\sum_{s \leq l}q[\delta_{m h} = s])
  logP_Dm <- array(0,c(M,K,K))
  logP_Dm_leq <- log(outer(rep(1,M),matrix(cumsum(rep(1/K,K)),K,K),"*"))
  
  # NB: the last column is always inactive (rho_{m K}=1 a.s.)
  a_rho_m <- matrix(1+rep(1,K-1),M,K-1,byrow=T)       
  b_rho_m <- matrix(alpha,M,K-1) + matrix(seq(K-1,1,by=-1),M,K-1,byrow=T)
  
  a_chi_m   <- matrix(a_chi+0.5*p_m,M,K)
  b_chi_m   <- matrix(b_chi+0.5*p_m*b_chi/(a_chi-1),M,K)
  
  m_Loadings_m <- v_Loadings_m <- m_phi_m <- v_phi_m <- list()
  a_s2_inv_m <- b_s2_inv_m <- m_mu_m <- v_mu_m <- list()
  logP_Zm <- logP_Zm_leq <- a_nu_m <- b_nu_m <- a_tau_m <- b_tau_m <- list()
  
  for(m in c(1:M)){
    
    m_phi_m[[m]] <- matrix(rnorm(n*K_Gm[m]),n,K_Gm[m])
    v_phi_m[[m]] <- diag(1,K_Gm[m],K_Gm[m])
   
    m_Loadings_m[[m]] <- matrix(NA,p_m[m],K+K_Gm[m])
    v_Loadings_m[[m]] <- array(NA,c(p_m[m],K+K_Gm[m],K+K_Gm[m]))
    
    a_s2_inv_m[[m]] <- rep(a_m[m]+0.5*n,p_m[m])
    b_s2_inv_m[[m]] <- rep(b_m[m]+0.5*n,p_m[m])
      
    m_mu_m[[m]] <- rep(0,p_m[m])
    v_mu_m[[m]] <- rep(0.05,p_m[m])
    
    # logP_Zm_leq[[m]][l,h] = log(q[z_{m h} \leq l]) = = log(\sum_{s \leq l}q[z_{m h} = s])
    logP_Zm_leq[[m]] <- matrix(log(cumsum(rep(1/K_Gm[m],K_Gm[m]))),K_Gm[m],K_Gm[m]) 
    
    # NB: the last column is always inactive (nu_{m K_m}=1 a.s.)
    a_nu_m[[m]] <- 1+rep(1,K_Gm[m]-1)       
    b_nu_m[[m]] <- alpha_loc[m]+seq(K_Gm[m]-1,1,by=-1) 
    
    a_tau_m[[m]]   <- rep(a_tau[m]+0.5*p_m[m],M,K_Gm[m])
    b_tau_m[[m]]   <- rep(b_tau[m]+0.5*p_m[m]*b_tau[m]/(a_tau[m]-1),M,K_Gm[m])                   
  }
  
  # ------ Auxiliary Variables Initialization ------ #
  
  # Auxiliary Variables for delta & delta_m updates
  
  # logP_Dm_leq_Diag[m,h] = log(q[\delta_{m,h} \leq h])
  logP_Dm_leq_Diag <- t(apply(logP_Dm_leq,1,diag))
  
  # logP_Dm_leq_Sum[h] = log(\sum_{m=1}^M q[\delta_{m,h} \leq h])
  logP_Dm_leq_Sum  <- colSums(logP_Dm_leq_Diag)
    
  # Auxiliary Variables for Loadings & Idiosyncratic Components
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
  
  # Log-Determinants for ELBO computation
  
  lD_Theta <- 0
  lD_eta <- 0
  lD_Loadings_m <- lapply(1:M, function(m) rep(0,p_m[m]))
  lD_phi_m <- lapply(1:M, function(m) 0)
  vb_elbo <- rep(NA,iterMax)
  
  # ------ VB Updates ----------------------------------------------------------
  
  t=0
  elbo_old = -Inf
  elbo_check = F
  
  # ------ print status ------- #
  print(paste0("VB Iteration n. ",t))
  
  while(!elbo_check & t<iterMax){
    
    t <- t+1
        
    # STEP 1: Loadings ---------------------------------------------------------
    
    ## 1.a Theta ---------------------------------------------------------------
    
    prob_chi <- exp(logP_D[2,])*(1-exp(logP_Dm_leq_Sum))
    
    E_chi <- prob_chi*(a_chi_vb/b_chi_vb) + (1-prob_chi)*rep(1/var_spike_theta,K)
    
    Q_Theta <- diag(E_chi,K,K) + (a_s2_inv/b_s2_inv)*(n*v_eta+etaTeta)
    r_Theta <- (a_s2_inv/b_s2_inv)*crossprod(m_eta,y-rep(m_mu_y,n))
    
    L_Theta <- chol(Q_Theta)
    lD_Theta <- 2 * sum(log(diag(L_Theta)))
    
    v_Theta <- chol2inv(L_Theta)
    m_Theta <- as.vector(backsolve(L_Theta, forwardsolve(t(L_Theta), r_Theta)))
    
    ## 1.b Lambda_m, Gamma_m ---------------------------------------------------
    
    for(m in c(1:M)){

      prob_chi_m <- (1-exp(logP_Dm_leq_Diag[m,]))*(1-exp(logP_D[1,]+colSums(logP_Dm_leq_Diag[-m,])))
      prob_tau_m <- (1-exp(diag(logP_Zm_leq[[m]])))

      E_chi_m    <- prob_chi_m*a_chi_m[m,]/b_chi_m[m,] + (1-prob_chi_m)*rep(1/var_spike[m],K)
      E_tau_m    <- prob_tau_m*a_tau_m[[m]]/b_tau_m[[m]] + (1-prob_tau_m)*rep(1/var_spike_loc[m],K_Gm[m])
      
      diag_prec_m <- diag(c(E_chi_m,E_tau_m),K+K_Gm[m],K+K_Gm[m])
      
      for(j in c(1:p_m[m])){
        
        Q_Loadings_mj <- diag_prec_m + (a_s2_inv_m[[m]][j]/b_s2_inv_m[[m]][j])*facTfac_m[[m]]
        r_Loadings_mj <- (a_s2_inv_m[[m]][j]/b_s2_inv_m[[m]][j])*crossprod(cbind(m_eta,m_phi_m[[m]]),X_m[[m]][,j]-rep(m_mu_m[[m]][j],n))
        
        L_Loadings_mj <- chol(Q_Loadings_mj)
        lD_Loadings_m[[m]][j] <- 2 * sum(log(diag(L_Loadings_mj)))
        
        v_Loadings_m[[m]][j,,] <- chol2inv(L_Loadings_mj)
        m_Loadings_m[[m]][j,]  <- backsolve(L_Loadings_mj, forwardsolve(t(L_Loadings_mj), r_Loadings_mj))
        
      }
    }
    
    # STEP 2: Intercepts -------------------------------------------------------
    
    ## 2.a mu_y ----------------------------------------------------------------
    
    if(t>10){
      v_mu_y <- 1/(n*(a_s2_inv/b_s2_inv)+prec0)
      m_mu_y <- (a_s2_inv/b_s2_inv)*v_mu_y*sum(y-m_eta%*%m_Theta)
    }
    
    ## 2.b mu_m ----------------------------------------------------------------
    
    if(t>10){
      for(m in 1:M){
        v_mu_m[[m]] <- 1/(n*(a_s2_inv_m[[m]]/b_s2_inv_m[[m]])+rep(prec0m[m],p_m[m]))
        m_mu_m[[m]] <- (a_s2_inv_m[[m]]/b_s2_inv_m[[m]])*v_mu_m[[m]]*
          colSums(X_m[[m]]-tcrossprod(cbind(m_eta,m_phi_m[[m]]),m_Loadings_m[[m]]))
      }
    }
    
    # STEP 3: Idiosyncratic components -----------------------------------------
    
    ## 3.a s2_inv --------------------------------------------------------------
    
    # a_s2_inv  <- a_sig+0.5*n   
    b_s2_inv  <- as.numeric( b_sig + 0.5*sum(y^2) - m_mu_y*sum(y) + 0.5*n*(v_mu_y+m_mu_y^2) -
        sum((y-rep(m_mu_y,n))*(m_eta%*%m_Theta)) + 0.5*sum(diag((n*v_eta+etaTeta)%*%(v_Theta+tcrossprod(m_Theta)))) )
    
    ## 3.b s2_inv_m ------------------------------------------------------------
    
    for(m in c(1:M)){
      
      # a_s2_inv_m[[m]] <- rep(a_m[m]+0.5*n,p_m[m])
      for(j in 1:p_m[m]){
        Load_Load_T <- v_Loadings_m[[m]][j,,] + tcrossprod(m_Loadings_m[[m]][j,])
        
        b_s2_inv_m[[m]][j] <- as.numeric( b_m[m] + 0.5*sum(X_m[[m]][,j]^2) -
            m_mu_m[[m]][j]*sum(X_m[[m]][,j]) + 0.5*n*(v_mu_m[[m]][j]+m_mu_m[[m]][j]^2) - 
            (X_m[[m]][,j]-rep(m_mu_m[[m]][j],n))%*%(cbind(m_eta,m_phi_m[[m]])%*%m_Loadings_m[[m]][j,]) +
             0.5*sum(diag(facTfac_m[[m]]%*%Load_Load_T)) )
      }
      
    }
    
    # STEP 4: Factors ----------------------------------------------------------
    
    res_m <- list()
    for(m in 1:M){
      res_m[[m]] <- X_m[[m]]-matrix(m_mu_m[[m]],n,p_m[m],byrow=T)
    }
    
    ## 4.a eta -----------------------------------------------------------------
    
    Q_eta <- diag(1,K,K) + (a_s2_inv/b_s2_inv)*(v_Theta+tcrossprod(m_Theta))
    r_eta <- (a_s2_inv/b_s2_inv)*tcrossprod(m_Theta,y-rep(m_mu_y,n))
    
    for(m in 1:M){
      s_L_m <- (a_s2_inv_m[[m]]/b_s2_inv_m[[m]])*m_Loadings_m[[m]][,1:K]
      
      Q_eta <- Q_eta + crossprod(m_Loadings_m[[m]][,1:K],s_L_m) +
        colSums((a_s2_inv_m[[m]]/b_s2_inv_m[[m]])*v_Loadings_m[[m]][,1:K,1:K]) 
      r_eta <- r_eta + t(res_m[[m]]%*%s_L_m) - tcrossprod(crossprod(s_L_m,m_Loadings_m[[m]][,-c(1:K)]) +
        colSums((a_s2_inv_m[[m]]/b_s2_inv_m[[m]])*v_Loadings_m[[m]][,1:K,-c(1:K)]), m_phi_m[[m]])
    }
    
    L_eta <- chol(Q_eta)
    lD_eta <- 2 * sum(log(diag(L_eta)))
    
    v_eta <- chol2inv(L_eta)
    m_eta <- t(backsolve(L_eta, forwardsolve(t(L_eta), r_eta)))
    
    ## 4.b phi_m ---------------------------------------------------------------
    for(m in c(1:M)){
      
      s_G_m <- (a_s2_inv_m[[m]]/b_s2_inv_m[[m]])*m_Loadings_m[[m]][,-c(1:K)]
      
      Q_phi_m <- diag(1,K_Gm[m],K_Gm[m]) +
        colSums((a_s2_inv_m[[m]]/b_s2_inv_m[[m]])*v_Loadings_m[[m]][,-c(1:K),-c(1:K)]) +
        crossprod(m_Loadings_m[[m]][,-c(1:K)],s_G_m)
      r_phi_m <- t(res_m[[m]]%*%s_G_m) - tcrossprod(crossprod(s_G_m,m_Loadings_m[[m]][,1:K]) +
        colSums((a_s2_inv_m[[m]]/b_s2_inv_m[[m]])*v_Loadings_m[[m]][,-c(1:K),1:K]), m_eta)
      
      L_phi_m <- chol(Q_phi_m)
      lD_phi_m[[m]] <- 2 * sum(log(diag(Q_phi_m)))
      
      v_phi_m[[m]] <- chol2inv(L_phi_m)
      m_phi_m[[m]] <- t(backsolve(L_phi_m, forwardsolve(t(L_phi_m), r_phi_m)))
    }
        
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
    
    # STEP 5: Latent indicatirs ------------------------------------------------
    
    ## 5.a Shared Component ----------------------------------------------------
        
    # efficient computation of multivariate Normal & Student-T (log)-pdf
    
    E_T_2 <- diag(v_Theta) + m_Theta^2
    V_T_2 <- 2*diag(v_Theta)^2 + 4*m_Theta^2*diag(v_Theta)
    
    logP_Spikes <- - 0.5*E_T_2/var_spike_theta - rep(0.5*log(2*pi*var_spike_theta),K)
    logP_Slabs  <- - (0.5+a_theta)*log(1+0.5*E_T_2^2/b_theta) +
      0.5*(0.5+a_theta)*(0.25*V_T_2^2/b_theta^2)/(1+0.5*E_T_2^2/b_theta)^2 +
      rep(lgamma(0.5+a_theta) - lgamma(a_theta) - 0.5*log(2*pi*b_theta),K)
    
    logP_T_Diff <- logP_Slabs - logP_Spikes
    
    logP_L_Diff <- matrix(0,nrow=M,K)
    for(m in 1:M){
      
      E_Lm_2 <- rowSums(apply(v_Loadings_m[[m]][,1:K,1:K],1,diag)) + colSums(m_Loadings_m[[m]][,1:K]^2)
      V_Lm_2 <- 2*rowSums(apply(v_Loadings_m[[m]][,1:K,1:K],1,diag)^2) +
        4*rowSums(apply(v_Loadings_m[[m]][,1:K,1:K],1,diag)*t(m_Loadings_m[[m]][,1:K]^2))
      
      logP_Spikes <- - 0.5*E_Lm_2/var_spike[m] - rep(0.5*p_m[m]*log(2*pi*var_spike[m]),K) 
      logP_Slabs  <- - (0.5*p_m[m]+a_chi[m])*log(1+0.5*E_Lm_2/b_chi[m]) +
        0.5*(0.5*p_m[m]+a_chi[m])*(0.25*V_Lm_2/b_chi[m]^2)/(1+0.5*E_Lm_2/b_chi[m])^2 + 
        rep(lgamma(0.5*p_m[m]+a_chi[m]) - lgamma(a_chi[m]) - 0.5*p_m[m]*log(2*pi*b_chi[m]),K)
      
      logP_L_Diff[m,] <- logP_Slabs - logP_Spikes
    }
    
    ### 5.a.1 delta ------------------------------------------------------------
    
    # Remarks:  (i)   a priori: $ \xi = P[delta_h = 0] )
    #           (ii)  logP_D[l,h] = log(q[\delta_h = l])
    
    logP_D <- matrix(c(digamma(a_xi_vb),digamma(b_xi_vb)),2,K)
    logP_D <- logP_D + c(0,1)%x%t(logP_T_Diff*(1-exp(logP_Dm_leq_Sum)))
    for(m in 1:M){
      logP_D <- logP_D + matrix(logP_L_Diff[m,]*(1-exp(logP_Dm_leq_Diag[m,])),2,K,byrow=T)*
        (1-c(1,0)%x%t(exp(apply(logP_Dm_leq_Diag[-m,],2,sum))))
    }
    
    # ----- "-Max" in Log-Sum-Exp Trick ----- #
    logP_D <- logP_D - matrix(apply(logP_D,2,max),2,K,byrow=T)     
    
    # ----- # Normalized Log-Probabilities ----- #
    logP_D <- logP_D - matrix(log(colSums(exp(logP_D))),2,K,byrow=T) 
    
    ### 5.a.2 delta_m ----------------------------------------------------------
    
    # Remark:  logP_Dm[m,l,h] = log(q[\delta_{m h} = l])
    
    I_lh_leq <- ( matrix(1:K,K,K) <= matrix(1:K,K,K,byrow=T) ) # Auxiliary Indicators
    
    for(m in 1:M){
      
      a_digamma <- c(digamma(a_rho_m[m,]) - digamma(a_rho_m[m,]+b_rho_m[m,]),0)
      b_digamma <- c(0,cumsum(digamma(b_rho_m[m,]) - digamma(a_rho_m[m,]+b_rho_m[m,])))
      
      logP_Dm_new <- a_digamma%x%t(rep(1,K)) + b_digamma%x%t(rep(1,K))
      logP_Dm_new <- logP_Dm_new + matrix(logP_T_Diff*exp(logP_D[2,]),K,K,byrow=T)*
        (1-I_lh_leq*matrix(exp(apply(logP_Dm_leq_Diag[-m,],2,sum)),K,K,byrow=T))
      logP_Dm_new <- logP_Dm_new + matrix(logP_L_Diff[m,],K,K,byrow=T)*(1-I_lh_leq)*
        (1-matrix(exp(logP_D[1,]+apply(logP_Dm_leq_Diag[-m,],2,sum)),K,K,byrow=T))
      for(mm in c(1:M)[-m]){
        logP_Dm_new <- logP_Dm_new + matrix(logP_L_Diff[mm,]*(1-exp(logP_Dm_leq_Diag[mm,])),K,K,byrow=T)*
          (1-I_lh_leq*matrix(exp(logP_D[1,]+colSumsStable(logP_Dm_leq_Diag[-c(m,mm),])),K,K,byrow=T))
      }
      
      # ----- "-Max" in Log-Sum-Exp Trick ----- #
      logP_Dm_new <- logP_Dm_new - matrix(apply(logP_Dm_new, 2, max),K,K,byrow=T)     
    
      # ----- Normalized Probabilities ----- #
      logP_Dm[m,,] <- logP_Dm_new - matrix(log(colSums(exp(logP_Dm_new))),K,K,byrow=T)
      
      # ----- Cumulative Probabilities ----- #
      # logP_Dm_leq[m,l,h] = log(q[\delta_{m,h} \leq l]) = log(\sum_{s \leq l} q[\delta_{m,h} = s])
      
      logP_Dm_cummax <- apply(logP_Dm_new,2,cummax)
      cummax_tensor  <- aperm(outer(logP_Dm_cummax,rep(1,K),"*"),c(1,3,2))
      logP_Dm_tensor <- outer(rep(1,K),logP_Dm_new,"*")
      
      logP_Dm_leq[m,,] <- logP_Dm_cummax +
        log(apply(apply(exp(logP_Dm_tensor-cummax_tensor),c(1,3),cumsum),3,diag)) -
        log(rep(1,K)%x%t(colSums(exp(logP_Dm_new))))
      
      # ----- Auxiliary probabilities ----- #
      logP_Dm_leq_Diag <- t(apply(logP_Dm_leq,1,diag))   # logP_Dm_leq_Diag[m,h] = log(q[\delta_{m,h} \leq h])
      logP_Dm_leq_Sum  <- colSums(logP_Dm_leq_Diag)      # logP_Dm_leq_Sum[h] = log(\sum_{m=1}^M q[\delta_{m,h} \leq h])
    }
    
    ## 5.b Specific Components -------------------------------------------------
    
    for(m in 1:M){

      # efficient computation of multivariate Normal (log)-pdf
      
      E_Gm_2 <- rowSums(apply(v_Loadings_m[[m]][,-c(1:K),-c(1:K)],1,diag)) + colSums(m_Loadings_m[[m]][,-c(1:K)]^2)
      V_Lm_2 <- 2*rowSums(apply(v_Loadings_m[[m]][,-c(1:K),-c(1:K)],1,diag)^2) +
        4*rowSums(apply(v_Loadings_m[[m]][,-c(1:K),-c(1:K)],1,diag)*t(m_Loadings_m[[m]][,-c(1:K)]^2))
      
      logP_G_Spikes <- - 0.5*E_Gm_2/var_spike_loc[m] - rep(0.5*p_m[m]*log(2*pi*var_spike_loc[m]),K_Gm[m]) 
      logP_G_Slabs  <- - (0.5*p_m[m]+a_tau[m])*log(1+0.5*E_Gm_2/b_tau[m]) +
        0.5*(0.5*p_m[m]+a_tau[m])*(0.25*V_Lm_2/b_tau[m]^2)/(1+0.5*E_Gm_2/b_tau[m])^2 + 
        rep(lgamma(0.5*p_m[m]+a_tau[m]) - lgamma(a_tau[m]) - 0.5*p_m[m]*log(2*pi*b_tau[m]),K_Gm[m])
      
      logP_G_Diff <- logP_G_Slabs - logP_G_Spikes
      
      ### 5.b.1 z_m (latent indicators) ----
      # Remark:  logP_Zm[[m]][l,h] = log(q[z_{m h} = l])
      
      a_digamma <- c(digamma(a_nu_m[[m]]) - digamma(a_nu_m[[m]]+b_nu_m[[m]]),0)
      b_digamma <- c(0,cumsum(digamma(b_nu_m[[m]]) - digamma(a_nu_m[[m]]+b_nu_m[[m]])))
      
      I_lh_gr <- matrix(1:K_Gm[m],K_Gm[m],K_Gm[m]) > matrix(1:K_Gm[m],K_Gm[m],K_Gm[m],byrow=T) 
      
      logP_Zm_new <- a_digamma%x%t(rep(1,K_Gm[m])) + b_digamma%x%t(rep(1,K_Gm[m])) +
        (rep(1,K_Gm[m])%x%t(logP_G_Diff))*I_lh_gr
      
      # ----- "-Max" in Log-Sum-Exp Trick ----- #
      logP_Zm_new <- logP_Zm_new - matrix(apply(logP_Zm_new, 2, max),K_Gm[m],K_Gm[m],byrow=T)     
      
      # ----- Normalized Probabilities ----- #
      logP_Zm[[m]] <- logP_Zm_new - matrix(log(colSums(exp(logP_Zm_new))),K_Gm[m],K_Gm[m],byrow=T)
      
      # ----- Cumulative Probabilities ----- #
      # logP_Zm_leq[[m]][l,h] = log(q[z_{m,h} \leq l]) = log(q[z_{m,h} \leq l]) = log(\sum_{s \leq l} q[z_{m,h} = s])
      
      logP_Zm_cummax <- apply(logP_Zm_new,2,cummax)
      cummax_tensor  <- aperm(outer(logP_Zm_cummax,rep(1,K_Gm[m]),"*"),c(1,3,2))
      logP_Zm_tensor <- outer(rep(1,K_Gm[m]),logP_Zm_new,"*")
      
      logP_Zm_leq[[m]] <- logP_Zm_cummax +
        log(apply(apply(exp(logP_Zm_tensor-cummax_tensor),c(1,3),cumsum),3,diag)) -
        log(rep(1,K_Gm[m])%x%t(colSums(exp(logP_Zm_new))))
      
    }
    # ----- Auxiliary probabilities ----- #
    logP_Zm_leq_Diag <- lapply(logP_Zm_leq,diag)
          
    # STEP 6: Stick Breaking Elements ------------------------------------------
    
    ## 6.a Theta ---------------------------------------------------------------
    
    ### 6.a.2 xi (spike and slab mixture weight) ----
    # Remark: a priori $ \xi = P[delta_h = 0] $
    
    a_xi_vb    <- a_xi+sum(exp(logP_D[1,]))
    b_xi_vb    <- b_xi+sum(exp(logP_D[2,]))
    
    ### 6.a.2 chi (response loadings precisions) ----
    
    # a_chi_vb <- a_theta + 0.5
    b_chi_vb   <- as.vector(b_theta + 0.5*(diag(v_Theta)+m_Theta^2))
    
    ## 6.b Lambda_m ------------------------------------------------------------
    
    ### 6.b.1 rho_m (stick breaking elements) ----
    for(m in 1:M){
      
    # a_rho_m <- matrix(1,M,K-1) + t(apply(aperm(apply(exp(logP_Dm_leq), c(1,2), cumsum),c(2,3,1)),1,diag))[,-K] -
    #   cbind(rep(0,M),t(apply(aperm(apply(exp(logP_Dm_leq), c(1,2), cumsum),c(2,3,1)),1,diag))[,-c(K-1,K)])
    # b_rho_m <- alpha%x%rep(1,K) + t(apply(aperm(apply(1-exp(logP_Dm_leq), c(1,2), cumsum),c(2,3,1)),1,diag))[,-K]
    
      a_rho_m[m,] <- 1 + diag(t(apply(exp(logP_Dm_leq[m,,]), 1, cumsum)))[-K] - c(0,diag(t(apply(exp(logP_Dm_leq[m,,]), 1, cumsum)))[-c(K-1,K)])
      b_rho_m[m,] <- alpha[m] + diag(t(apply(1-exp(logP_Dm_leq[m,,]), 1, cumsum)))[-K]
    }
      
    ### 6.b.2 chi_m (loading precisions) ----
    for(m in 1:M){
      # a_chi_m[m,] <- a_chi[m] + 0.5*p_m[m] 
      b_chi_m[m,] <- b_chi[m] + 0.5 * (rowSums(apply(v_Loadings_m[[m]][,1:K,1:K], 1, diag)) + colSums(m_Loadings_m[[m]][,1:K]^2) )
    }

    ## 6.c Gamma_m -------------------------------------------------------------
    
    ### 6.c.1 nu_m (stick breaking elements) ----
    for(m in 1:M){
      a_nu_m[[m]] <- 1 + diag(t(apply(exp(logP_Zm_leq[[m]]), 1, cumsum)))[-K_Gm[m]] - c(0,diag(t(apply(exp(logP_Zm_leq[[m]]), 1, cumsum)))[-c(K_Gm[m]-1,K_Gm[m])])
      b_nu_m[[m]] <- alpha_loc[m] + diag(t(apply(1-exp(logP_Zm_leq[[m]]), 1, cumsum)))[-K_Gm[m]]
    }
      
    ### 6.c.2 tau_m (loading precisions) ----  
    for(m in 1:M){
      # a_tau_m[[m]] <- a_tau[m] + 0.5*p_m[m] 
      b_tau_m[[m]] <- b_tau[m] + 0.5 * (rowSums(apply(v_Loadings_m[[m]][,-c(1:K),-c(1:K)], 1, diag)) + colSums(m_Loadings_m[[m]][,-c(1:K)]^2) )
    }
    
    # ------ print status ------- #
    if(t %% iter_print == 0){
      print(paste0("VB Iteration n. ",t))
    }
    
    vb_elbo[t] <- JAFAR_elbo_vb(t,y,X_m,n,M,p_m,K,K_Gm,
                                which_cusp='orig',
                                a_sig, b_sig,                
                                var_spike_theta,             
                                a_xi, b_xi,                  
                                a_m, b_m,                    
                                prec0, prec0m,               
                                var_spike, var_spike_loc,    
                                alpha, alpha_loc,            
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
                                a_chi_vb=a_chi_vb, b_chi_vb=b_chi_vb,
                                a_chi_m=a_chi_m, b_chi_m=b_chi_m,
                                a_tau_m=a_tau_m, b_tau_m=b_tau_m)
    
    elbo_check <- abs((vb_elbo[t]-elbo_old)/vb_elbo[t]) < rel_thr
    elbo_old   <- vb_elbo[t]
    
  }
  
  # ------ Saving Active Columns Probabilities ------- #
  
  P_J_active <- exp(logP_Dm_leq_Sum)
  for(m in 1:M){
    P_J_active <- P_J_active + (1-exp(logP_Dm_leq_Diag[m,]))*exp(logP_D[1,]+colSums(logP_Dm_leq_Diag[-m,]))
  }
  P_J_active <- 1 - P_J_active
  
  # ------ Predictions ------- #
  
  m_s2_inv_y <- (a_s2_inv/b_s2_inv)
  m_s2_inv_m <- lapply(1:M, function(m) (a_s2_inv_m[[m]]/b_s2_inv_m[[m]]))
  m_Lambda_m <- lapply(1:M, function(m) m_Loadings_m[[m]][,1:K])
  m_Gamma_m  <- lapply(1:M, function(m) m_Loadings_m[[m]][,-c(1:K)])
                                      
  pred_coeff <- get_pred_coeff_JAFAR(M,K,K_Gm,p_m,m_Theta,m_s2_inv_y,m_Lambda_m,m_Gamma_m,m_s2_inv_m,rescale_pred=F)
  pred_coeff_rescaled <- get_pred_coeff_JAFAR(M,K,K_Gm,p_m,m_Theta,m_s2_inv_y,m_Lambda_m,m_Gamma_m,m_s2_inv_m,rescale_pred=T)
    
  # ------ Output ------- #
  
  hyper_params = list(model='jafar_joint_cusp',
                      K=K, K_Gm=K_Gm, a_theta=a_theta, b_theta=b_theta, 
                      a_xi=a_xi, b_xi=b_xi, var_spike_theta=var_spike_theta, 
                      a_sig=a_sig, b_sig=b_sig, prec0=prec0, 
                      a_m=a_m, b_m=b_m, prec0m=prec0m,
                      a_chi=a_chi, b_chi=b_chi, var_spike=var_spike, alpha=alpha,
                      a_tau=a_tau, b_tau=b_tau, var_spike_loc=var_spike_loc, alpha_loc=alpha_loc)
  
  output = list(hyper_param=hyper_params,vb_elbo=vb_elbo[1:t],P_J_active=P_J_active,
                pred_coeff=pred_coeff,pred_coeff_rescaled=pred_coeff_rescaled,
                m_Theta=m_Theta,v_Theta=v_Theta,m_eta=m_eta,v_eta=v_eta,
                a_s2_inv=a_s2_inv,b_s2_inv=b_s2_inv,m_mu_y=m_mu_y,v_mu_y=v_mu_y,
                a_chi_vb=a_chi_vb,b_chi_vb=b_chi_vb,a_xi_vb=a_xi_vb,b_xi_vb=b_xi_vb,
                logP_D=logP_D,logP_Dm=logP_Dm,logP_Zm=logP_Zm,
                logP_Dm_leq=logP_Dm_leq,logP_Zm_leq=logP_Zm_leq,
                a_rho_m=a_rho_m,b_rho_m=b_rho_m,a_chi_m=a_chi_m,b_chi_m=b_chi_m,
                m_phi_m=m_phi_m,v_phi_m=v_phi_m,m_mu_m=m_mu_m,v_mu_m=v_mu_m, 
                m_Loadings_m=m_Loadings_m,v_Loadings_m=v_Loadings_m,
                a_s2_inv_m=a_s2_inv_m,b_s2_inv_m=b_s2_inv_m,
                a_nu_m=a_nu_m,b_nu_m=b_nu_m,a_tau_m=a_tau_m,b_tau_m=b_tau_m)
  
  return(output)
}











