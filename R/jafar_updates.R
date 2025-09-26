
update_loadings_R <- function(n, p_m, Ktot, X_m, facTfac, fac_m, mu_m, s2_inv_m, prec_m){
  
  new_Loadings_m <- matrix(0,p_m,Ktot)
  
  prior_prec_mj <- diag(prec_m,Ktot,Ktot)
  
  for(j in c(1:p_m)){
    
    Q_Loadings_mj <- prior_prec_mj+s2_inv_m[j]*facTfac
    r_Loadings_mj <- s2_inv_m[j]*crossprod(fac_m,X_m[,j]-rep(mu_m[j],n))
    
    L_Loadings_mj  <- chol(Q_Loadings_mj)
    Lr_Loadings_mj <- forwardsolve(t(L_Loadings_mj), r_Loadings_mj)
    
    mean_Loadings_mj <- backsolve(L_Loadings_mj, Lr_Loadings_mj)
    std_Loadings_mj  <- backsolve(L_Loadings_mj, rnorm(K+K_Gm))
    
    new_Loadings_m[j,] <- t(mean_Loadings_mj + std_Loadings_mj)
  }
  
  return(new_Loadings_m)
}

update_factors_collapsed <- function(n, M, K, K_Gm, res_m, s2_inv_m, Lambda_m, Gamma_m){
  
  phi_m_new <- s2_La_m <- s2_Ga_m <- D_m_chol <- list()
  
  # Shared quantities used for sampling both eta phi_m
  for(m in 1:M){
    s2_La_m[[m]]  <- s2_inv_m[[m]]*Lambda_m[[m]]
    s2_Ga_m[[m]]  <- s2_inv_m[[m]]*Gamma_m[[m]]
    D_m_chol[[m]] <- t(chol(diag(1.,K_Gm[m],K_Gm[m])+crossprod(Gamma_m[[m]],s2_Ga_m[[m]])))
  }
  
  ## 4.a eta ----------------------------------------------------------------- #
  
  Q_eta <- diag(1,K,K) 
  r_eta <- rep(0,K)
  
  for(m in 1:M){
    
    Ga_T_s2_La_m     <- crossprod(s2_Ga_m[[m]],Lambda_m[[m]])
    Dinv_GaT_s2_La_m <- backsolve(t(D_m_chol[[m]]),forwardsolve(D_m_chol[[m]], Ga_T_s2_La_m))
    
    Q_eta <- Q_eta + crossprod(Lambda_m[[m]],s2_La_m[[m]]) -
      crossprod(Ga_T_s2_La_m,Dinv_GaT_s2_La_m)
    r_eta <- r_eta + t(res_m[[m]]%*%s2_La_m[[m]]) -
      t((res_m[[m]]%*%s2_Ga_m[[m]])%*%Dinv_GaT_s2_La_m)
  }
  
  L_eta    <- chol(Q_eta)
  Lr_eta   <- forwardsolve(t(L_eta), r_eta)

  mean_eta <- backsolve(L_eta, Lr_eta)
  std_eta  <- backsolve(L_eta, matrix(rnorm(K*n),K,n))
  
  eta_new  <- t(mean_eta + std_eta)
  
  ## 4.b phi_m --------------------------------------------------------------- #
  for(m in c(1:M)){
    
    r_phi_m     <- t((res_m[[m]]-tcrossprod(eta_new,Lambda_m[[m]]))%*%s2_Ga_m[[m]])
    Lr_phi_m    <- forwardsolve(D_m_chol[[m]], r_phi_m)
    
    mean_phi_m  <- backsolve(t(D_m_chol[[m]]), Lr_phi_m)
    std_phi_m   <- backsolve(t(D_m_chol[[m]]), matrix(rnorm(K_Gm[m]*n),K_Gm[m],n))
    
    phi_m_new[[m]]  <- t(mean_phi_m + std_phi_m)
  }
  
  return(list(eta=eta_new, phi_m=phi_m_new))
}

update_factors_supervised <- function(n, M, K, K_Gm, res_m, s2_inv_m, Lambda_m, Gamma_m,
                                      res_y, s2_inv, Theta_all){
  
  
  Q_fact <- s2_inv*tcrossprod(Theta_all) + diag(1,K+sum(K_Gm),K+sum(K_Gm))
  r_fact <- s2_inv*tcrossprod(Theta_all,res_y)
  
  for(m in 1:M){
    idx_m <- K + c(0,cumsum(K_Gm)[-M])[m] + 1:K_Gm[m]
    
    s2_La_m  <- s2_inv_m[[m]]*Lambda_m[[m]]
    s2_Ga_m  <- s2_inv_m[[m]]*Gamma_m[[m]]
    
    Q_fact[1:K,1:K] <- Q_fact[1:K,1:K] + crossprod(Lambda_m[[m]],s2_La_m)
    r_fact[1:K,]    <- r_fact[1:K,] + t(res_m[[m]]%*%s2_La_m)
    
    Q_fact[idx_m,idx_m] <- Q_fact[idx_m,idx_m] + crossprod(Gamma_m[[m]],s2_Ga_m)
    r_fact[idx_m,]      <- r_fact[idx_m,] + t(res_m[[m]]%*%s2_Ga_m)
    
    Q_fact[1:K,idx_m] <- Q_fact[1:K,idx_m] + crossprod(Lambda_m[[m]],s2_Ga_m)
    # Q_fact[idx_m,1:K] <- Q_fact[idx_m,1:K] + crossprod(Gamma_m[[m]],s2_La_m)
    Q_fact[idx_m,1:K] <- t(Q_fact[1:K,idx_m])
  }
  
  L_fact    <- t(chol(Q_fact))
  Lr_fact   <- forwardsolve(L_fact, r_fact)
  
  mean_fact <- backsolve(t(L_fact), Lr_fact)
  std_fact  <- backsolve(t(L_fact), matrix(rnorm((K+sum(K_Gm))*n),K+sum(K_Gm),n))
  
  new_fact  <- t(mean_fact + std_fact)
  
  eta_new <- new_fact[,1:K]
  
  phi_m_new <- vector("list", M)
  for(m in 1:M){
    idx_m <- K + c(0,cumsum(K_Gm)[-M])[m] + 1:K_Gm[m]
    phi_m_new[[m]] <- new_fact[,idx_m]
  }
  
  return(list(eta=eta_new, phi_m=phi_m_new))
}


update_factors_jfr <- function(n, M, K, res_m, s2_inv_m, Lambda_m){
  
  Q_fact <- diag(1,K,K)
  r_fact <- rep(0,K)
  
  for(m in 1:M){
    s2_La_m <- s2_inv_m[[m]]*Lambda_m[[m]]
    Q_fact  <- Q_fact + crossprod(Lambda_m[[m]],s2_La_m)
    r_fact  <- r_fact + t(res_m[[m]]%*%s2_La_m)
  }
  
  L_fact    <- t(chol(Q_fact))
  Lr_fact   <- forwardsolve(L_fact, r_fact)
  
  mean_fact <- backsolve(t(L_fact), Lr_fact)
  std_fact  <- backsolve(t(L_fact), matrix(rnorm(K*n),K,n))
  
  new_fact  <- t(mean_fact + std_fact)
  
  return(new_fact)
}


update_dcusp_joint <- function(M,K,logP_diff,pi_m,tens_I){
  
  delta_m <- matrix(0, nrow = M, ncol = K)
  
  # prior log probabilities
  tens_W = make_tensor_W(log(pi_m))    
  
  for(h in 1:K){
    
    # un-normalized probabilities
    tens_logP = tens_W
    for(m in 1:M){
      tens_logP = tens_logP + aperm(tens_I[,,,h],get_perm(m,M))*logP_diff[m,h]
    }
    
    # normalized probability matrix   
    pr_D <- c(exp(tens_logP - max(tens_logP)))
    pr_D <- pr_D / sum(pr_D)
    
    # sampling
    sampled_flat_index <- sample(x = 1:length(pr_D), size = 1, prob = pr_D)
    original_indices <- arrayInd(sampled_flat_index, .dim = dim(tens_logP))
    
    delta_m[,h] <- as.numeric(original_indices)
  }
  
  return(delta_m)
}

update_dcusp_seq <- function(M,K,logP_diff,pi_m,delta_m){
  
  # auxiliary indicators
  I_lh <- ( matrix(1:K,K,K) > matrix(1:K,K,K,byrow=T) )
  
  # random update order
  shuffle_m <- sample(1:M,M,replace=F)
  for(m in shuffle_m){
    
    # auxiliary indicators
    I_m_lh <- matrix(apply(delta_m[-m,,drop=F],2,max)>c(1:K),K,K,byrow=T)
    
    # un-normalized probabilities
    logP_D <- matrix(log(pi_m[m,]),K,K)
    for(mm in 1:M){
      if(mm==m){
        logP_D <- logP_D + matrix(logP_diff[m,],K,K,byrow=T)*I_lh*(1-(1-I_m_lh))
      } else {
        I_mm_lh <- (1-I_lh)
        if(M>2){
          I_mm_lh <- I_mm_lh * (1-matrix(apply(delta_m[-c(m,mm),,drop=F],2,max)>c(1:K),K,K,byrow=T))
        }
        logP_D <- logP_D + matrix(logP_diff[mm,]*(delta_m[mm,]>c(1:K)),K,K,byrow=T)*(1-I_mm_lh)
      }
    }
    
    # normalized probability matrix   
    pr_D   <- exp(logP_D - matrix(apply(logP_D, 2, max),K,K,byrow=T))
    pr_D   <- pr_D/matrix(apply(pr_D, 2, sum),K,K,byrow=T)
    
    # sampling
    #| Remarks:
    #|    in Hmisc::rMultinom(probs,nsample), the h-th row of `probs` gives the
    #|    probabilities for the l classes among which the h-th variable is sampled
    #|    implementation: pr_D[l,h] = P[\delta_m[m,h] = l]
    delta_m[m,]  <- as.vector(Hmisc::rMultinom(t(pr_D),1))
  }
  
  return(delta_m)
}

update_cusp <- function(K,logP_Spikes,logP_Slabs,pi_m){
  
  # un-normalized (log)-probability matrices
  lonP_N <- matrix(logP_Spikes,K,K) 
  logP_T <- matrix(logP_Slabs,K,K)
  
  # un-normalized (log)-probability matrix
  lonP_D <- lonP_N
  lonP_D[upper.tri(lonP_D,diag=F)] <- logP_T[upper.tri(lonP_D,diag=F)]
  lonP_D <- lonP_D + t(matrix(log(pi_m),K,K))
  
  # normalized probability matrix
  max_pD <- matrix(apply(lonP_D, 1, max),K,K)
  pr_D   <- exp(lonP_D - max_pD)
  pr_Tot <- apply(pr_D, 1, sum)
  pr_D   <- pr_D/pr_Tot
  
  # sampling
  #| Remarks:
  #|    in Hmisc::rMultinom(probs,nsample), the h-th row of `probs` gives the
  #|    probabilities for the l classes among which the h-th variable is sampled
  delta_m  <- as.vector(Hmisc::rMultinom(pr_D,1))
  
  return(delta_m)
}











