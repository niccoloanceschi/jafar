
is.scalar <- function(x) is.atomic(x) && length(x) == 1L

gibbs_JAFAR_CUSP <- function(y, X_m, n, M, p_m,                     # input data
                             nBurnIn=5000, nMCMC=2500, nThin=10,    # MCMC params
                             Kmax=NULL, Kmax_m=NULL,                # latent dim bounds
                             K0=NULL, K0_m=NULL,                    # initial number of factors
                             prec0=NULL, prec0m=NULL,               # intercept prior precisions
                             a_m=NULL, b_m=NULL,                    # idiosyncratic noises in predictors
                             a_sig=NULL, b_sig=NULL,                # response noise in response
                             a_theta=NULL, b_theta=NULL,            # slab in response loadings variances
                             a_chi=NULL, b_chi=NULL,                # slab in 'shared' loadings variances
                             a_tau=NULL, b_tau=NULL,                # slab in 'specific' loadings variances
                             var_spike_theta=NULL,                  # spike value in response loadings variances
                             a_xi=NULL, b_xi=NULL,                  # mixture weight in response loadings variances
                             var_spike=NULL,                        # spike value in loadings variances
                             alpha=NULL, alpha_loc=NULL,            # stick breaking params
                             t0=-1, t1=-5e-4, t0_adapt=20,          # adaptation
                             get_latent_vars=FALSE,                 # output loadings & latent factor
                             rescale_pred=FALSE,                    # rescale to cor in imputing NAs
                             get_last_sample=FALSE,                 # output full last sample
                             binary_y=FALSE,                        # is the response binary or continuous?
                             seed=123
                             ){
  
  set.seed(seed)
  
  # MCMC params
  nEff <- floor((nMCMC - nBurnIn)/nThin) 
  iter_print <- nMCMC %/% 10 
  teff <- 1
  
  t_delay <- floor(t0_adapt/2)
  
  # bounds on number of latent factors
  if(is.null(Kmax)){Kmax=max(min(floor(log(p_m)*3)),K0)}
  if(is.null(Kmax_m)){if(is.null(K0_m)){Kmax_m=floor(log(p_m)*3)
    } else {Kmax_m=pmin(floor(log(p_m)*3),K0_m)}
  } else if(is.scalar(Kmax_m)){Kmax_m=rep(Kmax_m,M)}
  
  # initial number of latent factors
  if(is.null(K0)){K0=Kmax-1} else {K0=min(K0,Kmax-1)}
  if(is.null(K0_m)){K0_m=Kmax_m-1
  } else if (is.scalar(K0_m)){K0_m=pmin(rep(K0_m,M),Kmax_m-1)
  } else {K0_m=pmin(K0_m,Kmax_m-1)}
  
  # hyperparams response noise
  if(is.null(a_sig)){a_sig=3}
  if(is.null(b_sig)){b_sig=1}
  
  # hyperparams idiosyncratic noises across modalities
  if(is.null(a_m)){a_m=rep(3,M)} else if(is.scalar(a_m)){a_m=rep(a_m,M)}
  if(is.null(b_m)){b_m=rep(1,M)} else if(is.scalar(b_m)){b_m=rep(b_m,M)}
  
  # hyperparams intercepts for responses & across modalities
  if(is.null(prec0)){prec0=1}
  if(is.null(prec0m)){prec0m=rep(1,M)} else if(is.scalar(prec0m)){prec0m=rep(prec0m,M)}
  
  # hyperparams slab in response loadings variances
  if(is.null(a_theta)){a_theta=0.5}
  if(is.null(b_theta)){b_theta=0.1}
  
  # hyperparams spike in response loadings variances
  if(is.null(var_spike_theta)){var_spike_theta=0.005}
  
  # hyperparams mixture weight in response loadings variances
  if(is.null(a_xi)){a_xi=3}
  if(is.null(b_xi)){b_xi=2}
  
  # hyperparam spike value in loadings variances
  if(is.null(var_spike)){var_spike=rep(0.005,M)} else if(is.scalar(var_spike)){var_spike=rep(var_spike,M)}
  if(is.null(var_spike)){var_spike=rep(0.005,M)} else if(is.scalar(var_spike)){var_spike=rep(var_spike,M)}
  
  # hyperparams slab in 'shared' loadings variances
  if(is.null(a_chi)){a_chi=rep(0.5,M)} else if(is.scalar(a_chi)){a_chi=rep(a_chi,M)}
  if(is.null(b_chi)){b_chi=rep(0.1,M)} else if(is.scalar(b_chi)){b_chi=rep(b_chi,M)}
  
  # hyperparams slab in 'specific' loadings variances
  if(is.null(a_tau)){a_tau=rep(0.5,M)} else if(is.scalar(a_tau)){a_tau=rep(a_tau,M)}
  if(is.null(b_tau)){b_tau=rep(0.1,M)} else if(is.scalar(b_tau)){b_tau=rep(b_tau,M)}
  
  # hyperparam beta dist stick breaking
  if(is.null(alpha)){alpha=rep(5,M)} else if(is.scalar(alpha)){alpha=rep(alpha,M)}
  if(is.null(alpha_loc)){alpha_loc=rep(5,M)} else if(is.scalar(alpha_loc)){alpha_loc=rep(alpha_loc,M)}
  
  # ------ Output Variables ------ #
  
  Theta_MC  <- array(0.,c(nEff,Kmax))
  
  if(get_latent_vars){
    eta_MC    <- array(0.,c(nEff,n,Kmax))
  }
  
  s2_inv_MC <- array(NA,c(nEff))
  mu_y_MC   <- array(NA,c(nEff))
  
  var_y_MC  <- array(NA,c(nEff))
  
  K_MC        <- array(NA,nEff) # overall NUMBER OF retained shared elements
  K_T_eff_MC  <- array(NA,nEff) # NUMBER OF active elements in Theta
  K_Lm_eff_MC <- array(NA,c(nEff,M)) # NUMBER OF active columns in Lambda_m
  
  K_Gm_MC     <- array(NA,c(nEff,M)) # n. of retained columns in Gamma_m
  K_Gm_eff_MC <- array(NA,c(nEff,M)) # NUMBER OF active columns in Gamma_m
  
  Lambda_m_MC <- Gamma_m_MC <- phi_m_MC <- s2_inv_m_MC <- mu_m_MC <- list()
  Marg_Var_m_MC <- Cov_m_mean <- list()
  
  active_T_MC <- array(0,c(nEff,Kmax))
  active_L_MC <- array(0,c(nEff,Kmax,M))
  
  for(m in c(1:M)){
    if(get_latent_vars){
      Lambda_m_MC[[m]] <- array(0.,c(nEff,p_m[m],Kmax))
      Gamma_m_MC[[m]]  <- array(0.,c(nEff,p_m[m],Kmax_m[m])) 
      phi_m_MC[[m]]    <- array(0.,c(nEff,n,Kmax_m[m]))
    }
    
    Marg_Var_m_MC[[m]] <- array(NA,c(nEff,p_m[m]))
    s2_inv_m_MC[[m]]   <- array(NA,c(nEff,p_m[m]))
    mu_m_MC[[m]]       <- array(0.,c(nEff,p_m[m]))
    
    Cov_m_mean[[m]]  <- matrix(0,p_m[m],p_m[m])
  }
  
  if(binary_y){ y_MC <- matrix(NA,nEff,n) }
  
  # ------ Latent Utilities & Missing Data --------------------------------------------------------

  # Re-name Binary Response
  if(binary_y){
    y_obs <- y
    
    left_thr = rep(-Inf,n)
    right_thr = rep(Inf,n)
    
    left_thr[which(y_obs>0)] <- 0
    right_thr[which(y_obs<1)] <- 0
  }
  
  # Check presence of NA in multi-view predictors
  NA_in_X <- max(sapply(X_m,function(df) max(is.na(df))))
  
  if(NA_in_X){
    # identifying indexes of NA
    Xm0   <- X_m
    Xm_na <- lapply(X_m,function(df) is.na(unname(as.matrix(df))))
    
    # Containers for output of imputed NA values
    na_idx     <- lapply(Xm_na, function(df) apply(df,1,which))
    na_row_idx <- lapply(na_idx,function(ll) c(1:n)[sapply(ll,length)>0])
    na_idx     <- lapply(1:M, function(m) na_idx[[m]][na_row_idx[[m]]])
    Xm_MC      <- lapply(na_idx, function(ll) lapply(ll, function(vec) matrix(NA,nEff,length(vec))))
    
    # Initial imputation from medians
    for(m in 1:M){
      Xm_medians <- apply(X_m[[m]],2,median,na.rm=T)
      for(idx in 1:length(na_idx[[m]])){
        X_m[[m]][na_row_idx[[m]][idx],unlist(na_idx[[m]][idx])] <- Xm_medians[unlist(na_idx[[m]][idx])]
      }
    }
    
  }
  
  # ------ Initialization ------------------------------------------------------
  
  par_init <- gibbs_JAFAR_CUSP_init(y, X_m, n, M, p_m,          # input data
                                    K0, K0_m,                   # initial number of factors
                                    a_sig, b_sig,               # response noise
                                    a_theta, b_theta,           # slab in response loadings variances
                                    var_spike_theta,            # spike value in response loadings variances
                                    a_xi, b_xi,                 # mixture weight in response loadings variances
                                    a_m, b_m,                   # idiosyncratic noises
                                    prec0, prec0m,              # intercepts
                                    var_spike, var_spike_vb,    # spike value in loadings variances
                                    a_chi, b_chi,               # slab in 'shared' loadings variances
                                    a_tau, b_tau,               # slab in 'specific' loadings variances
                                    alpha, alpha_loc,           # beta dist stick breaking
                                    seed)
  
  get_env = environment()
  list2env(par_init,envir=get_env)
  rm(par_init)
  
  s2_Ga_m <- D_m_chol <- res_m <- eta_LaT <- phi_GaT <- list()
  
  if(binary_y){ 
    s2_inv <- 1
    mu_y   <- qnorm(mean(y_obs))
  }

  Theta   <- rep(0,K)
  
  # ------ Gibbs Sampler Updates -----------------------------------------------
  
  # ------ print status ------- #
  print(sprintf(fmt = "%10s%3s%2s", "[",0,"%]"))
  
  for(t in c(1:nMCMC)){
    
    # STEP 0: Latent Utilities & Missing Data Imputation -----------------------

    if(binary_y){
      linPred <- rep(mu_y,n) + c(eta%*%Theta)
      # y <- linPred + (2*y_obs-1)*truncnorm::rtruncnorm(n, a=-(2*y_obs-1)*linPred, b=rep(Inf,n))
      y <- truncnorm::rtruncnorm(n, a=left_thr, b=right_thr, mean=linPred, sd=1)
    }
    
    if(NA_in_X & t>t_delay){
      for(m in 1:M){
        mar_std_m = rep(1,p_m[m])
        if(rescale_pred){
          mar_std_m = sqrt(1/s2_inv_m[[m]] + rowSums(Lambda_m[[m]]^2) + rowSums(Gamma_m[[m]]^2))
        }
        Ga_m = Gamma_m[[m]]/mar_std_m
        s2_m = s2_inv_m[[m]]*(mar_std_m^2)
        La_m = Lambda_m[[m]]/mar_std_m
        
        X_m_tmp <- matrix(mu_m[[m]],n,p_m[m],byrow=T) +
          tcrossprod(eta,La_m) + tcrossprod(phi_m[[m]],Ga_m) +
          t(matrix(rnorm(n*p_m[m]),p_m[m],n)/sqrt(s2_m))
        X_m[[m]][Xm_na[[m]]] <- X_m_tmp[Xm_na[[m]]]
      }
    }
    
    # STEP 1: Loadings ---------------------------------------------------------
    
    etaTeta <- crossprod(eta)
    
    ## 1.a Theta & mu_y --------------------------------------------------------
    
    etaTeta_mu <- rbind(c(n,colSums(eta)),cbind(colSums(eta),etaTeta))
    
    Q_Theta      <- diag(c(prec0,chi),K+1,K+1)+s2_inv*etaTeta_mu
    r_Theta      <- s2_inv*crossprod(cbind(rep(1,n),eta),y)
    
    L_Theta      <- t(chol(Q_Theta))
    Lr_Theta     <- forwardsolve(L_Theta, r_Theta)
    
    mean_Theta   <- backsolve(t(L_Theta), Lr_Theta)
    std_Theta    <- backsolve(t(L_Theta), rnorm(K+1))
    
    mu_y  <- as.vector(mean_Theta[1] + std_Theta[1])[1]
    Theta <- as.vector(mean_Theta[-1] + std_Theta[-1])
    
    ## 1.b Lambda_m, Gamma_m ---------------------------------------------------
    
    for(m in c(1:M)){
      
      facTfac <- matrix(0,K+K_Gm[m],K+K_Gm[m])
      etaTphi <- crossprod(eta,phi_m[[m]])
      phiTphi <- crossprod(phi_m[[m]])
      
      facTfac[c(1:K),c(1:K)]   <- etaTeta
      facTfac[-c(1:K),-c(1:K)] <- phiTphi
      facTfac[c(1:K),-c(1:K)]  <- etaTphi
      facTfac[-c(1:K),c(1:K)]  <- t(etaTphi)
      
      new_Loadings <- update_loadings(n, p_m[m], K, K_Gm[m], X_m[[m]], 
                                      facTfac, eta, phi_m[[m]], mu_m[[m]],
                                      s2_inv_m[[m]], chi_m[m,], tau_m[[m]])
      
      Lambda_m[[m]] <- new_Loadings[,c(1:K),drop=F]
      Gamma_m[[m]]  <- new_Loadings[,-c(1:K),drop=F]
      
      eta_LaT[[m]] <- tcrossprod(eta,Lambda_m[[m]])
      phi_GaT[[m]] <- tcrossprod(phi_m[[m]],Gamma_m[[m]])
    }
    
    # STEP 2: Intercepts -------------------------------------------------------
    
    ## 2.b mu_m ----------------------------------------------------------------
    
    if(t>t_delay){
      for(m in 1:M){
        vec_m <- colMeans(X_m[[m]]) - colMeans(eta_LaT[[m]]) - colMeans(phi_GaT[[m]])
        mean_mu_m <- vec_m*s2_inv_m[[m]]/(s2_inv_m[[m]]+prec0m[m]/n)
        mu_m[[m]] <- mean_mu_m + sqrt(1/(prec0m[m]+n*s2_inv_m[[m]]))*rnorm(p_m[m])
      }
    }
    
    # STEP 3: Idiosyncratic components -----------------------------------------
    
    ## 3.a s2_inv --------------------------------------------------------------
    if(!binary_y){
      res_y <- y-rep(mu_y,n)-eta%*%Theta
      s2_inv <- rgamma(1,shape=a_sig+0.5*n,rate=b_sig+0.5*sum(res_y^2))
    }
    
    ## 3.b s2_inv_m ------------------------------------------------------------
    for(m in c(1:M)){
      res_m[[m]] <- X_m[[m]]-matrix(mu_m[[m]],n,p_m[m],byrow=T)
      vec_m      <- eta_LaT[[m]] + phi_GaT[[m]]
      s2_inv_m[[m]] <- rgamma(p_m[m],shape=a_m[m]+0.5*n,rate=1) * 
        1/(b_m[m]+0.5*colSums((res_m[[m]]-vec_m)^2))
    }
    
    # STEP 4: Factors ----------------------------------------------------------
    
    # Shared quantities used for sampling both eta phi_m
    for(m in 1:M){
      s2_Ga_m[[m]]  <- s2_inv_m[[m]]*Gamma_m[[m]]
      D_m_chol[[m]] <- t(chol(diag(1.,K_Gm[m],K_Gm[m])+crossprod(Gamma_m[[m]],s2_Ga_m[[m]])))
    }
          
    ## 4.a eta -----------------------------------------------------------------
    
    Q_eta <- diag(1,K,K) + s2_inv*tcrossprod(Theta)
    r_eta <- s2_inv*tcrossprod(Theta,y-rep(mu_y,n))

    for(m in 1:M){
      
      Ga_T_s2_La_m     <- crossprod(s2_Ga_m[[m]],Lambda_m[[m]])
      Dinv_GaT_s2_La_m <- backsolve(t(D_m_chol[[m]]),forwardsolve(D_m_chol[[m]], Ga_T_s2_La_m))
      
      Q_eta <- Q_eta + crossprod(Lambda_m[[m]],s2_inv_m[[m]]*Lambda_m[[m]]) -
        crossprod(Ga_T_s2_La_m,Dinv_GaT_s2_La_m)
      r_eta <- r_eta + crossprod(Lambda_m[[m]],s2_inv_m[[m]]*t(res_m[[m]])) -
        t((res_m[[m]]%*%s2_Ga_m[[m]])%*%Dinv_GaT_s2_La_m)
    }
      
    L_eta    <- t(chol(Q_eta))
    Lr_eta   <- forwardsolve(L_eta, r_eta)
    
    mean_eta <- backsolve(t(L_eta), Lr_eta)
    std_eta  <- backsolve(t(L_eta), matrix(rnorm(K*n),K,n))
    
    eta      <- t(mean_eta + std_eta)

    ## 4.b phi_m ---------------------------------------------------------------
    for(m in c(1:M)){
      
      r_phi_m     <- t((res_m[[m]]-tcrossprod(eta,Lambda_m[[m]]))%*%s2_Ga_m[[m]])
      Lr_phi_m    <- forwardsolve(D_m_chol[[m]], r_phi_m)
      
      mean_phi_m  <- backsolve(t(D_m_chol[[m]]), Lr_phi_m)
      std_phi_m   <- backsolve(t(D_m_chol[[m]]), matrix(rnorm(K_Gm[m]*n),K_Gm[m],n))
      
      phi_m[[m]]  <- t(mean_phi_m + std_phi_m)
    }
      
    # STEP 5: Latent indicatirs ------------------------------------------------
    
    ## 5.a Shared Component ----------------------------------------------------
    
    # efficient computation of multivariate Normal & Student-T (log)-pdf
    
    logP_diff <- matrix(0,nrow=M+1,K)
    
    for(m in 1:M){
      vec_Lm <- colSums(Lambda_m[[m]]^2)
      lonP_Spikes <- -0.5*vec_Lm/var_spike[m] - rep(0.5*p_m[m]*log(2*pi*var_spike[m]),K) 
      logP_Slabs  <- -(0.5*p_m[m]+a_chi[m])*log(1+0.5*vec_Lm/b_chi[m]) +
        rep(lgamma(0.5*p_m[m]+a_chi[m]) - lgamma(a_chi[m]) - 0.5*p_m[m]*log(2*pi*b_chi[m]),K)
      logP_diff[m,] <- logP_Slabs - lonP_Spikes
    }
    
    lonP_Spikes <- - 0.5*Theta^2/var_spike_theta - rep(0.5*log(2*pi*var_spike_theta),K)
    logP_Slabs  <- - (0.5+a_theta)*log(1+0.5*Theta^2/b_theta) +
      rep(lgamma(0.5+a_theta) - lgamma(a_theta) - 0.5*log(2*pi*b_theta),K)

    logP_diff[M+1,] <- logP_Slabs - lonP_Spikes
    
    #| Remarks:
    #|     (i)  a priori: $ P[\delta_h = 0] = \xi $
    #|    (ii)  a priori: $ P[\delta_m[m,h] = l] = xi_{m l} $
    #|   (iii)  implementation: pr_D[l+1,h] = P[\delta_h = l]
    #|    (iv)  implementation: pr_D[l,h] = P[\delta_m[m,h] = l]

    ## 5.a.1 delta ------ #
    
    delta_max_mm <- matrix(0,M,K)
    for(m in 1:M){
      delta_max_mm[m,] <- apply(delta_m[-m,,drop=F],2,max) > c(1:K)
    }
    
    # un-normalized probabilities
    logP_D <- matrix(c(log(xi),log(1-xi)),2,K)
    logP_D[1,] <- logP_D[1,] +
      colSums(logP_diff[1:M,]*(delta_m > matrix(1:K,M,K,byrow=T))*delta_max_mm)
    logP_D[2,] <- logP_D[2,] + logP_diff[M+1,]*(apply(delta_m,2,max)>c(1:K)) +
      colSums(logP_diff[1:M,]*(delta_m > matrix(1:K,M,K,byrow=T)))
    
    # normalized probabilities
    pr_D <- exp(logP_D - matrix(apply(logP_D,2,max),2,K,byrow=T))
    pr_D <- pr_D[2,] / colSums(pr_D)
    
    # sampling 
    # delta  <- rbinom(K,1,pr_D)
    delta  <- c(rbinom(K-1,1,pr_D[-K]),0) # | enforcing last inactive entry in theta
    
    # identifying active components
    active_T <- which(delta==1)
    K_T_eff <- length(active_T)
    
    ## 5.a.2 delta_m ------ #
    
    # auxiliary indicators
    I_lh <- ( matrix(1:K,K,K) > matrix(1:K,K,K,byrow=T) )
    
    for(m in 1:M){
      
      # auxiliary indicators
      I_m_lh <- matrix(apply(delta_m[-m,,drop=F],2,max)>c(1:K),K,K,byrow=T)
      I_0_lh <- matrix(delta==1,K,K,byrow=T)
      
      # un-normalized probabilities
      logP_D <- matrix(log(xi_m[m,]),K,K)
      logP_D <- logP_D + matrix(logP_diff[M+1,],K,K,byrow=T)*I_0_lh*(1-(1-I_m_lh)*(1-I_lh))
      for(mm in 1:M){
        if(mm==m){
          logP_D <- logP_D + matrix(logP_diff[m,],K,K,byrow=T)*I_lh*(1-(1-I_m_lh)*(1-I_0_lh))
        } else {
          I_mm_lh <- (1-I_0_lh)*(1-I_lh)
          if(M>2){I_mm_lh <- I_mm_lh * (1-matrix(apply(delta_m[-c(m,mm),,drop=F],2,max)>c(1:K),K,K,byrow=T))}
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
      delta_m[m,]  <- as.vector(Hmisc::rMultinom(t(pr_D),1))
      
      # identifying active columns
      K_Lm_eff[m] <- 0
      active_L[[m]]  <- which( delta_m[m,] > c(1:K) ) # active columns within retained ones
      
      if(length(active_L[[m]])>0){K_Lm_eff[m] <- length(active_L[[m]])}
    }
    
    ## 5.b Specific Components -------------------------------------------------
    
    ### 5.b.1 z_m ---- #
    
    for(m in 1:M){
      
      if(K_Gm[m]>1){
        
        # efficient computation of multivariate Normal & Student-T (log)-pdf
        vec_G  <- colSums(Gamma_m[[m]][,1:K_Gm[m],drop=F]^2)
        lonP_N <- matrix( -0.5*vec_G/var_spike[m] - rep(0.5*p_m[m]*log(2*pi*var_spike[m]),K_Gm[m]) ,K_Gm[m],K_Gm[m]) 
        logP_T <- matrix( -(0.5*p_m[m]+a_tau[m])*log(1+0.5*vec_G/b_tau[m]) +
                            rep(lgamma(0.5*p_m[m]+a_tau[m]) - lgamma(a_tau[m]) - 0.5*p_m[m]*log(2*pi*b_tau[m]),K_Gm[m]) ,K_Gm[m],K_Gm[m])
        
        # un-normalized (log)-probability matrix
        lonP_Z <- lonP_N
        lonP_Z[upper.tri(lonP_Z,diag=F)] <- logP_T[upper.tri(lonP_Z,diag=F)]
        lonP_Z <- lonP_Z + t(matrix(log(w_m[[m]]),K_Gm[m],K_Gm[m]))
        
        # normalized probability matrix
        max_pZ <- matrix(apply(lonP_Z, 1, max),K_Gm[m],K_Gm[m])
        pr_z   <- exp(lonP_Z - max_pZ)
        pr_Tot <- apply(pr_z, 1, sum)
        pr_z   <- pr_z/pr_Tot
        
        # sampling
        z_m[[m]]  <- as.vector(Hmisc::rMultinom(pr_z,1))
      } else {
        z_m[[m]]  <- as.vector(1)
      }
      
      # identifying active columns
      K_Gm_eff[m] <- 0
      active_G[[m]]  <- which( z_m[[m]] > c(1:K_Gm[m]) ) # active columns among retained ones
      if(length(active_G[[m]])>0){K_Gm_eff[m] <- length(active_G[[m]])}
    }
    
    # STEP 6: Stick Breaking Elements ------------------------------------------
    
    ## 6.a Theta ---------------------------------------------------------------
    
    ### 6.a.1 xi (spike and slab mixture weight) ----
    # Remark: a priori $ \xi = P[delta_h = 0] $
    
    xi <- rbeta(1, shape1 = a_xi+K-sum(delta), shape2 = b_xi+sum(delta))
    
    ### 6.a.2 chi (response loadings precisions) ----
    chi <- rep(1./var_spike_theta,K)
    if(length(active_T)>0){
      chi[active_T] <- rgamma(length(active_T),shape=a_theta+0.5,rate=1) *
        1./(b_theta + 0.5*Theta[active_T]^2)
    }
    
    ## 6.b Lambda_m ------------------------------------------------------------
    
    ### 6.b.1 rho_m, xi_m (stick breaking elements) ----
    for(m in 1:M){
      count_eq <- colSums(as.matrix(delta_m[m,] == matrix(1:K,K,K,byrow=T)))
      count_gr <- rev(cumsum(rev(c(count_eq[-1],0))))
      
      rho_m[m,] <- c(rbeta(K-1, shape1 = 1+count_eq[-K], shape2 = alpha[m]+count_gr[-K]),1.)
      xi_m[m,]  <- rho_m[m,]*c(1,cumprod(1-rho_m[m,-K]))
    }  
    
    ### 6.b.2 chi_m (loading precisions) ----
    for(m in 1:M){
      chi_m[m,] <- rep(1./var_spike[m],K)
      if(length(active_L[[m]])>0){
        
        chi_m[m,][active_L[[m]]] <- rgamma(length(active_L[[m]]),shape=a_chi[m]+0.5*p_m[m],rate=1) *
          1./(b_chi[m] + 0.5 * colSums(Lambda_m[[m]][,active_L[[m]],drop=F]^2))
      }
    }
    
    ## 6.c Gamma_m -------------------------------------------------------------
    
    ### 6.c.1 nu_m, w_m (stick breaking elements) ----
    for(m in 1:M){
      count_eq <- colSums(as.matrix(z_m[[m]] == t(c(1:K_Gm[m])*matrix(1,K_Gm[m],K_Gm[m]))))
      count_gr <- rev(cumsum(rev(c(count_eq[-1],0))))
      
      nu_m[[m]] <- c(rbeta(K_Gm[m]-1, shape1 = 1+count_eq[-K_Gm[m]], shape2 = alpha_loc[m]+count_gr[-K_Gm[m]]),1.)
      w_m[[m]]  <- nu_m[[m]]*c(1,cumprod(1-nu_m[[m]][-K_Gm[m]]))
    }
    
    ### 6.c.2 tau_m (loading precisions) ----
    for(m in 1:M){
      tau_m[[m]] <- rep(1./var_spike[m],K_Gm[m])
      if(length(active_G[[m]])>0){
        
        tau_m[[m]][active_G[[m]]] <- rgamma(length(active_G[[m]]),shape=a_tau[m]+0.5*p_m[m],rate=1) *
          1./(b_tau[m] + 0.5 * colSums(Gamma_m[[m]][,active_G[[m]],drop=F]^2))
      }
    }
    
    # STEP 7: Adaptation -------------------------------------------------------
    
    if((t>t0_adapt) & (runif(1) < exp(t0 + t1*t))){
      
      ## 7.a K -----------------------------------------------------------------
      
      active_J <- which(colSums(rbind(delta_m > matrix(1:K,M,K,byrow=T),delta > 0)) > 1)
      K_eff    <- length(active_J)
      
      if(K_eff < K-1){
        
        K <- K_eff + 1
        
        eta   <- cbind(eta[,active_J,drop=F],rnorm(n))
        
        for(m in 1:M){
          Lambda_m[[m]] <- cbind(Lambda_m[[m]][,active_J,drop=F],sqrt(var_spike[m])*rnorm(p_m[m]))
        }
        delta_m <- cbind(delta_m[,active_J,drop=F],rep(1,M))
        rho_m   <- cbind(rho_m[,active_J,drop=F],rep(1,M))
        xi_m    <- cbind(xi_m[,active_J,drop=F],rep(1,M)-rowSums(xi_m[,active_J,drop=F]))
        chi_m   <- cbind(chi_m[,active_J,drop=F],1./var_spike)
        
        Theta <- c(Theta[active_J],sqrt(var_spike_theta)*rnorm(1))
        delta <- c(delta[active_J],0)
        chi   <- c(chi[active_J],1./var_spike_theta)
        
      } else if (K < Kmax){
        
        K <- K + 1
        
        eta   <- cbind(eta,rnorm(n))
        
        for(m in 1:M){
          Lambda_m[[m]] <- cbind(Lambda_m[[m]],sqrt(var_spike[m])*rnorm(p_m[m]))
        }
        delta_m <- cbind(delta_m,rep(1,M))
        rho_m   <- cbind(rho_m[,-(K-1)],rbeta(1,shape1=1,shape2=alpha),rep(1,M))
        xi_m    <- rho_m * cbind(rep(1,M),t(apply(1-rho_m[,-K],1,cumprod)))
        chi_m   <- cbind(chi_m,1./var_spike)
        
        Theta <- c(Theta,sqrt(var_spike_theta)*rnorm(1))
        delta <- c(delta,0)
        chi   <- c(chi,1./var_spike_theta)
        
      }
      
      # 7.b K_Gm --------------------------------------------------------------
      
      for(m in 1:M){
        
        if(K_Gm_eff[m]==0){
          
          K_Gm[m] <- 1
          
          phi_m[[m]]   <- matrix(rnorm(n),n,1)
          
          Gamma_m[[m]] <- matrix(sqrt(var_spike[m])*rnorm(p_m[m]),p_m[m],1)
          w_m[[m]]   <- c(1.)
          tau_m[[m]] <- c(1./var_spike[m])
          
        } else if (K_Gm_eff[m] < K_Gm[m]-1) {
          
          K_Gm[m] <- K_Gm_eff[m] + 1
          
          phi_m[[m]]   <- cbind(phi_m[[m]][,active_G[[m]],drop=F],rnorm(n))
          
          Gamma_m[[m]] <- cbind(Gamma_m[[m]][,active_G[[m]],drop=F],sqrt(var_spike[m])*rnorm(p_m[m]))
          w_m[[m]]   <- c(w_m[[m]][active_G[[m]]],1-sum(w_m[[m]][active_G[[m]]]))
          tau_m[[m]] <- c(tau_m[[m]][active_G[[m]]],1./var_spike[m])
          
        } else if (K_Gm[m] < Kmax_m[m]) {
          
          K_Gm[m] <- K_Gm[m] + 1
          
          phi_m[[m]]   <- cbind(phi_m[[m]],rnorm(n))
          
          Gamma_m[[m]] <- cbind(Gamma_m[[m]],sqrt(var_spike[m])*rnorm(p_m[m]))
          nu_m[[m]]  <- c(nu_m[[m]][-(K_Gm[m]-1)],rbeta(1,shape1=1,shape2=alpha_loc[m]),1.)
          w_m[[m]]   <- nu_m[[m]]*c(1,cumprod(1-nu_m[[m]][-K_Gm[m]]))
          tau_m[[m]] <- c(tau_m[[m]],1./var_spike[m])
          
        }
      }
    }
    
    # ------ saving outputs (after thinning) ------ #
    if((t %% nThin == 0) & (t > nBurnIn)) {
      
      Theta_MC[teff,1:K] <- Theta
      
      active_T_MC[teff,active_T] <- 1
      
      K_T_eff_MC[teff] <- K_T_eff
      
      if(get_latent_vars){
        eta_MC[teff,,1:K]  <- eta
      }
      
      s2_inv_MC[teff] <- s2_inv
      mu_y_MC[teff]   <- mu_y
      
      var_y_MC[teff]  <- 1/s2_inv + sum(Theta^2) 
      
      K_MC[teff]         <- K # overall NUMBER of retained shared columns
      K_Lm_eff_MC[teff,] <- K_Lm_eff # NUMBER of active columns in Lambda_m
      
      K_Gm_MC[teff,]     <- K_Gm     # n. of retained columns in Gamma_m
      K_Gm_eff_MC[teff,] <- K_Gm_eff # NUMBER of active columns in Gamma_m
      
      for(m in c(1:M)){
        if(get_latent_vars){
          Lambda_m_MC[[m]][teff,,1:K] <- Lambda_m[[m]]
          Gamma_m_MC[[m]][teff,,1:K_Gm[m]]  <- Gamma_m[[m]]
          phi_m_MC[[m]][teff,,1:K_Gm[m]]    <- phi_m[[m]]
        }
        
        active_L_MC[teff,active_L[[m]],m] <- 1
        
        Marg_Var_m_MC[[m]][teff,] <- 1/s2_inv_m[[m]] + rowSums(Lambda_m[[m]]^2) + rowSums(Gamma_m[[m]]^2)
        s2_inv_m_MC[[m]][teff,]   <- s2_inv_m[[m]]
        mu_m_MC[[m]][teff,]       <- mu_m[[m]]
        
        Cov_m_mean[[m]] <- Cov_m_mean[[m]] +
          ( tcrossprod(Lambda_m[[m]]) + tcrossprod(Gamma_m[[m]]) + diag(1/s2_inv_m[[m]]) ) / nEff
      }
      
      if(NA_in_X){
        for(m in 1:M){
          for(idx in 1:length(na_idx[[m]])){
            Xm_MC[[m]][[idx]][teff,] <- X_m[[m]][na_row_idx[[m]][idx],unlist(na_idx[[m]][idx])]
          }
        }
      }

      if(binary_y){y_MC[teff,] <- y}
      
      teff = teff + 1
    }
    
    # ------ print status ------- #
    if(t %% iter_print == 0){
      print(sprintf(fmt = "%10s%3s%2s", "[",(t%/%iter_print)*10,"%]"))
    }
  }
  
  hyper_params = list(model='jafar_joint_cusp', 
                      rescale_pred=rescale_pred, 
                      vb_init=vb_init, vb_rank_only=vb_rank_only,
                      nBurnIn=nBurnIn, nMCMC=nMCMC, nThin=nThin, 
                      Kmax=Kmax, Kmax_m=Kmax_m, K0=K0, K0_m=K0_m, 
                      seed=seed, t0=t0, t1=t1, t0_adapt=t0_adapt,
                      a_theta=a_theta, b_theta=b_theta, a_xi=a_xi,        
                      b_xi=b_xi, var_spike_theta=var_spike_theta, 
                      a_sig=a_sig, b_sig=b_sig, prec0=prec0, 
                      a_m=a_m, b_m=b_m, prec0m=prec0m,
                      var_spike=var_spike, var_spike_vb=var_spike_vb,
                      a_chi=a_chi, b_chi=b_chi, alpha=alpha,
                      a_tau=a_tau, b_tau=b_tau, alpha_loc=alpha_loc)
  
  output = list(K=K_MC,K_T_eff=K_T_eff_MC,K_Lm_eff=K_Lm_eff_MC,
                K_Gm=K_Gm_MC,K_Gm_eff=K_Gm_eff_MC,
                Theta=Theta_MC,var_y=var_y_MC,Cov_m_mean=Cov_m_mean,
                s2_inv=s2_inv_MC,mu_y=mu_y_MC,s2_inv_m=s2_inv_m_MC,mu_m=mu_m_MC,
                Marg_Var_m=Marg_Var_m_MC,hyper_param=hyper_params,
                active_Lm=active_L_MC,active_T=active_T_MC)

  if(binary_y){output$y_lat <- y_MC}
  
  if(NA_in_X){
    output$Xm_MC=Xm_MC
    output$na_idx=na_idx
    output$na_row_idx=na_row_idx
  }
  
  if(get_latent_vars){
    output$eta=eta_MC
    output$Lambda_m=Lambda_m_MC
    output$Gamma_m=Gamma_m_MC
    output$phi_m=phi_m_MC
  }
  
  if(get_last_sample){
    output$last_sample = list(Theta=Theta,eta=eta,s2_inv=s2_inv,mu_y=mu_y,
                              chi=chi,delta=delta,xi=xi,Lambda_m=Lambda_m,
                              Gamma_m=Gamma_m,phi_m=phi_m,s2_inv_m=s2_inv_m,
                              mu_m=mu_m,z_m=z_m,nu_m=nu_m,w_m=w_m,tau_m=tau_m,
                              delta_m=delta_m,rho_m=rho_m,xi_m=xi_m,chi_m=chi_m)
  }
  
  return(output)
}
              











