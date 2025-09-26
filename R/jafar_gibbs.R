
#' Gibbs Sampler for JAFAR
#'
#' @description
#' Fits a Joint Additive FActor Regression (JAFAR) model using Gibbs sampling.
#' Variation across multiple data-views is explained via shared and study-specific latent factors.
#' Default and optional outputs include posterior means of the induced covariances,
#' posterior samples of residual variances, latent factors, and factor loadings. 
#' Supports parallel computation and tempered updates to limit rank estimation in extreme large-p-small-n settings.
#'
#' @param X_m Multi-view input data. Rows should correspond to samples, columns to features.
#'  (list of length \code{M}; m-th element: matrix n x p_m[m]).
#' @param y Vector responses (of length n). Set to NULL for unsupervised mode (default: NULL).
#' @param yBinary Logical, indicating if the response(s) are binary (default: FALSE).
#' @param K0 Upper bound to numbers of shared latent factors (optional)
#'  If NULL, \code{K0} is set to \code{3*log(max(p_m))}
#' @param K0_m Upper bounds to numbers of view-specific latent factors (optional)
#'  Length should equal \code{length(X_m)}. If NULL, \code{K0[m]} is set to \code{3*log(max_(p_m[m]))}
#' @param tMCMC Total number of MCMC iterations (default: 20000).
#' @param tBurnIn Number of burn-in iterations (default: 15000).
#' @param tThin Thinning interval for saving samples (default: 10).
#' @param hyperparams List of hyperparameters for the D-CUSP prior distributions.
#'  Missing hyperparameters are replaced by defaults encoded in \code{\link{jafar_set_hyperparameters}}.
#' @param get_latent_vars Logical, whether to return latent factors and loading matrices (default: TRUE).
#' @param get_last_sample Logical, whether to return only the last sample of the MCMC chain (default: FALSE).
#' @param parallel Logical, whether to use parallel computation for the loadings' update (default: TRUE).
#' @param tempered Logical, temperature parameter for tempered sampling (default: FALSE, no tempering).
#' @param rescale_pred Logical, whether to rescale loadings when computing response predictions (default: FALSE).
#'
#' @note
#' - Ensure that all matrices in \code{X_m} have the same number of rows (subjects).  
#' - Missing data in \code{X_m} are allowed as NA and imputed in the MCMC.
#'
#' @return A list containing posterior samples, latent variables (if requested), and other relevant model outputs.
#'
#' @details
#'  The number of samples in output is \code{tEff=(tMCMC-tBurnIn)%/%tThin}.
#'  The output list contains:
#'  \itemize{
#'    \item{\code{K}}{Number shared latent factors (vector of length \code{tEff}).}
#'    \item{\code{K_Gm}}{Number view-specific latent factors (matrix \code{tEff x M}).}
#'    \item{\code{K_Lm_eff}}{Numbers of shared factors active in each view (matrix \code{tEff x M}).}
#'    \item{\code{K_Gm_eff}}{Numbers of specific factors active in each view (matrix \code{tEff x M}).}
#'    \item{\code{active_Lm}}{Binary indicators of shared factors activity across views (binary array \code{tEff x K x M}s).}
#'    \item{\code{Cov_m_mean}}{Posterior mean of the covariance matrix for each dataset (list of length \code{M}; m-th element: matrix \code{p_m[m] x p_m[m]}).}
#'    \item{\code{Marg_Var_m}}{Marginal variances of features (list of length \code{M}; m-th element: matrix \code{tEff x p_m[m]}).}
#'    \item{\code{s2_inv_m}}{Inverse residual variances across views (list of length \code{M}; m-th element: matrix \code{tEff x p_m[m]}).}
#'    \item{\code{mu_m}}{Features intercepts across views (list of length \code{M}; m-th element: matrix \code{tEff x p_m[m]}).}
#'    \item{\code{hyper_param}}{List of hyperparameters used for the model, including user-specified values and defaults ones were missing.}
#'  }
#'  If \code{is_supervised = TRUE}:
#'  \itemize{   
#'    \item{\code{K_T_eff}}{Numbers of shared factors active in the response (vector of length \code{tEff}).}
#'    \item{\code{K_Tm_eff}}{Numbers of specific factors active in the response (matrix \code{tEff x M}).}
#'    \item{\code{active_T}}{Binary indicators of shared factors activity in the response (binary matrix \code{tEff x K}).}
#'    \item{\code{active_Tm}}{Binary indicators of specific factors activity in the response (list of length \code{M}; m-th element: matrix \code{tEff x K_Gm[m]}).}
#'    \item{\code{s2_inv}}{Response inverse residual variances (vector of length \code{tEff}).}
#'    \item{\code{mu_y}}{Response intercept (vector of length \code{tEff}).}
#'    \item{\code{Theta}}{Response loadings on shared factors (matrix \code{tEff x K}).}
#'    \item{\code{Theta_m}}{Response loadings on specific factors (list of length \code{M}; m-th element: matrix \code{tEff x K_Gm[m]}).}
#'    \item{\code{y_MC}}{Latent probit utilities (matrix \code{tEff x n}). (only if \code{yBinary = TRUE}).}
#'  }
#'  If \code{get_latent_vars = TRUE}:
#'  \itemize{   
#'    \item{\code{Lambda_m}}{Loadings matrices on shared factors (list of length \code{M}; m-th element: array \code{tEff x p_m[m] x K}).}
#'    \item{\code{Gamma_m}}{Loadings matrices on view-specific factors (list of length \code{M}; m-th element: array \code{tEff x p_m[m] x K_Gm[m]}).}
#'    \item{\code{eta}}{Shared latent factors (array \code{tEff x n x K}).}
#'    \item{\code{phi_m}}{View-specific latent factors (list of length \code{M}; m-th element: array \code{tEff x n x K_Gm[m]}).}
#'  }
#'  If the input matrices \code{X_m} contain missing values:
#'  \itemize{
#'    \item{\code{Xm_MC}}{Posterior samples of imputed values for missing entries.
#'      A list of length \code{M}; the m-th element is itself a list (one per feature with missingness),
#'      each containing an \code{tEff × n_miss} matrix of imputed values across MCMC iterations.}
#'    \item{\code{na_idx}}{List of length \code{M}; the m-th element gives the column indices of missing entries in \code{X_m[[m]]}.}
#'    \item{\code{na_row_idx}}{List of length \code{M}; the m-th element gives the corresponding row indices of missing entries in \code{X_m[[m]]}.}
#' }
#'  If \code{get_last_sample = TRUE}:
#'  \itemize{   
#'    \item{\code{last_sample}}{List of posterior values of all parameters at the last MCMC iteration, including latent factors, loadings, residual variances, and hyperparameters.}
#' }
#'
#' @export
#'
gibbs_jafar <- function(X_m, y=NULL, yBinary=F, K0=NULL, K0_m=NULL, 
                        tMCMC=20000, tBurnIn=15000, tThin=10, 
                        hyperparams = list(),
                        get_latent_vars=TRUE, get_last_sample=FALSE,
                        parallel=TRUE, tempered=FALSE, 
                        rescale_pred=FALSE){
  
  get_env <- environment()
  
  is_supervised <- !is.null(y)
  if(is_supervised){get_latent_vars=T}
  
  which_prior='d-cusp'
  print(paste0('Running ',c('un-', '')[1+is_supervised]
               ,'supervised jafar under ',which_prior,' prior'))
  
  M = length(X_m) # number of modalities
  p_m = sapply(X_m,ncol) # number of features per modality
  n = nrow(X_m[[1]]) # number of samples
  
  if(M<2){stop('jafar requires at least 2 data views. Fir jfr if you have just one.')}
  
  # MCMC parameters
  tEff <- floor((tMCMC - tBurnIn)/tThin) 
  iter_print <- tMCMC %/% 10 
  teff <- 1
  
  # Latent dimensions
  if(is.null(K0)){K0<-floor(3*log(max(p_m)))}
  if(is.null(K0_m)){K0_m <- floor(3*log(p_m))}
  if(is.scalar(K0_m)){K0_m <- rep(K0_m)}
  
  # Set hyperparameters
  hyperparams <- jafar_set_hyperparameters(hyperparams, M, is_supervised)
  list2env(hyperparams, envir=get_env)
  set.seed(seed)
  
  # Output Variables
  output_containers <- jafar_initialize_output(tEff, n, K0, K0_m, M, p_m,
                                               get_latent_vars, is_supervised, yBinary)
  list2env(output_containers, envir=get_env)
  rm(output_containers)
  K_MC <- rep(NA,tMCMC%/%tThin)
  K_Gm_MC <- matrix(NA,tMCMC%/%tThin,M)
  
  # Initialization 
  par_init <- jafar_initialize_sampler(n, M, p_m, K0, K0_m, hyperparams,
                                       is_supervised, yBinary)
  list2env(par_init,envir=get_env)
  rm(par_init)
  
  # Auxiliary variables
  res_m <- eta_LaT <- phi_GaT <- list()
  
  # Likelihood Tempering Power
  inv_pow_Gm <- inv_pow_Lm <- rep(1,M)
  if(tempered){
    inv_pow_Gm <- pmin(rep(n/M,M),p_m)/p_m
    inv_pow_Lm <- pmin(rep(n,M),p_m)/p_m
  }
  
  # Rescaling DP concentrations
  alpha_L = alpha_L / sqrt(M)
  alpha_G = alpha_G / sqrt(M)
    
  # ------ Latent Utilities & Missing Data --------------------------------------------------------
  
  # Check presence of NA in omics data
  impute_na <- max(sapply(X_m,function(df) max(is.na(df))))
  
  if(impute_na){
    # identifying indexes of NA
    Xm0   <- X_m
    Xm_na <- lapply(X_m,function(df) is.na(unname(as.matrix(df))))
    
    # Containers for output of imputed NA values
    na_idx     <- lapply(Xm_na, function(df) apply(df,1,which))
    na_row_idx <- lapply(na_idx,function(ll) c(1:n)[sapply(ll,length)>0])
    na_idx     <- lapply(1:M, function(m) na_idx[[m]][na_row_idx[[m]]])
    Xm_MC      <- lapply(na_idx, function(ll) lapply(ll, function(vec) matrix(NA,tEff,length(vec))))
    
    # Initial imputation from medians
    for(m in 1:M){
      Xm_medians <- apply(X_m[[m]],2,median,na.rm=T)
      for(idx in 1:length(na_idx[[m]])){
        X_m[[m]][na_row_idx[[m]][idx],unlist(na_idx[[m]][idx])] <- Xm_medians[unlist(na_idx[[m]][idx])]
      }
    }
  }
  
  if(is_supervised & yBinary){
     y_obs <- y
     
     left_thr = rep(-Inf,n)
     right_thr = rep(Inf,n)
     
     left_thr[which(y_obs>0)] <- 0
     right_thr[which(y_obs<1)] <- 0
  }
  
  # ------ Gibbs Sampler Updates -----------------------------------------------
  
  # ------ print status ------- #
  cat(sprintf("%10s[%3d%%] K=%02d %s\n", "", 0, K,
      paste(sprintf("K%d=%02d",1:M,K_Gm),collapse=" ")))
  
  for(t in c(1:tMCMC)){
    
    # 0: Latent Utilities & Missing Data Imputation ----------------------------
    
    if(impute_na & t>t0_adapt){
      for(m in 1:M){
        mar_std_m = rep(1,p_m[m])
        if(rescale_pred){
          mar_std_m = sqrt(1/s2_inv_m[[m]] + rowSums(Lambda_m[[m]]^2) + rowSums(Gamma_m[[m]]^2))
        }
        s2_m = s2_inv_m[[m]]*(mar_std_m^2)
        La_m = Lambda_m[[m]]/mar_std_m
        Ga_m = Gamma_m[[m]]/mar_std_m
        
        X_m_tmp <- matrix(mu_m[[m]],n,p_m[m],byrow=T) +
          t(matrix(rnorm(n*p_m[m]),p_m[m],n)/sqrt(s2_m)) + 
          tcrossprod(eta,La_m) + tcrossprod(phi_m[[m]],Ga_m)
        X_m[[m]][Xm_na[[m]]] <- X_m_tmp[Xm_na[[m]]]
      }
    }
    
    if(is_supervised & yBinary){
      Theta_all = c(Theta,unlist(Theta_m))
      linPred <- rep(mu_y,n) + as.vector(cbind(eta,do.call(cbind, phi_m))%*%Theta_all)
      # y <- linPred + (2*y_obs-1)*truncnorm::rtruncnorm(n, a=-(2*y_obs-1)*linPred, b=rep(Inf,n))
      y <- truncnorm::rtruncnorm(n, a=left_thr, b=right_thr, mean=linPred, sd=1)
    }
    
    # 1: Loadings --------------------------------------------------------------
    
    etaTeta <- crossprod(eta)
    
    for(m in c(1:M)){
      
      facTfac <- matrix(0,K+K_Gm[m],K+K_Gm[m])
      etaTphi <- crossprod(eta,phi_m[[m]])
      phiTphi <- crossprod(phi_m[[m]])
      
      facTfac[c(1:K),c(1:K)]   <- etaTeta
      facTfac[-c(1:K),-c(1:K)] <- phiTphi
      facTfac[c(1:K),-c(1:K)]  <- etaTphi
      facTfac[-c(1:K),c(1:K)]  <- t(etaTphi)
      
      if(parallel){
        new_Loadings <- update_loadings_parallel(n, p_m[m], K+K_Gm[m], X_m[[m]], facTfac, 
                                                cbind(eta,phi_m[[m]]), mu_m[[m]],
                                                s2_inv_m[[m]], c(chi_m[m,],tau_m[[m]]))
      } else {
        new_Loadings <- update_loadings(n, p_m[m], K+K_Gm[m], X_m[[m]], facTfac, 
                                        cbind(eta,phi_m[[m]]), mu_m[[m]],
                                        s2_inv_m[[m]], c(chi_m[m,],tau_m[[m]]))
      }
      
      Lambda_m[[m]] <- new_Loadings[,c(1:K),drop=F]
      Gamma_m[[m]]  <- new_Loadings[,-c(1:K),drop=F]
      
      eta_LaT[[m]] <- tcrossprod(eta,Lambda_m[[m]])
      phi_GaT[[m]] <- tcrossprod(phi_m[[m]],Gamma_m[[m]])
    }
    
    if(is_supervised){
      fac_mu <- cbind(rep(1,n),eta,do.call(cbind, phi_m))
      facTfac_mu <- crossprod(fac_mu)
      
      Q_Theta      <- s2_inv*facTfac_mu + diag(c(prec0,psi,unlist(psi_m)),1+K+sum(K_Gm),1+K+sum(K_Gm))
      r_Theta      <- s2_inv*crossprod(fac_mu,y)
      
      L_Theta      <- t(chol(Q_Theta))
      Lr_Theta     <- forwardsolve(L_Theta, r_Theta)
      
      mean_Theta   <- backsolve(t(L_Theta), Lr_Theta)
      std_Theta    <- backsolve(t(L_Theta), rnorm(1+K+sum(K_Gm)))
      
      mu_y      <- as.vector(mean_Theta[1] + std_Theta[1])[1]
      Theta_all <- as.vector(mean_Theta[-1] + std_Theta[-1])
      
      Theta     <- Theta_all[1:K]
      Theta_m   <- lapply(1:M, function(m) Theta_all[K+c(0,cumsum(K_Gm)[-M])[m]+1:K_Gm[m]])
    }
    
    # 2: Intercepts ------------------------------------------------------------
    
    if(t>t0_adapt){
      for(m in 1:M){
        vec_m <- colMeans(X_m[[m]]) - colMeans(eta_LaT[[m]]) - colMeans(phi_GaT[[m]])
        mean_mu_m <- vec_m*s2_inv_m[[m]]/(s2_inv_m[[m]]+prec0m[m]/n)
        mu_m[[m]] <- mean_mu_m + sqrt(1/(prec0m[m]+n*s2_inv_m[[m]]))*rnorm(p_m[m])
      }
    }
    
    # 3: Idiosyncratic components ----------------------------------------------
    
    for(m in c(1:M)){
      res_m[[m]] <- X_m[[m]]-matrix(mu_m[[m]],n,p_m[m],byrow=T)
      vec_m      <- eta_LaT[[m]] + phi_GaT[[m]]
      s2_inv_m[[m]] <- rgamma(p_m[m],shape=a_m[m]+0.5*n,rate=1) * 
        1/(b_m[m]+0.5*colSums((res_m[[m]]-vec_m)^2))
    }
    
    if(is_supervised & !yBinary){
      res_y <- y - rep(mu_y,n) - as.vector(cbind(eta,do.call(cbind, phi_m))%*%Theta_all)
      s2_inv <- rgamma(1,shape=a_sig+0.5*n,rate=b_sig+0.5*sum(res_y^2))
    }

    # 4: Factors ---------------------------------------------------------------

    if(!is_supervised){
      fact_new <- update_factors_collapsed(n, M, K, K_Gm, res_m, s2_inv_m, Lambda_m, Gamma_m)
    } else {
      fact_new <- update_factors_supervised(n, M, K, K_Gm, res_m, s2_inv_m, Lambda_m, Gamma_m,
                                            y-rep(mu_y,n), s2_inv, Theta_all)
    }
    eta <- fact_new$eta
    phi_m <- fact_new$phi_m
    
    # 5: Latent indicators -----------------------------------------------------

    ## 5.1: Shared Component ---------------------------------------------------

    # efficient computation of multivariate Normal & Student-T (log)-pdf

    logP_diff <- logP_Spikes <- logP_Slabs <- matrix(0,nrow=M,ncol=K)

    for(m in 1:M){
      vec_Lm <- colSums(Lambda_m[[m]]^2)
      logP_Spikes[m,] <- -0.5*vec_Lm/var_spike[m] - rep(0.5*p_m[m]*log(2*pi*var_spike[m]),K)
      logP_Slabs[m,]  <- -(0.5*p_m[m]+a_chi[m])*log(1+0.5*vec_Lm/b_chi[m]) +
        rep(lgamma(0.5*p_m[m]+a_chi[m]) - lgamma(a_chi[m]) - 0.5*p_m[m]*log(2*pi*b_chi[m]),K)
      
      # Likelihood Tempering Power
      logP_Spikes[m,] <- logP_Spikes[m,] * inv_pow_Lm[m]
      logP_Slabs[m,] <- logP_Slabs[m,] * inv_pow_Lm[m]
      
      logP_diff[m,] <- logP_Slabs[m,] - logP_Spikes[m,]
    }

    #| Remark:
    #|    (i)   a priori: $ P[\delta_Lm[m,h] = l] = xi_{m l} $

    delta_Lm <- update_dcusp_seq(M,K,logP_diff,pi_m,delta_Lm)

    # identifying active columns
    for(m in 1:M){
      K_Lm_eff[m]   <- 0
      active_L[[m]] <- which( delta_Lm[m,] > c(1:K) ) # active columns within retained ones
      if(length(active_L[[m]])>0){K_Lm_eff[m] <- length(active_L[[m]])}
    }
    
    ## 5.2: Specific Components ------------------------------------------------

    for(m in 1:M){

      delta_Gm[[m]]  <- as.vector(1)
      if(K_Gm[m]>1){
        # efficient computation of multivariate Normal & Student-T (log)-pdf
        vec_G  <- colSums(Gamma_m[[m]][,1:K_Gm[m],drop=F]^2)
        logP_spike <- -0.5*vec_G/var_spike[m] - rep(0.5*p_m[m]*log(2*pi*var_spike[m]),K_Gm[m])
        logP_slab <- -(0.5*p_m[m]+a_chi[m])*log(1+0.5*vec_G/b_chi[m]) +
          rep(lgamma(0.5*p_m[m]+a_chi[m]) - lgamma(a_chi[m]) - 0.5*p_m[m]*log(2*pi*b_chi[m]),K_Gm[m])

        # Likelihood Tempering
        logP_spike <- logP_spike * inv_pow_Gm[m]
        logP_slab <- logP_slab * inv_pow_Gm[m]
        
        delta_Gm[[m]] <- update_cusp(K_Gm[m],logP_spike,logP_slab,w_m[[m]])
      }
    }

    # identifying active columns
    for(m in 1:M){
      K_Gm_eff[m]   <- 0
      active_G[[m]] <- which( delta_Gm[[m]] > c(1:K_Gm[m]) ) # active columns among retained ones
      if(length(active_G[[m]])>0){K_Gm_eff[m] <- length(active_G[[m]])}
    }
    
    ## 5.3: Response Component -------------------------------------------------
    
    if(is_supervised){
      
      #| Remark: a priori: $ P[\zeta_{m h} = 0] = \xi_m $
       
      # Spike-Slab Prior for Theta
      logP_spike <- - 0.5*Theta_all^2/var_spike_y - rep(0.5*log(2*pi*var_spike_y),K+sum(K_Gm))
      logP_slab  <- - (0.5+a_theta)*log(1+0.5*Theta_all^2/b_theta) +
        rep(lgamma(0.5+a_theta) - lgamma(a_theta) - 0.5*log(2*pi*b_theta),K+sum(K_Gm))
      
      # Prior probabilities
      logP_spike <- logP_spike + log(c(rep(xi,K),rep(xi_m,K_Gm)))
      logP_slab  <- logP_slab  + log(1-c(rep(xi,K),rep(xi_m,K_Gm)))
      
      # Normalized probabilities
      lopP_max <- pmax(logP_spike,logP_slab)
      pr_D <- exp(logP_slab-lopP_max) / (exp(logP_spike-lopP_max)+exp(logP_slab-lopP_max))
      
      # Sampling 
      zeta_all  <- rbinom(K+sum(K_Gm),1,pr_D)
      
      zeta   <- zeta_all[1:K]
      zeta_m <- lapply(1:M, function(m) zeta_all[K+c(0,cumsum(K_Gm)[-M])[m]+1:K_Gm[m]])
      
      # Identifying active components
      active_T <- which(zeta==1)
      K_T_eff  <- length(active_T)
      for(m in 1:M){
        active_Tm[[m]] <- which(zeta_m[[m]]==1)
        K_Tm_eff[m]    <- length(active_Tm[[m]])
      }
    }

    # 6: Stick Breaking Elements -----------------------------------------------

    ## 6.1: Shared Component 
    for(m in 1:M){
      count_eq <- colSums(as.matrix(delta_Lm[m,] == matrix(1:K,K,K,byrow=T)))
      count_gr <- rev(cumsum(rev(c(count_eq[-1],0))))

      rho_m[m,] <- c(rbeta(K-1, shape1 = 1+count_eq[-K], shape2 = alpha_L[m]+count_gr[-K]),1.)
      pi_m[m,]  <- rho_m[m,]*c(1,cumprod(1-rho_m[m,-K]))
    }

    ## 6.2: Specific Components
    for(m in 1:M){
      count_eq <- colSums(as.matrix(delta_Gm[[m]] == t(c(1:K_Gm[m])*matrix(1,K_Gm[m],K_Gm[m]))))
      count_gr <- rev(cumsum(rev(c(count_eq[-1],0))))

      nu_m[[m]] <- c(rbeta(K_Gm[m]-1, shape1 = 1+count_eq[-K_Gm[m]], shape2 = alpha_G[m]+count_gr[-K_Gm[m]]),1.)
      w_m[[m]]  <- nu_m[[m]]*c(1,cumprod(1-nu_m[[m]][-K_Gm[m]]))
    }
    
    ## 6.3: Response
    if(is_supervised){
      xi <- rbeta(1, shape1 = a_xi+K-sum(zeta), shape2 = b_xi+sum(zeta))
      for(m in 1:M){
        xi_m[m] <- rbeta(1, shape1 = a_xi+K_Gm[m]-sum(zeta_m[[m]]),
                         shape2 = b_xi+sum(zeta_m[[m]]))
      }
    }

    # 7: Precision of Loading's Prior ------------------------------------------

    for(m in 1:M){
      chi_m[m,] <- rep(1./var_spike[m],K)
      if(length(active_L[[m]])>0){
        chi_m[m,][active_L[[m]]] <- rgamma(length(active_L[[m]]),shape=a_chi[m]+0.5*p_m[m],rate=1) *
        1./(b_chi[m] + 0.5 * colSums(Lambda_m[[m]][,active_L[[m]],drop=F]^2))
      }
    }

    for(m in 1:M){
      tau_m[[m]] <- rep(1./var_spike[m],K_Gm[m])
      if(length(active_G[[m]])>0){
        tau_m[[m]][active_G[[m]]] <- rgamma(length(active_G[[m]]),shape=a_chi[m]+0.5*p_m[m],rate=1) *
          1./(b_chi[m] + 0.5 * colSums(Gamma_m[[m]][,active_G[[m]],drop=F]^2))
      }
    }
    
    if(is_supervised){
      psi   <- rep(1./var_spike_y,K)
      psi_m <- lapply(1:M, function(m) rep(1./var_spike_y,K_Gm[m]))

      psi0 <- rgamma(1,shape=a_theta+0.5*sum(zeta_all),rate=1) *
        1./(b_theta + 0.5*sum(Theta_all^2*zeta_all))

      psi[active_T] <- psi0
      for(m in 1:M){
        psi_m[[m]][active_Tm[[m]]] <- psi0
      }
    }

    # 8: Adaptation ------------------------------------------------------------

    if((t>t0_adapt) & (t<tBurnIn) & (runif(1) < exp(t0 + t1*t))){

      ## 8.a K -----------------------------------------------------------------

      active_J <- which(colSums(delta_Lm > matrix(1:K,M,K,byrow=T)) > 1)
      K_eff <- length(active_J)
      
      if(K_eff < K-1){
        
        # TODO: check
        if(K_eff==0){
          K_eff = 1
          active_J = c(1)
        }

        K <- K_eff + 1

        eta <- cbind(eta[,active_J,drop=F],rnorm(n))

        for(m in 1:M){
          Lambda_m[[m]] <- cbind(Lambda_m[[m]][,active_J,drop=F],sqrt(var_spike[m])*rnorm(p_m[m]))
        }
        delta_Lm <- cbind(delta_Lm[,active_J,drop=F],rep(1,M))
        rho_m    <- cbind(rho_m[,active_J,drop=F],rep(1,M))
        pi_m     <- cbind(pi_m[,active_J,drop=F],rep(1,M)-rowSums(pi_m[,active_J,drop=F]))
        chi_m    <- cbind(chi_m[,active_J,drop=F],1./var_spike)
        
        if(is_supervised){
          Theta <- c(Theta[active_J],sqrt(var_spike_y)*rnorm(1))
          zeta <- c(zeta[active_J],0)
          psi   <- c(psi[active_J],1./var_spike_y)
        }
      } else if (K < K0){

        K <- K + 1

        eta <- cbind(eta,rnorm(n))

        for(m in 1:M){
          Lambda_m[[m]] <- cbind(Lambda_m[[m]],sqrt(var_spike[m])*rnorm(p_m[m]))
        }
        delta_Lm <- cbind(delta_Lm,rep(1,M))
        rho_m    <- cbind(rho_m[,-(K-1)],rbeta(1,shape1=1,shape2=alpha_L),rep(1,M))
        pi_m     <- rho_m * cbind(rep(1,M),t(apply(1-rho_m[,-K],1,cumprod)))
        chi_m    <- cbind(chi_m,1./var_spike)
        
        if(is_supervised){
          Theta <- c(Theta,sqrt(var_spike_y)*rnorm(1))
          zeta <- c(zeta,0)
          psi   <- c(psi,1./var_spike_y)
        }
      }

      # 8.b K_Gm ---------------------------------------------------------------

      for(m in 1:M){
        
        if(K_Gm_eff[m]==0){

          K_Gm[m] <- 1

          phi_m[[m]]   <- matrix(rnorm(n),n,1)

          Gamma_m[[m]] <- matrix(sqrt(var_spike[m])*rnorm(p_m[m]),p_m[m],1)
          w_m[[m]]     <- c(1.)
          tau_m[[m]]   <- c(1./var_spike[m])
          
          if(is_supervised){
            Theta_m[[m]] <- c(sqrt(var_spike_y)*rnorm(1))
            zeta_m[[m]]  <- c(0)
            psi_m[[m]]   <- c(1./var_spike_y)
          }
        } else if (K_Gm_eff[m] < K_Gm[m]-1) {

          K_Gm[m] <- K_Gm_eff[m] + 1

          phi_m[[m]]   <- cbind(phi_m[[m]][,active_G[[m]],drop=F],rnorm(n))

          Gamma_m[[m]] <- cbind(Gamma_m[[m]][,active_G[[m]],drop=F],sqrt(var_spike[m])*rnorm(p_m[m]))
          w_m[[m]]     <- c(w_m[[m]][active_G[[m]]],1-sum(w_m[[m]][active_G[[m]]]))
          tau_m[[m]]   <- c(tau_m[[m]][active_G[[m]]],1./var_spike[m])
          
          if(is_supervised){
            Theta_m[[m]] <- c(Theta_m[[m]][active_G[[m]]],sqrt(var_spike_y)*rnorm(1))
            zeta_m[[m]]  <- c(zeta_m[[m]][active_G[[m]]],0)
            psi_m[[m]]   <- c(psi_m[[m]][active_G[[m]]],1./var_spike_y)
          }
        } else if (K_Gm[m] < K0_m[m]) {

          K_Gm[m] <- K_Gm[m] + 1

          phi_m[[m]]   <- cbind(phi_m[[m]],rnorm(n))

          Gamma_m[[m]] <- cbind(Gamma_m[[m]],sqrt(var_spike[m])*rnorm(p_m[m]))
          nu_m[[m]]    <- c(nu_m[[m]][-(K_Gm[m]-1)],rbeta(1,shape1=1,shape2=alpha_G[m]),1.)
          w_m[[m]]     <- nu_m[[m]]*c(1,cumprod(1-nu_m[[m]][-K_Gm[m]]))
          tau_m[[m]]   <- c(tau_m[[m]],1./var_spike[m])
          
          if(is_supervised){
            Theta_m[[m]] <- c(Theta_m[[m]],sqrt(var_spike_y)*rnorm(1))
            zeta_m[[m]]  <- c(zeta_m[[m]],0)
            psi_m[[m]]   <- c(psi_m[[m]],1./var_spike_y)
          }
        }
      }
    }

    # 9. Saving outputs --------------------------------------------------------
    
    if((t %% tThin == 0)) {
      K_MC[t%/%tThin]         <- K # overall NUMBER of retained shared columns
      K_Gm_MC[t%/%tThin,]     <- K_Gm     # n. of retained columns in Gamma_m
    }
    
    if((t %% tThin == 0) & (t > tBurnIn)) {

      if(get_latent_vars){
        eta_MC[teff,,1:K]  <- eta
        for(m in c(1:M)){
          Lambda_m_MC[[m]][teff,,1:K] <- Lambda_m[[m]]
          Gamma_m_MC[[m]][teff,,1:K_Gm[m]]  <- Gamma_m[[m]]
          phi_m_MC[[m]][teff,,1:K_Gm[m]]    <- phi_m[[m]]
        }
      }

      K_Lm_eff_MC[teff,] <- K_Lm_eff # NUMBER of active columns in Lambda_m

      K_Gm_eff_MC[teff,] <- K_Gm_eff # NUMBER of active columns in Gamma_m

      for(m in c(1:M)){
        
        active_L_MC[teff,active_L[[m]],m] <- 1

        Marg_Var_m_MC[[m]][teff,] <- 1/s2_inv_m[[m]] + rowSums(Lambda_m[[m]]^2) + rowSums(Gamma_m[[m]]^2)
        s2_inv_m_MC[[m]][teff,]   <- s2_inv_m[[m]]
        mu_m_MC[[m]][teff,]       <- mu_m[[m]]

        Cov_m_mean[[m]] <- Cov_m_mean[[m]] +
          ( tcrossprod(Lambda_m[[m]]) + tcrossprod(Gamma_m[[m]]) + diag(1/s2_inv_m[[m]]) ) / tEff
      }

      if(impute_na){
        for(m in 1:M){
          for(idx in 1:length(na_idx[[m]])){
            Xm_MC[[m]][[idx]][teff,] <- X_m[[m]][na_row_idx[[m]][idx],unlist(na_idx[[m]][idx])]
          }
        }
      }

      if(is_supervised){
        
        Theta_MC[teff,1:K] <- Theta
        active_T_MC[teff,active_T] <- 1
        K_T_eff_MC[teff] <- K_T_eff
        
        for(m in 1:M){
          Theta_m_MC[[m]][teff,1:K_Gm[m]] <- Theta_m[[m]]
          active_Tm_MC[[m]][teff,active_Tm[[m]]] <- 1
          K_Tm_eff_MC[teff,m] <- K_Tm_eff[m]
        }
        
        s2_inv_MC[teff] <- s2_inv
        mu_y_MC[teff]   <- mu_y
        
        var_y_MC[teff]  <- 1/s2_inv + sum(Theta^2) + sum(sapply(1:M,function(m) sum(Theta_m[[m]]^2)))
        
        if(yBinary){y_MC[teff,] <- y}
      }

      teff = teff + 1
    }
    
    # ------ print status ------- #
    if(t %% iter_print == 0){
      cat(sprintf("%10s[%3d%%] K=%02d %s\n","", (t%/%iter_print)*10, K,
          paste(sprintf("K%d=%02d",1:M,K_Gm),collapse=" ")))
    }
  }
  
  # Trim Samples ---------------------------------------------------------------
  
  idx_BurnIn = c(1:(tBurnIn%/%tThin))
  
  Kmax = max(K_MC[-idx_BurnIn])
  Kmax_m = apply(K_Gm_MC[-idx_BurnIn,],2,max)
  
  active_L_MC <- active_L_MC[,1:Kmax,]

  if(get_latent_vars){
    eta_MC <- eta_MC[,,1:Kmax]
    for(m in c(1:M)){
      Lambda_m_MC[[m]] <- Lambda_m_MC[[m]][,,1:Kmax]
      Gamma_m_MC[[m]]  <- Gamma_m_MC[[m]][,,1:Kmax_m[m]]
      phi_m_MC[[m]]    <- phi_m_MC[[m]][,,1:Kmax_m[m]]
    }
  }
  
  if(is_supervised){
    Theta_MC<- Theta_MC[,1:Kmax]
    active_T_MC <- active_T_MC[,1:Kmax]
    for(m in 1:M){
      Theta_m_MC[[m]] <- Theta_m_MC[[m]][,1:Kmax_m[m]]
      active_Tm_MC[[m]] <- active_Tm_MC[[m]][,1:Kmax_m[m]]
    }
  }
  
  # Output ---------------------------------------------------------------------
  
  output = list(K=K_MC,K_Lm_eff=K_Lm_eff_MC, K_Gm=K_Gm_MC,K_Gm_eff=K_Gm_eff_MC,
                Cov_m_mean=Cov_m_mean,s2_inv_m=s2_inv_m_MC,mu_m=mu_m_MC,
                Marg_Var_m=Marg_Var_m_MC,active_Lm=active_L_MC)
  
  extra_params = list(model=paste0('jafar_',which_prior),
                      rescale_pred=rescale_pred, K0=K0, K0_m=K0_m, 
                      tBurnIn=tBurnIn, tMCMC=tMCMC, tThin=tThin,
                      tempered=tempered)
  hyperparams = c(extra_params, hyperparams)
  
  output$hyper_param=hyperparams
  
  if(is_supervised){
    output_y <- list(Theta=Theta_MC,Theta_m=Theta_m_MC,s2_inv=s2_inv_MC,
                     mu_y=mu_y_MC,active_T=active_T_MC,
                     active_Tm=active_Tm_MC,K_T_eff=K_T_eff_MC,
                     K_Tm_eff=K_Tm_eff_MC)
    if(yBinary){output_y$y_MC=y_MC}
    output = c(output, output_y)
  }
  
  if(impute_na){
    output_na = list(Xm_MC=Xm_MC,na_idx=na_idx,na_row_idx=na_row_idx)
    output = c(output, output_na)
  }
  
  if(get_latent_vars){
    output_latent = list(eta=eta_MC,Lambda_m=Lambda_m_MC,Gamma_m=Gamma_m_MC,phi_m=phi_m_MC)
    output = c(output, output_latent)
  }
  
  if(get_last_sample){
    output$last_sample = list(eta=eta,Lambda_m=Lambda_m,Gamma_m=Gamma_m,
                              phi_m=phi_m,s2_inv_m=s2_inv_m,mu_m=mu_m,
                              delta_Gm=delta_Gm,nu_m=nu_m,w_m=w_m,tau_m=tau_m,
                              delta_Lm=delta_Lm,rho_m=rho_m,pi_m=pi_m,chi_m=chi_m)
    if(is_supervised){
      output_last_y <- list(Theta=Theta,Theta_m=Theta_m,s2_inv=s2_inv,mu_y=mu_y,
                            zeta=zeta,zeta_m=zeta_m,psi=psi,psi_m=psi_m)
      output$last_sample = c(output$last_sample,output_last_y)
    }
  }
  
  return(output)
}












