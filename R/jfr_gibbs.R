
#' Gibbs sampler for \code{jfr}
#'
#' @description
#' Fits a Bayesian Joint Factor Regression (\code{jfr}) model using Gibbs sampling.
#' Variation across multiple data views is explained by a single set of global latent factors.
#' The model can be fitted in both unsupervised and supervised settings. 
#' Default and optional outputs include posterior means of the induced covariances,
#' posterior samples of residual variances, latent factors, and factor loadings. 
#' Supports parallel computation and tempered loading updates to limit rank estimation in extreme large-p-small-n settings.
#'
#' @references
#' Anceschi N., Ferrari F., Dunson D. B., & Mallick H. (2025). 
#' *Bayesian Joint Additive Factor Models for Multiview Learning*.
#' <https://arxiv.org/abs/2406.00778>
#' 
#' Legramanti S., Durante D., & Dunson D. B. (2020).
#' *Bayesian cumulative shrinkage for infinite factorizations*. 
#' Biometrika, 107(3), 745-752. <https://doi.org/10.1093/biomet/asaa008>
#'
#' @param X_m Multi-view input data, pre-processed via \code{\link{preprocess_X}}.
#'  List of length \code{M}; m-th element: matrix \code{n x p_m[m]}.
#'  Rows should correspond to samples, columns to features.
#' @param y Vector of responses (length \code{n}) pre-processed via \code{\link{preprocess_y}}.
#'  Set to \code{NULL} for unsupervised mode (default: \code{NULL}).
#' @param yBinary Logical, indicating if the response(s) are binary (default: \code{FALSE}).
#' @param K0 Upper bound to numbers of latent factors (optional)
#'  If \code{NULL}, \code{K0} is set to \code{3*log(max(p_m))}
#' @param tMCMC Total number of MCMC iterations (default: 20000).
#' @param tBurnIn Number of burn-in iterations (default: 15000).
#' @param tThin Thinning interval for saving samples (default: 10).
#' @param hyperparams List of hyperparameters for the \code{I-CUSP} prior distributions.
#'  Missing hyperparameters are replaced by default values encoded in \code{\link{set_hyperparameters}}.
#' @param get_latent_vars Return latent factors and loading matrices (logical, default: \code{TRUE}).
#' @param get_last_sample Return the last sample of the MCMC chain (logical, default: \code{FALSE}).
#' @param parallel Use parallel computation for the loadings update (logical, default: \code{TRUE}).
#' @param tempered Use tempered full-conditional for the loadings matrices (logical, default: \code{FALSE}).
#' @param rescale_pred Rescale loadings when computing response predictions (logical, default: \code{FALSE}).
#'
#' @details
#' - Ensure that all matrices in \code{X_m} have the same number of rows (subjects).  
#' - Missing data in \code{X_m} are allowed as \code{NA} and imputed in the MCMC.
#'
#' @return A list containing posterior samples, latent variables (if requested), and other relevant model outputs.
#'
#' @note
#'  All posterior samples are reported only after burn-in, except for \code{K}.
#'  The number of samples after thinning is \code{tFull=tMCMC%/%tThin} and \code{tEff=(tMCMC-tBurnIn)%/%tThin} 
#'  for the full chain and post burn-in, respectively. \cr
#'  The output list includes:
#'  \itemize{
#'    \item{\code{K}: Number shared latent factors (vector of length \code{tFull}).}
#'    \item{\code{K_Lm_eff}: Numbers of latent factors active in each view (matrix \code{tEff x M}).}
#'    \item{\code{active_Lm}: Binary indicators of latent factors activity across views (binary array \code{tEff x K x M}).}
#'    \item{\code{Cov_m_mean}: Posterior mean of the covariance matrix for each dataset (list of length \code{M}; m-th element: matrix \code{p_m[m] x p_m[m]}).}
#'    \item{\code{Marg_Var_m}: Marginal variances of features (list of length \code{M}; m-th element: matrix \code{tEff x p_m[m]}).}
#'    \item{\code{s2_inv_m}: Inverse residual variances across views (list of length \code{M}; m-th element: matrix \code{tEff x p_m[m]}).}
#'    \item{\code{mu_m}: Features intercepts across views (list of length \code{M}; m-th element: matrix \code{tEff x p_m[m]}).}
#'    \item{\code{hyper_param}: List of hyperparameters used for the model, including user-specified values and defaults ones were missing.}
#'  }
#'  If \code{is_supervised = TRUE}:
#'  \itemize{   
#'    \item{\code{K_T_eff}: Numbers of latent factors active in the response (vector of length \code{tEff}).}
#'    \item{\code{active_T}: Binary indicators of latent factors activity in the response (binary matrix \code{tEff x K}).}
#'    \item{\code{s2_inv}: Response inverse residual variances (vector of length \code{tEff}).}
#'    \item{\code{mu_y}: Response intercept (vector of length \code{tEff}).}
#'    \item{\code{Theta}: Response loadings on latent factors (matrix \code{tEff x K}).}
#'    \item{\code{y_MC}: Latent probit utilities (matrix \code{tEff x n}). (only if \code{yBinary = TRUE}).}
#'  }
#'  If \code{get_latent_vars = TRUE}:
#'  \itemize{   
#'    \item{\code{Lambda_m}: Loadings matrices on latent factors (list of length \code{M}; m-th element: array \code{tEff x p_m[m] x K}).}
#'    \item{\code{eta}: Latent factors (array \code{tEff x n x K}).}
#'  }
#'  If the input matrices \code{X_m} contain missing values:
#'  \itemize{
#'    \item{\code{Xm_MC}: Posterior samples of imputed values for missing entries.
#'      A list of length \code{M}; the m-th element is itself a list (one per feature with missingness),
#'      each containing an \code{tEff × n_miss} matrix of imputed values across MCMC iterations.}
#'    \item{\code{na_idx}: List of length \code{M}; the m-th element gives the column indices of missing entries in \code{X_m[[m]]}.}
#'    \item{\code{na_row_idx}: List of length \code{M}; the m-th element gives the corresponding row indices of missing entries in \code{X_m[[m]]}.}
#' }
#'  If \code{get_last_sample = TRUE}:
#'  \itemize{   
#'    \item{\code{last_sample}: List of posterior values of all parameters at the last MCMC iteration, including latent factors, loadings, residual variances, and hyperparameters.}
#'  }
#'
#' @export
#'
gibbs_jfr <- function(X_m, y=NULL, yBinary=F, K0=NULL, 
                      tMCMC=20000, tBurnIn=15000, tThin=10, 
                      hyperparams = list(), 
                      get_latent_vars=TRUE, get_last_sample=FALSE,
                      parallel=TRUE, tempered=FALSE, 
                      rescale_pred=FALSE){
  
  get_env <- environment()
  
  is_supervised <- !is.null(y)
  if(is_supervised){get_latent_vars=T}
  
  M = length(X_m) # number of modalities
  p_m = sapply(X_m,ncol) # number of features per modality
  n = nrow(X_m[[1]]) # number of samples
  
  # MCMC parameters
  tEff <- floor((tMCMC - tBurnIn)/tThin) 
  iter_print <- tMCMC %/% 10 
  teff <- 1
  
  # Latent dimensions
  if(is.null(K0)){K0<-floor(5*log(sum(p_m)))}
  
  # Set hyperparameters
  hyperparams <- set_hyperparameters(hyperparams, M, is_supervised)
  hyperparams <- hyperparams[!names(hyperparams) %in% c('alpha_loc')]
  list2env(hyperparams, envir=get_env)
  set.seed(seed)
  
  # Output Variables
  output_containers <- jafar_initialize_output(tEff, n, K0, rep(0,M), M, p_m,
                                               get_latent_vars, is_supervised, yBinary)
  output_containers <- output_containers[!names(output_containers) %in% c('K_Gm_MC', 'K_Gm_eff_MC')]
  if(get_latent_vars){
    output_containers <- output_containers[!names(output_containers) %in% c('Gamma_m_MC', 'phi_m_MC')]
  }
  if(is_supervised){
    names_to_remove <- c('K_Tm_eff_MC', 'active_Tm_MC', 'Theta_m_MC')
    output_containers <- output_containers[!names(output_containers) %in% names_to_remove]
  }
  list2env(output_containers, envir=get_env)
  rm(output_containers)
  K_MC <- rep(NA,tMCMC%/%tThin)
  
  # Initialization 
  par_init <- jafar_initialize_sampler(n, M, p_m, K0, rep(0,M), hyperparams,
                                       is_supervised, yBinary)
  names_to_remove <- c('delta_Gm', 'nu_m', 'w_m', 'tau_m', 'Gamma_m', 'phi_m','active_G', 'K_Gm_eff', 'K_Gm')
  par_init <- par_init[!names(par_init) %in% names_to_remove]
  if(is_supervised){
    names_to_remove <- c('K_Tm_eff', 'active_Tm', 'zeta_m', 'psi_m', 'Theta_m', 'xi_m')
    par_init <- par_init[!names(par_init) %in% names_to_remove]
  }
  list2env(par_init,envir=get_env)
  rm(par_init)
  
  # Auxiliary variables
  res_m <- eta_LaT <- list()
  
  # Likelihood Tempering Power
  inv_pow_Lm <- rep(1,M)
  if(tempered){
    inv_pow_Lm <- pmin(rep(n/M,M),p_m)/p_m
  }
  
  # Rescaling CUSP concentration
  alpha_L = alpha_L / sqrt(M)
  
  # ------ Latent Utilities & Missing Data -------------------------------------
  
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
  cat(sprintf("%10s[%3d%%] K=%02d\n", "", 0, K))
  
  for(t in c(1:tMCMC)){
    
    # 0: Latent Utilities & Missing Data Imputation ----------------------------
    
    if(impute_na & t>t0_adapt){
      for(m in 1:M){
        mar_std_m = rep(1,p_m[m])
        if(rescale_pred){
          mar_std_m = sqrt(1/s2_inv_m[[m]] + rowSums(Lambda_m[[m]]^2))
        }
        s2_m = s2_inv_m[[m]]*(mar_std_m^2)
        La_m = Lambda_m[[m]]/mar_std_m
        
        X_m_tmp <- matrix(mu_m[[m]],n,p_m[m],byrow=T) +
          t(matrix(rnorm(n*p_m[m]),p_m[m],n)/sqrt(s2_m)) +
          tcrossprod(eta,La_m)
        X_m[[m]][Xm_na[[m]]] <- X_m_tmp[Xm_na[[m]]]
      }
    }
    
    if(is_supervised & yBinary){
      linPred <- rep(mu_y,n) + c(eta%*%Theta)
      # y <- linPred + (2*y_obs-1)*truncnorm::rtruncnorm(n, a=-(2*y_obs-1)*linPred, b=rep(Inf,n))
      y <- truncnorm::rtruncnorm(n, a=left_thr, b=right_thr, mean=linPred, sd=1)
    }
    
    # 1: Loadings --------------------------------------------------------------
    
    etaTeta <- crossprod(eta)
    
    for(m in c(1:M)){
      
      if(parallel){
        new_Loadings <- update_loadings_parallel(n, p_m[m], K, X_m[[m]], etaTeta, 
                                                 eta, mu_m[[m]], s2_inv_m[[m]], chi_m[m,])
      } else {
        new_Loadings <- update_loadings(n, p_m[m], K, X_m[[m]], etaTeta, 
                                        eta, mu_m[[m]], s2_inv_m[[m]], chi_m[m,])
      }
      
      Lambda_m[[m]] <- new_Loadings[,c(1:K),drop=F]
      
      eta_LaT[[m]] <- tcrossprod(eta,Lambda_m[[m]])
    }
    
    if(is_supervised){
      fac_mu <- cbind(rep(1,n),eta)
      facTfac_mu <- crossprod(fac_mu)
      
      Q_Theta      <- s2_inv*facTfac_mu + diag(c(prec0,psi),1+K,1+K)
      r_Theta      <- s2_inv*crossprod(fac_mu,y)
      
      L_Theta      <- t(chol(Q_Theta))
      Lr_Theta     <- forwardsolve(L_Theta, r_Theta)
      
      mean_Theta   <- backsolve(t(L_Theta), Lr_Theta)
      std_Theta    <- backsolve(t(L_Theta), rnorm(1+K))
      
      mu_y      <- as.vector(mean_Theta[1] + std_Theta[1])[1]
      Theta     <- as.vector(mean_Theta[-1] + std_Theta[-1])
    }
    
    # 2: Intercepts ------------------------------------------------------------
    
    if(t>t0_adapt){
      for(m in 1:M){
        vec_m <- colMeans(X_m[[m]]) - colMeans(eta_LaT[[m]])
        mean_mu_m <- vec_m*s2_inv_m[[m]]/(s2_inv_m[[m]]+prec0m[m]/n)
        mu_m[[m]] <- mean_mu_m + sqrt(1/(prec0m[m]+n*s2_inv_m[[m]]))*rnorm(p_m[m])
      }
    }
    
    # 3: Idiosyncratic components ----------------------------------------------
    
    for(m in c(1:M)){
      res_m[[m]] <- X_m[[m]]-matrix(mu_m[[m]],n,p_m[m],byrow=T)
      s2_inv_m[[m]] <- rgamma(p_m[m],shape=a_m[m]+0.5*n,rate=1) * 
        1/(b_m[m]+0.5*colSums((res_m[[m]]-eta_LaT[[m]])^2))
    }
    
    if(is_supervised & !yBinary){
      res_y <- y - rep(mu_y,n) - as.vector(eta%*%Theta)
      s2_inv <- rgamma(1,shape=a_sig+0.5*n,rate=b_sig+0.5*sum(res_y^2))
    }
    
    # 4: Factors ---------------------------------------------------------------
    
    if(!is_supervised){
      eta <- update_factors_jfr(n, M, K, res_m, s2_inv_m, Lambda_m)
    } else {
      eta <- update_factors_jfr(n, M+1, K, c(res_m,list(matrix(y-mu_y,ncol=1))),
                                c(s2_inv_m, list(s2_inv)),
                                c(Lambda_m, list(matrix(Theta,nrow=1))))
    }
    
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
    
    for(m in 1:M){
      delta_Lm[m,] <- update_cusp(K,logP_Spikes[m,],logP_Slabs[m,],pi_m[m,])
    }
    
    # identifying active columns
    for(m in 1:M){
      K_Lm_eff[m]   <- 0
      active_L[[m]] <- which( delta_Lm[m,] > c(1:K) ) # active columns within retained ones
      if(length(active_L[[m]])>0){K_Lm_eff[m] <- length(active_L[[m]])}
    }
    
    ## 5.2: Response Component -------------------------------------------------
    
    if(is_supervised){
      
      #| Remark: a priori: $ P[\delta_h = 0] = \xi $
      
      # Spike-Slab Prior for Theta
      logP_spike <- - 0.5*Theta^2/var_spike_y - rep(0.5*log(2*pi*var_spike_y),K)
      logP_slab  <- - (0.5+a_theta)*log(1+0.5*Theta^2/b_theta) +
        rep(lgamma(0.5+a_theta) - lgamma(a_theta) - 0.5*log(2*pi*b_theta),K)
      
      # Prior probabilities
      logP_spike <- logP_spike + log(xi)
      logP_slab  <- logP_slab  + log(1-xi)
      
      # Normalized probabilities
      lopP_max <- pmax(logP_spike,logP_slab)
      pr_D <- exp(logP_slab-lopP_max) / (exp(logP_spike-lopP_max)+exp(logP_slab-lopP_max))
      
      # Sampling 
      zeta  <- rbinom(K,1,pr_D)
      
      # Identifying active components
      active_T <- which(zeta==1)
      K_T_eff  <- length(active_T)
    }
    
    # 6: Stick Breaking Elements -----------------------------------------------
    
    for(m in 1:M){
      count_eq <- colSums(as.matrix(delta_Lm[m,] == matrix(1:K,K,K,byrow=T)))
      count_gr <- rev(cumsum(rev(c(count_eq[-1],0))))
      
      rho_m[m,] <- c(rbeta(K-1, shape1 = 1+count_eq[-K], shape2 = alpha_L[m]+count_gr[-K]),1.)
      pi_m[m,]  <- rho_m[m,]*c(1,cumprod(1-rho_m[m,-K]))
    }
    
    if(is_supervised){
      xi <- rbeta(1, shape1 = a_xi+K-sum(zeta), shape2 = b_xi+sum(zeta))
    }
    
    # 7: Precision of Loading's Prior ------------------------------------------
    
    for(m in 1:M){
      chi_m[m,] <- rep(1./var_spike[m],K)
      if(length(active_L[[m]])>0){ 
        chi_m[m,][active_L[[m]]] <- rgamma(length(active_L[[m]]),shape=a_chi[m]+0.5*p_m[m],rate=1) *
          1./(b_chi[m] + 0.5 * colSums(Lambda_m[[m]][,active_L[[m]],drop=F]^2))
      }
    }
    
    if(is_supervised){
      psi   <- rep(1./var_spike_y,K)
      psi0 <- rgamma(1,shape=a_theta+0.5*sum(zeta),rate=b_theta + 0.5*sum(Theta^2*zeta))
      psi[active_T] <- psi0
    }
    
    # 8: Adaptation ------------------------------------------------------------
    
    if((t>t0_adapt) & (t<tBurnIn) & (runif(1) < exp(t0 + t1*t))){
      
      active_J <- unique(unlist(active_L))
      K_eff <- length(active_J)
      
      if(K_eff < K-1){
        
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
    }
    
    # 9. Saving outputs --------------------------------------------------------
    
    if((t %% tThin == 0)) {
      K_MC[t%/%tThin] <- K # overall NUMBER of retained shared columns
    }
    
    if((t %% tThin == 0) & (t > tBurnIn)) {
      
      if(get_latent_vars){
        eta_MC[teff,,1:K]  <- eta
        for(m in c(1:M)){
          Lambda_m_MC[[m]][teff,,1:K] <- Lambda_m[[m]]
        }
      }
      
      K_Lm_eff_MC[teff,] <- K_Lm_eff # NUMBER of active columns in Lambda_m
      
      for(m in c(1:M)){
        
        active_L_MC[teff,active_L[[m]],m] <- 1
        
        Marg_Var_m_MC[[m]][teff,] <- 1/s2_inv_m[[m]] + rowSums(Lambda_m[[m]]^2)
        s2_inv_m_MC[[m]][teff,]   <- s2_inv_m[[m]]
        mu_m_MC[[m]][teff,]       <- mu_m[[m]]
        
        Cov_m_mean[[m]] <- Cov_m_mean[[m]] +
          ( tcrossprod(Lambda_m[[m]]) + diag(1/s2_inv_m[[m]]) ) / tEff
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
        
        s2_inv_MC[teff] <- s2_inv
        mu_y_MC[teff]   <- mu_y
        
        var_y_MC[teff]  <- 1/s2_inv + sum(Theta^2) 
        
        if(yBinary){y_MC[teff,] <- y}
      }
      
      teff = teff + 1
    }
    
    # ------ print status ------- #
    if(t %% iter_print == 0){
      cat(sprintf("%10s[%3d%%] K=%02d\n", "", (t%/%iter_print)*10, K))
    }
  }
  
  # Trim Samples ---------------------------------------------------------------
  
  idx_BurnIn = c(1:(tBurnIn%/%tThin))
  
  Kmax = max(K_MC[-idx_BurnIn])  

  active_L_MC <- active_L_MC[,1:Kmax,]
  
  if(get_latent_vars){
    eta_MC <- eta_MC[,,1:Kmax]
    for(m in c(1:M)){
      Lambda_m_MC[[m]] <- Lambda_m_MC[[m]][,,1:Kmax]
    }
  }
  
  if(is_supervised){
    Theta_MC<- Theta_MC[,1:Kmax]
    active_T_MC <- active_T_MC[,1:Kmax]
  }
  
  # Output ---------------------------------------------------------------------
  
  output = list(K=K_MC,K_Lm_eff=K_Lm_eff_MC,
                Cov_m_mean=Cov_m_mean,s2_inv_m=s2_inv_m_MC,mu_m=mu_m_MC,
                Marg_Var_m=Marg_Var_m_MC,active_Lm=active_L_MC)
  
  extra_params = list(model='jfr_i-cusp',
                      rescale_pred=rescale_pred, K0=K0, 
                      tBurnIn=tBurnIn, tMCMC=tMCMC, tThin=tThin,
                      tempered=tempered)
  hyperparams = c(extra_params, hyperparams)
  
  output$hyper_param=hyperparams
  
  if(is_supervised){
    output_y <- list(Theta=Theta_MC,s2_inv=s2_inv_MC,mu_y=mu_y_MC,
                     active_T=active_T_MC,K_T_eff=K_T_eff_MC)
    if(yBinary){output_y$y_MC=y_MC}
    output = c(output, output_y)
  }
  
  if(impute_na){
    output_na = list(Xm_MC=Xm_MC,na_idx=na_idx,na_row_idx=na_row_idx)
    output = c(output, output_na)
  }
  
  if(get_latent_vars){
    output_latent = list(eta=eta_MC,Lambda_m=Lambda_m_MC)
    output = c(output, output_latent)
  }
  
  if(get_last_sample){
    output$last_sample = list(eta=eta,Lambda_m=Lambda_m,s2_inv_m=s2_inv_m,mu_m=mu_m,
                              delta_Lm=delta_Lm,rho_m=rho_m,pi_m=pi_m,chi_m=chi_m)
    if(is_supervised){
      output_last_y <- list(Theta=Theta,s2_inv=s2_inv,mu_y=mu_y,zeta=zeta,psi=psi)
      output$last_sample = c(output$last_sample,output_last_y)
    }
  }
  
  return(output)
}




