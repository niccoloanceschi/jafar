
library(reticulate)
source_python("/Users/nico/Documents/GitHub/JAFAR_all/inspecting_indCUSP/get_array_indicators.py")
# Rcpp::sourceCpp('/Users/nico/Documents/GitHub/JAFAR_all/inspecting_indCUSP/get_array_indicators.cpp')

make_tensor_W_vectorized <- function(...) {
  vectors <- list(...)
  out <- Reduce(function(x, y) outer(x, y, "+"), vectors)
  return(out)
}

make_tensor_W <- function(log_prob_mat){
  do.call(make_tensor_W_vectorized, asplit(log_prob_mat, MARGIN = 1))
}

get_perm <- function(m,M){
  vec = 1:M
  vec[1] = m
  vec[m] = 1
  return(vec)
}

is.scalar <- function(x){ is.atomic(x) && length(x) == 1L }

jafar_set_hyperparameters <- function(hyperparams_list, M, is_supervised=F) {
  
  # Define default values for each hyperparameter
  mcmc_defaults <- list(
    seed = 123, # random seed for reproducibility
    t0 = -1, # adaptation log-prob intercept
    t1 = -5e-4, # adaptation log-prob slope
    t0_adapt = 20 # adaptation start
  )
  
  defaults <- list(
    a_m = 3, # shape of inv.gamma prior of idiosyncratic noise
    b_m = 1, # rate of inv.gamma prior of idiosyncratic noise
    prec0m = 1/0.5, # precision of normal prior on intercepts
    var_spike = 0.005, # variance of normal spike in cusps
    a_chi = 0.5, # shape of inv.gamma slab hyperprior in cusps
    b_chi = 0.1, # rate of inv.gamma slab hyperprior in cusps
    alpha = 5, # DP concentration parameter in shared component
    alpha_loc = 5  # DP concentration parameter in specific components
  )
  
  full_names = c(names(mcmc_defaults),names(defaults))
  
  if(is_supervised){
    resp_defaults <- list(
      a_sig = 3, # shape of inv.gamma prior of idiosyncratic noise
      b_sig = 1, # rate of inv.gamma prior of idiosyncratic noise
      prec0 = 1/0.5, # precision of normal prior on intercepts
      var_spike_y = 0.005, # variance of normal spike
      a_theta = 0.5, # shape of inv.gamma slab hyperprior
      b_theta = 1, # rate of inv.gamma slab hyperprior
      a_xi = 0.5, # shape1 of beta prior on mixture weight in response loadings variances
      b_xi = 0.5 # shape2 of beta prior on mixture weight in response loadings variances
    )
    full_names <- c(full_names, names(resp_defaults))
  }
  
  # Loop through defaults and assign if not present in the input list
  for (param_name in full_names) {
    if (is.null(hyperparams_list[[param_name]])) {
      hyperparams_list[[param_name]] <- defaults[[param_name]]
    }
    
    # Ensure vectorization if it's a scalar and needs to be repeated M times
    if (is.scalar(hyperparams_list[[param_name]]) && length(hyperparams_list[[param_name]]) == 1) {
      if (param_name %in% names(defaults)) {
        hyperparams_list[[param_name]] <- rep(hyperparams_list[[param_name]], M)
      }
    }
  }
  return(hyperparams_list)
}

jafar_initialize_output <- function(nEff, n, K0, K0_m, M, p_m,
                              get_latent_vars=F, is_supervised=F, yBinary=F) {
  
  vars_list <- list(
    K_MC        = array(NA,nEff), # overall NUMBER OF retained shared elements
    K_Lm_eff_MC = array(NA,c(nEff,M)), # NUMBER OF active columns in Lambda_m
    
    K_Gm_MC     = array(NA,c(nEff,M)), # n. of retained columns in Gamma_m
    K_Gm_eff_MC = array(NA,c(nEff,M)), # NUMBER OF active columns in Gamma_m
    
    s2_inv_m_MC = lapply(1:M, function(m) array(NA,c(nEff,p_m[m]))),
    mu_m_MC = lapply(1:M, function(m) array(0.,c(nEff,p_m[m]))),
    Marg_Var_m_MC = lapply(1:M, function(m) array(NA,c(nEff,p_m[m]))),
    Cov_m_mean = lapply(1:M, function(m) matrix(0,p_m[m],p_m[m])),
    
    active_L_MC = array(0,c(nEff,K0,M))
  )
  
  if(get_latent_vars){
    latent_vars_list <- list(
      Lambda_m_MC = lapply(1:M, function(m) array(0.,c(nEff,p_m[m],K0))),
      Gamma_m_MC = lapply(1:M, function(m) array(0.,c(nEff,p_m[m],K0_m[m]))),
      phi_m_MC = lapply(1:M, function(m) array(0.,c(nEff,n,K0_m[m]))),
      eta_MC = array(0.,c(nEff,n,K0))
    )
    vars_list <- c(vars_list, latent_vars_list)
  }
    
  if(is_supervised){
    resp_vars_list <- within(list(), {
      Theta_MC  <- array(0.,c(nEff,K0))
      Theta_m_MC <- lapply(1:M, function(m) array(0.,c(nEff,K0_m[m])))
      
      s2_inv_MC <- array(NA,c(nEff))
      mu_y_MC   <- array(NA,c(nEff))
      
      var_y_MC  <- array(NA,c(nEff))
      
      active_T_MC = array(0,c(nEff,K0))
      active_Tm_MC = lapply(1:M, function(m) array(0.,c(nEff,K0_m[m])))
      
      K_T_eff_MC  <- array(NA,nEff)
      K_Tm_eff_MC  <- array(NA,c(nEff,M))
      
      if(binary_y){ y_MC <- matrix(NA,nEff,n) }
    }) 
    vars_list <- c(vars_list, resp_vars_list)
  }
  
  return(vars_list)
}

jafar_initialize_sampler <- function(n, M, p_m, K0, K0_m,                   
                               hyperparams_list, is_supervised=F, yBinary=F) {
  
  param_init <- within(list(), {
    K <- K0
    K_Lm_eff <- rep(K0,M) # NUMBER of active columns in Lambda_m
    
    K_Gm     <- K0_m      # n. of retained columns in Gamma_m
    K_Gm_eff <- K0_m      # NUMBER of active columns in Gamma_m
    
    active_L <- active_G <- list() # active columns among retained ones
    
    eta     <- matrix(rnorm(n*K0),n,K0)
    
    Lambda_m <- Gamma_m <- phi_m <- s2_inv_m <- mu_m <- list()
    delta_Gm <- nu_m <- w_m <- tau_m <- list()
    
    delta_Lm <- matrix(K0,M,K0) 
    
    rho_m   <- matrix(NA,M,K0)
    pi_m    <- matrix(1/K0,M,K0) 
    chi_m   <- matrix(1,M,K0)
    
    for(m in c(1:M)){
      active_L[[m]] <- c(1:K0)
      active_G[[m]] <- c(1:K0_m[m])
      
      Lambda_m[[m]] <- matrix(NA,p_m[m],K0)
      Gamma_m[[m]]  <- matrix(NA,p_m[m],K0_m[m])
      phi_m[[m]]    <- matrix(rnorm(K0_m[m]*n),n,K0_m[m])
      
      s2_inv_m[[m]] <- rep(0.5,p_m[m])
      mu_m[[m]]     <- rep(0,p_m[m])
      
      delta_Gm[[m]] <- rep(K0_m[m],K0_m[m])
      
      nu_m[[m]]  <- rep(NA,K0_m[m])
      w_m[[m]]   <- rep(1/K0_m[m],K0_m[m])
      tau_m[[m]] <- rep(1,K0_m[m])
    }  
  })
  
  param_init <- param_init[!names(param_init) %in% c('m')]
  
  if(is_supervised){
    resp_param_init <- within(list(), {
      Theta   <- rep(0,K0)
      Theta_m <- lapply(1:M, function(m) rep(0,K0_m[m]))
      
      s2_inv  <- 0.5
      mu_y    <- 0.
      
      xi      <- 0.5
      xi_m    <- rep(0.5,M)
      
      psi     <- rep(1,K0)
      psi_m   <- lapply(1:M, function(m) rep(1,K0_m[m]))
      
      zeta    <- rep(1,K0)
      zeta_m  <- lapply(1:M, function(m) rep(1,K0_m[m]))
      
      active_T  <- rep(0,K0) 
      active_Tm <- lapply(1:M, function(m) rep(0,K0_m[m]))
      
      K_T_eff  <- 0
      K_Tm_eff <- rep(0,M)
    })
    param_init = c(param_init, resp_param_init)
  }
    
  return(param_init)
}

