
is.scalar <- function(x){ is.atomic(x) && length(x) == 1L }

#' Set the hyperparameters for \code{jafar} and \code{jfr}
#'
#' @description Helper function to set hyperparameters for \code{\link{gibbs_jfr}} and \code{\link{gibbs_jafar}}.
#'  Supports both unsupervised and supervised settings.
#'    
#' @details Missing hyperparameters are assigned default values.
#'
#' @param hyperparams_list Named list of model hyperparameters.
#' @param M Integer, number of data-views.
#' @param is_supervised Running supervised model (logical, default: \code{FALSE}).
#'
#' @return A named list of hyperparameters with defaults filled in where missing.
#'    Scalar values are replicated \code{M} times where necessary.
#'
#' @note
#' Default hyperparameters include:
#' \itemize{
#'   \item{\code{seed}: random seed for reproducibility (default: 123).}
#'   \item{\code{t0, t1, t0_adapt}: adaptation parameters for MCMC (default: \code{t0=-1, t1=-5e-4, t0_adapt=200}).}
#'   \item{\code{a_m, b_m}: shape and rate of inverse-gamma prior for idiosyncratic noise in each view.}
#'      Scalars of vectors of length \code{M} (default: \code{a_m[m]=3, b_m[m]=1}).
#'   \item{\code{prec0m}: precision of normal prior on intercepts.}
#'      Scalar of vector of length \code{M} (default: \code{prec0m[m]=2}).
#'   \item{\code{var_spike}: variance of normal spike in cusps.}
#'      Scalar of vector of length \code{M} (default: \code{var_spike[m]=0.005}).
#'   \item{\code{a_chi, b_chi}: hyperparameters for slab inverse-gamma prior in cusps.
#'      Scalars of vectors of length \code{M} (default: \code{a_chi[m]=0.5, b_chi[m]=0.1}).}
#'   \item{\code{alpha_L, alpha_G}: Dirichlet process concentration parameters giving the expected number of factors, shared and local.
#'      Scalars of vectors of length \code{M} (default: \code{alpha_L[m]=5, alpha_G[m]=5}).}
#' }
#' If \code{is_supervised = TRUE}, additional hyperparameters for the response model are
#' \itemize{
#'   \item{\code{a_sig, b_sig}: shape and rate of inverse-gamma prior for idiosyncratic noise (default: \code{a_sig=3, b_sig=1}).}
#'   \item{\code{prec0}: precision of normal prior on intercept (default: \code{prec0=2}).}
#'   \item{\code{var_spike_y}: variance of normal spike (default: \code{var_spike_y=0.005}).}
#'   \item{\code{a_theta, b_theta}: hyperparameters for slab inverse-gamma prior in the slab (default: \code{a_theta=0.5, b_theta=0.1}).}
#'   \item{\code{a_xi, b_xi}: shape parameters for beta prior on mixture weight in response loadings (default: \code{a_xi=3, b_xi=2}).}
#' }
#'
#' @export
#' 
set_hyperparameters <- function(hyperparams_list, M, is_supervised=FALSE) {
  
  # Define default values for each hyperparameter
  mcmc_defaults <- list(
    seed = 123, 
    t0 = -0.5, 
    t1 = -5e-4, 
    t0_adapt = 20 
  )
  
  defaults <- list(
    a_m = 3, 
    b_m = 1, 
    prec0m = 1/0.5, 
    var_spike = 0.005, 
    a_chi = 0.5, 
    b_chi = 0.1, 
    alpha_L = 5, 
    alpha_G = 5
  )
  
  full_names = c(names(mcmc_defaults), names(defaults))
  full_values = c(mcmc_defaults, defaults)
  
  if(is_supervised){
    resp_defaults <- list(
      a_sig = 3, 
      b_sig = 1, 
      prec0 = 1/0.5, 
      var_spike_y = 0.005, 
      a_theta = 0.5, 
      b_theta = 1, 
      a_xi = 3, 
      b_xi = 2
    )
    full_names <- c(full_names, names(resp_defaults))
    full_values <- c(full_values, resp_defaults)
  }
  
  for (param_name in full_names) {
    if (is.null(hyperparams_list[[param_name]])) {
      hyperparams_list[[param_name]] <- full_values[[param_name]]
    }
    
    if (is.scalar(hyperparams_list[[param_name]]) && length(hyperparams_list[[param_name]]) == 1) {
      if (param_name %in% names(defaults)) {
        hyperparams_list[[param_name]] <- rep(hyperparams_list[[param_name]], M)
      }
    }
  }
  
  return
  
  return(hyperparams_list)
}


jafar_initialize_output <- function(tEff, n, K0, K0_m, M, p_m,
                              get_latent_vars=F, is_supervised=F, yBinary=F) {
  
  vars_list <- list(
    K_MC        = array(NA,tEff), # overall NUMBER OF retained shared elements
    K_Lm_eff_MC = array(NA,c(tEff,M)), # NUMBER OF active columns in Lambda_m
    
    K_Gm_MC     = array(NA,c(tEff,M)), # n. of retained columns in Gamma_m
    K_Gm_eff_MC = array(NA,c(tEff,M)), # NUMBER OF active columns in Gamma_m
    
    s2_inv_m_MC = lapply(1:M, function(m) array(NA,c(tEff,p_m[m]))),
    mu_m_MC = lapply(1:M, function(m) array(0.,c(tEff,p_m[m]))),
    Marg_Var_m_MC = lapply(1:M, function(m) array(NA,c(tEff,p_m[m]))),
    Cov_m_mean = lapply(1:M, function(m) matrix(0,p_m[m],p_m[m])),
    
    active_L_MC = array(0,c(tEff,K0,M))
  )
  
  if(get_latent_vars){
    latent_vars_list <- list(
      Lambda_m_MC = lapply(1:M, function(m) array(0.,c(tEff,p_m[m],K0))),
      Gamma_m_MC = lapply(1:M, function(m) array(0.,c(tEff,p_m[m],K0_m[m]))),
      phi_m_MC = lapply(1:M, function(m) array(0.,c(tEff,n,K0_m[m]))),
      eta_MC = array(0.,c(tEff,n,K0))
    )
    vars_list <- c(vars_list, latent_vars_list)
  }
    
  if(is_supervised){
    resp_vars_list <- within(list(), {
      Theta_MC  <- array(0.,c(tEff,K0))
      Theta_m_MC <- lapply(1:M, function(m) array(0.,c(tEff,K0_m[m])))
      
      s2_inv_MC <- array(NA,c(tEff))
      mu_y_MC   <- array(NA,c(tEff))
      
      var_y_MC  <- array(NA,c(tEff))
      
      active_T_MC = array(0,c(tEff,K0))
      active_Tm_MC = lapply(1:M, function(m) array(0.,c(tEff,K0_m[m])))
      
      K_T_eff_MC  <- array(NA,tEff)
      K_Tm_eff_MC  <- array(NA,c(tEff,M))
      
      if(yBinary){ y_MC <- matrix(NA,tEff,n) }
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

