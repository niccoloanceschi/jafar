
#' Initialize all unknown/latent random variables in the Gibbs sampler for JAFAR under the D-CUSP prior
#'
#' @param n number of observations
#' @param M number of views
#' @param p_m views dimensions (length M)
#' @param K0 Initial number of latent factors in shared components
#' @param K0_m Initial number of latent factors in view-specific components (length M)
#' @param a_sig Shape parameter of inverse-gamma prior on response noise
#' @param b_sig Scale parameter of inverse-gamma prior on response noise
#' @param a_theta Shape parameter in slab element of prior for response loadings
#' @param b_theta Scale parameter in slab element of prior for response loadings
#' @param var_spike_theta Variance parameter in spike element of prior for response-loadings
#' @param a_xi Shape1 parameters in beta of prior on response-loadings spike and slab weights
#' @param b_xi Shape2 parameters in beta of prior on response-loadings spike and slab weights
#' @param a_m Shape parameters of inverse-gamma prior on predictors idiosyncratic components
#' @param b_m Scale parameters of inverse-gamma prior on predictors idiosyncratic components
#' @param prec0 Prior precision for response intercept
#' @param prec0m Prior precisions for predictors intercepts
#' @param var_spike Variance parameter in spike element of prior for predictors-loadings
#' @param a_chi Shape parameters in slab element of prior fo shared-component loadings
#' @param b_chi Scale parameters in slab element of prior for shared-component loadings
#' @param alpha Stick-breaking parameter in shared-component loadings
#' @param alpha_loc Stick-breaking parameter in view-specific loadings
#' @param seed See for random number generation
#' @return A list of initial values for all unknown/latent random variables in the Gibbs sampler for JAFAR under the D-CUSP prior
#'
gibbs_JAFAR_CUSP_init <- function(n, M, p_m,          
                                  K0, K0_m,                   
                                  a_sig, b_sig,               
                                  a_theta, b_theta,           
                                  var_spike_theta,            
                                  a_xi, b_xi,                 
                                  a_m, b_m,                  
                                  prec0, prec0m,              
                                  var_spike,                  
                                  a_chi, b_chi,               
                                  alpha, alpha_loc,       
                                  seed) {
  
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

