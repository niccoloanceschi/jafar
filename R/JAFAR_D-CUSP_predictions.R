
#' Compute the induced regression coefficients in the linear predictor of y|X=x for a single MCMC iteration
#'
#' @param M number of views
#' @param K Number of latent factors in shared components
#' @param K Number of latent factors in view-specific components (length M)
#' @param p_m views dimensions (length M)
#' @param Theta Response loadings 
#' @param s2_inv_y Response noise precision
#' @param Theta Response loadings s2_inv_y
#' @param Lambda_m Shared-component loadings
#' @param Gamma_m View-specific loadings
#' @param s2_inv_m Precisions of idiosyncratic components in multi-view features
#' @param rescale_pred Rescale loadings to induce correlation matrix 
#' @return A list with the values of regression coefficients for each view in the linear predictor of y|X=x for a single MCMC iteration
#'
get_coeff_JAFAR <- function(M,K,K_m,p_m,Theta,s2_inv_y,Lambda_m,Gamma_m,s2_inv_m,rescale_pred=FALSE){
  beta_m <- list()
  C_inv <- diag(1,K,K)
  
  for(m in 1:M){
    
    mar_std_m = rep(1,p_m[m])
    if(rescale_pred){
      mar_std_m = sqrt(1/s2_inv_m[[m]] + rowSums(Lambda_m[[m]]^2) + rowSums(Gamma_m[[m]]^2))
    }
    Ga_m = Gamma_m[[m]]/mar_std_m
    s2_m = s2_inv_m[[m]]*(mar_std_m^2)
    La_m = Lambda_m[[m]]/mar_std_m
    
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

#' Compute the induced regression coefficients in the linear predictor of y|X=x for all MCMC iteration
#'
#' @param risMCMC Output of the Gibbs Sampler for JAFAR under the D-CUSP prior
#' @param rescale_pred Rescale loadings to induce correlation matrix 
#' @return A list with the values of regression coefficients for each view in the linear predictor of y|X=x for all MCMC iteration
#' 
#' @export
#' 
coeff_JAFAR <- function(risMCMC,rescale_pred=FALSE){
  
  M = length(risMCMC$mu_m)
  p_m = sapply(risMCMC$mu_m,ncol)
  
  tMCMC = dim(risMCMC$K_Gm)[1]
  iter_print <- tMCMC %/% 10 
  
  pred_coeff_mcmc <- sapply(1:M, function(m) matrix(NA,p_m[m],tMCMC))
  
  print(" - Computing Coefficients - ")
  print(sprintf(fmt = "%10s%3s%2s", "[",0,"%]"))
  
  for(t in 1:tMCMC){
    
    Ga_m_t <- lapply(1:M, function(m) risMCMC$Gamma_m[[m]][t,,1:risMCMC$K_Gm[t,m]])
    La_m_t <- lapply(risMCMC$Lambda_m, function(df) df[t,,1:risMCMC$K[t]])
    mu_m_t <- lapply(risMCMC$mu_m, function(df) df[t,])
    s2_m_t <- lapply(risMCMC$s2_inv_m, function(df) df[t,])
    
    Theta_t <- risMCMC$Theta[t,1:risMCMC$K[t]]
    
    coeff_t <- get_coeff_JAFAR(M,risMCMC$K[t],risMCMC$K_Gm[t,],p_m,
                              Theta_t,risMCMC$s2_inv[t],La_m_t,Ga_m_t,
                              s2_m_t,rescale_pred=rescale_pred)
    
    for(m in 1:M){
      pred_coeff_mcmc[[m]][,t] <- coeff_t$pred_coeff_m[[m]]
    }
    
    if(t %% iter_print == 0){
      print(sprintf(fmt = "%10s%3s%2s", "[",(t%/%iter_print)*10,"%]"))
    }
  }
  
  return(pred_coeff_mcmc)
  
}

#' Compute the linear predictor of y|X=x in out-of-sample observations for a single MCMC iteration via sampling of the latent factors
#'
#' @param Xpred Multi-view features in out-of-sample observations
#' @param nPred Number of out-of-sample observations
#' @param M number of views
#' @param p_m views dimensions (length M)
#' @param K Number of latent factors in shared components
#' @param K_m Number of latent factors in view-specific components (length M)
#' @param Theta Response loadings 
#' @param Lambda_m Shared-component loadings
#' @param Gamma_m View-specific loadings
#' @param mu_y Response Intercepts
#' @param s2_inv_y Response noise precision
#' @param mu_m Intercepts of idiosyncratic components in multi-view features
#' @param s2_inv_m Precisions of idiosyncratic components in multi-view features
#' @param rescale_pred Rescale loadings to induce correlation matrix 
#' @return The values of linear predictor of y|X=x in out-of-sample observations for a single MCMC iteration
#'
y_cond_pred_JAFAR_sampling <- function(Xpred,nPred,M,p_m,K,K_m,
                                       Theta,Lambda_m,Gamma_m,mu_y,s2_inv_y,mu_m,s2_inv_m,
                                       rescale_pred=FALSE){
  
  Q_eta <- diag(1,K,K)
  r_eta <- rep(0,K)
  
  for(m in 1:M){
    
    if(!is.null(Xpred[[m]])){
      
      mar_std_m = rep(1,p_m[m])
      if(rescale_pred){
        mar_std_m = sqrt(1/s2_inv_m[[m]] + rowSums(as.matrix(Lambda_m[[m]]^2)) + rowSums(as.matrix(Gamma_m[[m]]^2)))
      }
      Ga_m = Gamma_m[[m]]/mar_std_m
      s2_m = s2_inv_m[[m]]*(mar_std_m^2)
      La_m = Lambda_m[[m]]/mar_std_m
      
      s2_La_m <- s2_m*La_m
      s2_Ga_m <- s2_m*Ga_m
      GaT_s2_La_m = crossprod(s2_Ga_m,La_m)
      
      D_m_inv = chol(diag(1.,K_m[m],K_m[m])+crossprod(Ga_m,s2_Ga_m))
      D_GaT_s2_La_m = backsolve(D_m_inv,forwardsolve(t(D_m_inv),GaT_s2_La_m))
      
      res_m <- Xpred[[m]]-rep(1,nPred)%x%t(mu_m[[m]])
      
      Q_eta <- Q_eta + crossprod(La_m,s2_La_m) - crossprod(GaT_s2_La_m,D_GaT_s2_La_m)
      r_eta <- r_eta + t( res_m%*%(s2_La_m - s2_Ga_m%*%D_GaT_s2_La_m) )
    }
  }
  
  L_eta    <- t(chol(Q_eta))
  Lr_eta   <- forwardsolve(L_eta, r_eta)
  
  mean_eta <- backsolve(t(L_eta), Lr_eta)
  std_eta  <- backsolve(t(L_eta), matrix(rnorm(K*nPred),K,nPred))
  
  eta_mc      <- t(mean_eta + std_eta)
  
  lin_pred_y <- rep(mu_y,nPred)+eta_mc%*%Theta
  
  return(lin_pred_y)
}


#' Compute the linear predictor of y|X=x in out-of-sample observations for a single MCMC iteration via Exact Expression & with missing data in features
#'
#' @param Xpred Multi-view features in out-of-sample observations
#' @param nPred Number of out-of-sample observations
#' @param M number of views
#' @param p_m views dimensions (length M)
#' @param K Number of latent factors in shared components
#' @param K_m Number of latent factors in view-specific components (length M)
#' @param Theta Response loadings 
#' @param Lambda_m Shared-component loadings
#' @param Gamma_m View-specific loadings
#' @param mu_y Response Intercepts
#' @param s2_inv_y Response noise precision
#' @param mu_m Intercepts of idiosyncratic components in multi-view features
#' @param s2_inv_m Precisions of idiosyncratic components in multi-view features
#' @param na_row_idx view-wise list of observations with missing entries
#' @param na_idx view-wise list of missing feature in each observations
#' @param rescale_pred Rescale loadings to induce correlation matrix 
#' @return The values of linear predictor of y|X=x in out-of-sample observations for a single MCMC iteration
#'
y_cond_pred_JAFAR_NA <- function(Xpred,nPred,M,p_m,K,K_m,
                                 Theta,Lambda_m,Gamma_m,
                                 mu_y,s2_inv_y,mu_m,s2_inv_m,
                                 na_row_idx,na_idx,
                                 rescale_pred=FALSE){
  
  lin_pred_y <- rep(mu_y,nPred)
  
  Ga_m <- La_m <- s2_Ga_m <- s2_La_m <- C_inv_m <- list()
  GaT_s2_Ga_m <- GaT_s2_La_m <- LaT_s2_La_m <- D_GaT_s2_La_m <- list()
  
  for(m in 1:M){
    
    if(!is.null(Xpred[[m]])){
      
      mar_std_m = rep(1,p_m[m])
      if(rescale_pred){
        mar_std_m = sqrt(1/s2_inv_m[[m]] + rowSums(Lambda_m[[m]]^2) + rowSums(Gamma_m[[m]]^2))
      }
      s2_m      <- s2_inv_m[[m]]*(mar_std_m^2)
      Ga_m[[m]] <- Gamma_m[[m]]/mar_std_m
      La_m[[m]] <- Lambda_m[[m]]/mar_std_m
      
      s2_Ga_m[[m]] <- s2_m*Ga_m[[m]]
      s2_La_m[[m]] <- s2_m*La_m[[m]]
      
      GaT_s2_Ga_m[[m]] <- crossprod(s2_Ga_m[[m]],Ga_m[[m]])
      GaT_s2_La_m[[m]] <- crossprod(s2_Ga_m[[m]],La_m[[m]])
      LaT_s2_La_m[[m]] <- crossprod(s2_La_m[[m]],La_m[[m]])
      
      D_m_chol0 <- chol(diag(1.,K_m[m],K_m[m]) + GaT_s2_Ga_m[[m]])
      
      D_GaT_s2_La_m[[m]] <- backsolve(D_m_chol0,forwardsolve(t(D_m_chol0),GaT_s2_La_m[[m]]))
      
      C_inv_m[[m]] <- LaT_s2_La_m[[m]] - crossprod(GaT_s2_La_m[[m]],D_GaT_s2_La_m[[m]])
    }
  }
  
  for(i in 1:nPred){
    
    Ci_inv  <- diag(1,K,K)
    La_D_Xi <- rep(0.,K)
    
    for(m in 1:M){
      
      if(!is.null(Xpred[[m]])){
        
        idx_i <- match(i,na_row_idx[[m]])
        
        if(is.na(idx_i)){
          Ci_inv  <- Ci_inv + C_inv_m[[m]]
          
          # La_D_Xi <- La_D_Xi + colSums(s2_La_m[[m]]*(Xpred[[m]][i,]-mu_m[[m]])) -
          #   colSums(D_GaT_s2_La_m[[m]]*colSums(s2_Ga_m[[m]]*(Xpred[[m]][i,]-mu_m[[m]])))
          
          La_D_Xi <- La_D_Xi + crossprod(s2_La_m[[m]],Xpred[[m]][i,]-mu_m[[m]]) -
            crossprod(D_GaT_s2_La_m[[m]],crossprod(s2_Ga_m[[m]],Xpred[[m]][i,]-mu_m[[m]]))
          
        } else {
          
          idx_m_i <- unlist(na_idx[[m]][idx_i])
          
          D_m_inv <- chol( diag(1.,K_m[m],K_m[m]) + GaT_s2_Ga_m[[m]] -
                             crossprod(s2_Ga_m[[m]][idx_m_i,,drop=F],Ga_m[[m]][idx_m_i,,drop=F]) ) 
          
          GaT_s2_La_m_NA   <- GaT_s2_La_m[[m]] -
            crossprod(s2_Ga_m[[m]][idx_m_i,,drop=F],La_m[[m]][idx_m_i,,drop=F]) 
          
          D_GaT_s2_La_m_NA <- backsolve(D_m_inv,forwardsolve(t(D_m_inv),GaT_s2_La_m_NA))
          
          Ci_inv  <- Ci_inv + LaT_s2_La_m[[m]] -
            crossprod(s2_La_m[[m]][idx_m_i,,drop=F],La_m[[m]][idx_m_i,,drop=F]) -
            crossprod(D_GaT_s2_La_m_NA,GaT_s2_La_m_NA)
          
          La_D_Xi <- La_D_Xi + crossprod(s2_La_m[[m]][-idx_m_i,,drop=F],
                                         Xpred[[m]][i,-idx_m_i]-mu_m[[m]][-idx_m_i]) -
            crossprod(D_GaT_s2_La_m_NA,crossprod(s2_Ga_m[[m]][-idx_m_i,,drop=F],
                                                 Xpred[[m]][i,-idx_m_i]-mu_m[[m]][-idx_m_i]))
        }
      }
    }
    
    C_chol = chol(Ci_inv)
    Theta_C = backsolve(C_chol,forwardsolve(t(C_chol),Theta[1:K,drop=F]))
    
    lin_pred_y[i] <- lin_pred_y[i] + sum(Theta_C*La_D_Xi)
  }
  
  return(lin_pred_y)
}

#' Compute the linear predictor of y|X=x in out-of-sample observations for all MCMC iteration
#'
#' @param Xpred Multi-view features in out-of-sample observations
#' @param risMCMC Output of the Gibbs Sampler for JAFAR under the D-CUSP prior
#' @param rescale_pred Rescale loadings to induce correlation matrix 
#' @return Linear predictors of y|X=x in out-of-sample observations for all MCMC iteration
#' 
#' @export
#' 
y_pred_JAFAR <- function(Xpred,risMCMC,rescale_pred=FALSE){
  
  M = length(Xpred)
  nPred = unlist(sapply(Xpred,nrow))[1]
  p_m = sapply(Xpred,ncol)
  
  NA_in_X <- max(sapply(Xpred,function(df) max(is.na(df))))
  
  if(NA_in_X){
    get_NA_pred <- get_NA_X(Xpred)
    na_row_idx  <- get_NA_pred$na_row_idx
    na_idx      <- get_NA_pred$na_idx
  }
  
  tMCMC = dim(risMCMC$K_Gm)[1]
  iter_print <- tMCMC %/% 10 
  
  lin_pred_mcmc  <- matrix(NA,nPred,tMCMC)
  
  print(" - Computing Response Predictions - ")
  print(sprintf(fmt = "%10s%3s%2s", "[",0,"%]"))
  
  for(t in 1:tMCMC){
    
    Ga_m_t  <- lapply(1:M, function(m) matrix(risMCMC$Gamma_m[[m]][t,,1:risMCMC$K_Gm[t,m]],ncol=risMCMC$K_Gm[t,m]))
    La_m_t  <- lapply(risMCMC$Lambda_m, function(df) matrix(df[t,,1:risMCMC$K[t]],ncol=risMCMC$K[t]))
    mu_m_t  <- lapply(risMCMC$mu_m, function(df) df[t,])
    s2_m_t  <- lapply(risMCMC$s2_inv_m, function(df) df[t,])
    Theta_t <- risMCMC$Theta[t,1:risMCMC$K[t]]
    
    if(NA_in_X){
      lin_pred_mcmc[,t] <- y_cond_pred_JAFAR_NA(Xpred,nPred,M,p_m,
                                                risMCMC$K[t],risMCMC$K_Gm[t,],
                                                Theta_t,La_m_t,Ga_m_t,
                                                risMCMC$mu_y[t],risMCMC$s2_inv[t],
                                                mu_m_t,s2_m_t,
                                                rescale_pred=rescale_pred,
                                                na_row_idx,na_idx)
    } else {
      lin_pred_mcmc[,t] <- y_cond_pred_JAFAR_sampling(Xpred,nPred,M,p_m,
                                                      risMCMC$K[t],risMCMC$K_Gm[t,],
                                                      Theta_t,La_m_t,Ga_m_t,
                                                      risMCMC$mu_y[t],risMCMC$s2_inv[t],
                                                      mu_m_t,s2_m_t,
                                                      rescale_pred=rescale_pred)
    }
    
    if(t %% iter_print == 0){
      print(sprintf(fmt = "%10s%3s%2s", "[",(t%/%iter_print)*10,"%]"))
    }
  }
  
  return(lin_pred_mcmc)
  
}





















