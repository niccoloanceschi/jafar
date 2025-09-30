
# response coefficients ----

jfr_coeff_y_t <- function(Xpred,nPred,M,p_m,K,Lambda_m,Theta,mu_y,
                         s2_inv_y,mu_m,s2_inv_m,rescale_pred=F){
  
  lin_pred <- rep(mu_y,nPred)
  
  Q_eta <- diag(1,K,K)
  r_eta_m <- coeff_m <- list()
  
  for(m in 1:M){
    
    mar_std_m = rep(1,p_m[m])
    if(rescale_pred){
      mar_std_m = sqrt(1/s2_inv_m[[m]] + rowSums(as.matrix(Lambda_m[[m]]^2)))
    }
    s2_m = s2_inv_m[[m]]*(mar_std_m^2)
    La_m = Lambda_m[[m]]/mar_std_m
    
    s2_La_m <- s2_m*La_m
    
    Q_eta <- Q_eta + crossprod(La_m,s2_La_m) 
    r_eta_m[[m]] <- t(s2_La_m)
  }
  
  L_eta    <- t(chol(Q_eta))
  
  for(m in 1:M){
    coeff_m[[m]] <- as.vector(Theta%*%backsolve(t(L_eta),forwardsolve(L_eta, r_eta_m[[m]])))
    
    if(!is.null(Xpred[[m]])){
      res_m <- Xpred[[m]]-rep(1,nPred)%x%t(mu_m[[m]])
      res_m[is.na(res_m)] <- 0
      
      lin_pred <- lin_pred + as.vector(res_m%*%coeff_m[[m]])
    }
  }
  
  return(list(lin_pred=lin_pred,coeff=coeff_m))
}

jafar_coeff_y_t <- function(Xpred,nPred,M,p_m,K,K_m,Lambda_m,Gamma_m,Theta,
                            Theta_m,mu_y,s2_inv_y,mu_m,s2_inv_m,rescale_pred=F){
  
  lin_pred <- rep(mu_y,nPred)
  coeff_m <- lapply(1:M, function(m) rep(0,p_m[m]))
  
  Q_eta <- diag(1,K,K)
  r_eta_m <- list()
  
  for(m in 1:M){
    
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
    D_Theta_m = backsolve(D_m_inv,forwardsolve(t(D_m_inv),Theta_m[[m]]))
    
    coeff_m[[m]] <- coeff_m[[m]] + as.vector(s2_Ga_m%*%D_Theta_m)
    
    Theta <- Theta - as.vector(Theta_m[[m]]%*%D_GaT_s2_La_m)
    
    Q_eta <- Q_eta + crossprod(La_m,s2_La_m) - crossprod(GaT_s2_La_m,D_GaT_s2_La_m)
    
    r_eta_m[[m]] <- t(s2_La_m - s2_Ga_m%*%D_GaT_s2_La_m)
  }
  
  L_eta    <- t(chol(Q_eta))
  LtL_Theta <- backsolve(t(L_eta),forwardsolve(L_eta,Theta))
  
  for(m in 1:M){
    coeff_m[[m]] <- coeff_m[[m]] + as.vector(LtL_Theta%*%r_eta_m[[m]])
    
    if(!is.null(Xpred[[m]])){
      res_m <- Xpred[[m]]-rep(1,nPred)%x%t(mu_m[[m]])
      res_m[is.na(res_m)] <- 0
      
      lin_pred <- lin_pred + as.vector(res_m%*%coeff_m[[m]])
    }
  }

  return(list(lin_pred=lin_pred,coeff=coeff_m))
}

#' Response predictions for \code{jafar} and \code{jfr}
#' 
#' @description Compute induced regression coefficients and predicted responses either in-sample or out-of-sample.
#'
#' @param Xpred A list of \code{M} features matrices, the m-th one of dimension \code{nPred x p_m[m]} or possibly missing (i.e. \code{X_m[[m]]=NULL}).
#' @param risMCMC Output of \code{\link{gibbs_jafar}} or \code{\link{gibbs_jfr}} containing posterior samples.
#' @param rescale_pred Rescale loadings when computing response predictions (logical, default: \code{FALSE}).
#'
#' @return A list containing posterior samples of the predicted responses (matrix \code{tEff x nPred}),
#'    and of the induced regression coefficients for each view (list of length \code{M}; m-th element: \code{tEff x p_m[m]}).
#'
#' @export
#'
predict_y <- function(Xpred,risMCMC,rescale_pred=FALSE){
  
  M = length(risMCMC$mu_m)
  p_m = sapply(risMCMC$mu_m,ncol)
  
  nPred = unlist(sapply(Xpred,nrow))[1]
  
  tMCMC = nrow(risMCMC$K_Lm_eff)
  iter_print <- max(1,tMCMC %/% 10)  
  
  mu_MC  <- matrix(NA,tMCMC,nPred)
  coeff_MC <- lapply(1:M, function(m) matrix(NA,tMCMC,p_m[m]))
  
  t_vals = 1:(risMCMC$hyper_param$tBurnIn%/%risMCMC$hyper_param$tThin)
  
  risMCMC$K = risMCMC$K[-t_vals]
  risMCMC$K_Gm = risMCMC$K_Gm[-t_vals,]
  
  print(" - Computing Response Predictions - ")
  print(sprintf(fmt = "%10s%3s%2s", "[",0,"%]"))
  
  for(t in 1:tMCMC){
    
    La_m_t <- lapply(risMCMC$Lambda_m, function(df) df[t,,1:risMCMC$K[t]])
    Theta_t   <- risMCMC$Theta[t,1:risMCMC$K[t]]
    
    mu_m_t <- lapply(risMCMC$mu_m, function(df) df[t,])
    s2_m_t <- lapply(risMCMC$s2_inv_m, function(df) df[t,])
    
    if(grepl('jafar',risMCMC$hyper_param$model)){
      Ga_m_t <- lapply(1:M, function(m) risMCMC$Gamma_m[[m]][t,,1:risMCMC$K_Gm[t,m]])
      Theta_m_t <- lapply(1:M, function(m) risMCMC$Theta_m[[m]][t,1:risMCMC$K_Gm[t,m]])
      
      pred_t <- jafar_coeff_y_t(Xpred,nPred,M,p_m,risMCMC$K[t],risMCMC$K_Gm[t,],
                               La_m_t,Ga_m_t,Theta_t,Theta_m_t,
                               risMCMC$mu_y[t],risMCMC$s2_inv[t],mu_m_t,s2_m_t,
                               rescale_pred=rescale_pred)   
    } else {
      pred_t <- jfr_coeff_y_t(Xpred,nPred,M,p_m,risMCMC$K[t],
                             La_m_t,Theta_t,
                             risMCMC$mu_y[t],risMCMC$s2_inv[t],mu_m_t,s2_m_t,
                             rescale_pred=rescale_pred)   
    }
    
    mu_MC[t,]  <- pred_t$lin_pred
    for(m in 1:M){
      coeff_MC[[m]][t,]  <- pred_t$coeff[[m]]
    }
    
    if(t %% iter_print == 0){
      print(sprintf(fmt = "%10s%3s%2s", "[",(t%/%iter_print)*10,"%]"))
    }
  }
  
  return(list(mean=mu_MC,coeff=coeff_MC))
  
}


# response predictions ----

jfr_pred_y_t <- function(Xpred,nPred,M,p_m,K,Lambda_m,Theta,mu_y,
                         s2_inv_y,mu_m,s2_inv_m,rescale_pred=F){
  
  lin_pred <- rep(mu_y,nPred)
  
  Q_eta <- diag(1,K,K)
  r_eta <- matrix(0,K,nPred)
  
  for(m in 1:M){
    
    mar_std_m = rep(1,p_m[m])
    if(rescale_pred){
      mar_std_m = sqrt(1/s2_inv_m[[m]] + rowSums(as.matrix(Lambda_m[[m]]^2)))
    }
    s2_m = s2_inv_m[[m]]*(mar_std_m^2)
    La_m = Lambda_m[[m]]/mar_std_m
    
    s2_La_m <- s2_m*La_m
    
    Q_eta <- Q_eta + crossprod(La_m,s2_La_m) 
    
    if(!is.null(Xpred[[m]])){
      
      res_m <- Xpred[[m]]-rep(1,nPred)%x%t(mu_m[[m]])
      res_m[is.na(res_m)] <- 0
    
      r_eta <- r_eta + t(res_m%*%s2_La_m)
    }
  }
  
  L_eta    <- t(chol(Q_eta))
  Lr_eta   <- backsolve(t(L_eta),forwardsolve(L_eta, r_eta))
  
  lin_pred <- lin_pred + as.vector(crossprod(Lr_eta,Theta))
  
  return(lin_pred)
}

jafar_pred_y_t <- function(Xpred,nPred,M,p_m,K,K_m,Lambda_m,Gamma_m,Theta,
                           Theta_m,mu_y,s2_inv_y,mu_m,s2_inv_m,rescale_pred=F){
  
  lin_pred <- rep(mu_y,nPred)
  
  Q_eta <- diag(1,K,K)
  r_eta <- matrix(0,K,nPred)
  
  for(m in 1:M){
    
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
    D_Theta_m = backsolve(D_m_inv,forwardsolve(t(D_m_inv),Theta_m[[m]]))
    
    Theta <- Theta - as.vector(Theta_m[[m]]%*%D_GaT_s2_La_m)
    
    Q_eta <- Q_eta + crossprod(La_m,s2_La_m) - crossprod(GaT_s2_La_m,D_GaT_s2_La_m)
    
    if(!is.null(Xpred[[m]])){
      
      res_m <- Xpred[[m]]-rep(1,nPred)%x%t(mu_m[[m]])
      res_m[is.na(res_m)] <- 0
      
      lin_pred <- lin_pred + as.vector(res_m%*%(s2_Ga_m%*%D_Theta_m))
      
      r_eta <- r_eta + t( res_m%*%(s2_La_m - s2_Ga_m%*%D_GaT_s2_La_m) )
    }
  }
  
  L_eta    <- t(chol(Q_eta))
  Lr_eta   <- backsolve(t(L_eta),forwardsolve(L_eta, r_eta))
  
  lin_pred <- lin_pred + as.vector(crossprod(Lr_eta,Theta))
  
  return(lin_pred)
}

#' Response predictions for \code{jafar} and \code{jfr}
#' 
#' @description Compute predicted responses, either in-sample or out-of-sample.
#'
#' @param Xpred A list of \code{M} features matrices, the m-th one of dimension \code{nPred x p_m[m]} or possibly missing (i.e. \code{X_m[[m]]=NULL}).
#' @param risMCMC Output of \code{\link{gibbs_jafar}} or \code{\link{gibbs_jfr}} containing posterior samples.
#' @param rescale_pred Rescale loadings when computing response predictions (logical, default: \code{FALSE}).
#'
#' @return A list containing posterior samples of the predicted responses (matrix \code{tEff x nPred}).
#'
#' @export
#'
predict_y_raw <- function(Xpred,risMCMC,rescale_pred=FALSE){
  
  M = length(risMCMC$mu_m)
  p_m = sapply(risMCMC$mu_m,ncol)
  
  nPred = unlist(sapply(Xpred,nrow))[1]
  
  tMCMC = nrow(risMCMC$K_Lm_eff)
  iter_print <- max(1,tMCMC %/% 10)  
  
  mu_MC  <- matrix(NA,tMCMC,nPred)
  
  t_vals = 1:(risMCMC$hyper_param$tBurnIn%/%risMCMC$hyper_param$tThin)
  
  risMCMC$K = risMCMC$K[-t_vals]
  risMCMC$K_Gm = risMCMC$K_Gm[-t_vals,]
  
  print(" - Computing Response Predictions - ")
  print(sprintf(fmt = "%10s%3s%2s", "[",0,"%]"))
  
  for(t in 1:tMCMC){
    
    La_m_t  <- lapply(risMCMC$Lambda_m, function(df) df[t,,1:risMCMC$K[t]])
    Theta_t <- risMCMC$Theta[t,1:risMCMC$K[t]]
    
    mu_m_t <- lapply(risMCMC$mu_m, function(df) df[t,])
    s2_m_t <- lapply(risMCMC$s2_inv_m, function(df) df[t,])
    
    if(grepl('jafar',risMCMC$hyper_param$model)){
      Ga_m_t <- lapply(1:M, function(m) risMCMC$Gamma_m[[m]][t,,1:risMCMC$K_Gm[t,m]])
      Theta_m_t <- lapply(1:M, function(m) risMCMC$Theta_m[[m]][t,1:risMCMC$K_Gm[t,m]])
      
      pred_t <- jafar_pred_y_t(Xpred,nPred,M,p_m,risMCMC$K[t],risMCMC$K_Gm[t,],
                               La_m_t,Ga_m_t,Theta_t,Theta_m_t,
                               risMCMC$mu_y[t],risMCMC$s2_inv[t],mu_m_t,s2_m_t,
                               rescale_pred=rescale_pred)   
    } else {
      pred_t <- jfr_pred_y_t(Xpred,nPred,M,p_m,risMCMC$K[t],
                               La_m_t,Theta_t,
                               risMCMC$mu_y[t],risMCMC$s2_inv[t],mu_m_t,s2_m_t,
                               rescale_pred=rescale_pred)   
    }
    
    mu_MC[t,]  <- pred_t
    
    if(t %% iter_print == 0){
      print(sprintf(fmt = "%10s%3s%2s", "[",(t%/%iter_print)*10,"%]"))
    }
  }
  
  return(list(mean=mu_MC))
  
}

# views predictions ----

jfr_pred_X_t <- function(Xpred,nPred,M,p_m,K,Lambda_m,mu_m,s2_inv_m,
                           rescale_pred=rescale_pred){
  
  Xm_pred <- res_m <- s2_m <- La_m <- s2_La_m <- list()
  
  for(m in 1:M){
    
    if(!is.null(Xpred[[m]])){
      
      mar_std_m = rep(1,p_m[m])
      if(rescale_pred){
        mar_std_m = sqrt(1/s2_inv_m[[m]] + rowSums(as.matrix(Lambda_m[[m]]^2)) )
      }
      s2_m[[m]] = s2_inv_m[[m]]*(mar_std_m^2)
      La_m[[m]] = Lambda_m[[m]]/mar_std_m
      
      s2_La_m[[m]] <- s2_m[[m]]*La_m[[m]]
      
      res_m[[m]] <- Xpred[[m]]-rep(1,nPred)%x%t(mu_m[[m]])
    }
  }
  
  for(m in 1:M){
    Q_eta <- diag(1,K,K)
    r_eta <- rep(0,K)
    
    for(mm in c(1:M)[-m]){
      Q_eta <- Q_eta + crossprod(La_m[[mm]],s2_La_m[[mm]]) 
      r_eta <- r_eta + t( res_m[[mm]]%*%s2_La_m[[mm]] )
    }
    
    L_eta    <- t(chol(Q_eta))
    Lr_eta   <- forwardsolve(L_eta, r_eta)
    
    mean_eta <- backsolve(t(L_eta), Lr_eta)
    std_eta  <- backsolve(t(L_eta), matrix(rnorm(K*nPred),K,nPred))
    
    eta_mc   <- t(mean_eta + std_eta)
    
    Xm_pred[[m]] <- rep(1,nPred)%x%t(mu_m[[m]]) +
      tcrossprod(eta_mc,La_m[[m]]) +
      t(sqrt(1/s2_m[[m]])*matrix(rnorm(p_m[m]*nPred),nrow=p_m[m],ncol=nPred))
  }
  
  return(Xm_pred)
}

jafar_pred_X_t <- function(Xpred,nPred,M,p_m,K,K_m,Lambda_m,Gamma_m,mu_m,s2_inv_m,
                           rescale_pred=rescale_pred){
  
  Xm_pred <- res_m <- s2_m <- La_m <- Ga_m <- s2_La_m <- s2_Ga_m <- GaT_s2_La_m <- D_GaT_s2_La_m <- list()
  
  for(m in 1:M){
    
    if(!is.null(Xpred[[m]])){
      
      mar_std_m = rep(1,p_m[m])
      if(rescale_pred){
        mar_std_m = sqrt(1/s2_inv_m[[m]] + rowSums(as.matrix(Lambda_m[[m]]^2)) + rowSums(as.matrix(Gamma_m[[m]]^2)))
      }
      Ga_m[[m]] = Gamma_m[[m]]/mar_std_m
      s2_m[[m]] = s2_inv_m[[m]]*(mar_std_m^2)
      La_m[[m]] = Lambda_m[[m]]/mar_std_m
      
      s2_La_m[[m]] <- s2_m[[m]]*La_m[[m]]
      s2_Ga_m[[m]] <- s2_m[[m]]*Ga_m[[m]]
      GaT_s2_La_m[[m]] = crossprod(s2_Ga_m[[m]],La_m[[m]])
      
      D_m_inv = chol(diag(1.,K_m[m],K_m[m])+crossprod(Ga_m[[m]],s2_Ga_m[[m]]))
      D_GaT_s2_La_m[[m]] = backsolve(D_m_inv,forwardsolve(t(D_m_inv),GaT_s2_La_m[[m]]))
      
      res_m[[m]] <- Xpred[[m]]-rep(1,nPred)%x%t(mu_m[[m]])
    }
  }
  
  for(m in 1:M){
    Q_eta <- diag(1,K,K)
    r_eta <- rep(0,K)
    
    for(mm in c(1:M)[-m]){
      Q_eta <- Q_eta + crossprod(La_m[[mm]],s2_La_m[[mm]]) - crossprod(GaT_s2_La_m[[mm]],D_GaT_s2_La_m[[mm]])
      r_eta <- r_eta + t( res_m[[mm]]%*%(s2_La_m[[mm]] - s2_Ga_m[[mm]]%*%D_GaT_s2_La_m[[mm]]) )
    }
    
    L_eta    <- t(chol(Q_eta))
    Lr_eta   <- forwardsolve(L_eta, r_eta)
    
    mean_eta <- backsolve(t(L_eta), Lr_eta)
    std_eta  <- backsolve(t(L_eta), matrix(rnorm(K*nPred),K,nPred))
    
    eta_mc   <- t(mean_eta + std_eta)
    phi_mc   <- matrix(rnorm(nPred*K_m[m]),nPred,K_m[m])
    
    Xm_pred[[m]] <- rep(1,nPred)%x%t(mu_m[[m]]) +
      tcrossprod(eta_mc,La_m[[m]]) +
      tcrossprod(phi_mc,Ga_m[[m]]) +
      t(sqrt(1/s2_m[[m]])*matrix(rnorm(p_m[m]*nPred),nrow=p_m[m],ncol=nPred))
  }
  
  return(Xm_pred)
}

predict_X <- function(Xpred,risMCMC,rescale_pred=FALSE){
  
  M = length(risMCMC$mu_m)
  p_m = sapply(risMCMC$mu_m,ncol)
  
  nPred = unlist(sapply(Xpred,nrow))[1]
  
  tMCMC = nrow(risMCMC$K_Lm_eff)
  iter_print <- max(1,tMCMC %/% 10)  
  
  Xm_mcmc <- lapply(1:M, function(m) array(NA,c(nPred,p_m[m],tMCMC)))
  
  t_vals = 1:(risMCMC$hyper_param$tBurnIn%/%risMCMC$hyper_param$tThin)
  
  risMCMC$K = risMCMC$K[-t_vals]
  risMCMC$K_Gm = risMCMC$K_Gm[-t_vals,]
  
  print(" - Computing Omics Predictions - ")
  print(sprintf(fmt = "%10s%3s%2s", "[",0,"%]"))
  
  for(t in 1:tMCMC){
    
    La_m_t <- lapply(risMCMC$Lambda_m, function(df) df[t,,1:risMCMC$K[t]])
    mu_m_t <- lapply(risMCMC$mu_m, function(df) df[t,])
    s2_m_t <- lapply(risMCMC$s2_inv_m, function(df) df[t,])
    
    if(grepl('jafar',risMCMC$hyper_param$model)){
      Ga_m_t <- lapply(1:M, function(m) risMCMC$Gamma_m[[m]][t,,1:risMCMC$K_Gm[t,m]])
      
      pred_t <- jafar_pred_X_t(Xpred,nPred,M,p_m,risMCMC$K[t],risMCMC$K_Gm[t,],
                                La_m_t,Ga_m_t,mu_m_t,s2_m_t,rescale_pred=rescale_pred)
    } else {
      pred_t <- jfr_pred_X_t(Xpred,nPred,M,p_m,risMCMC$K[t],
                               La_m_t,mu_m_t,s2_m_t,rescale_pred=rescale_pred) 
    }
    
    for(m in 1:M){Xm_mcmc[[m]][,,t] <- pred_t[[m]]}
    
    if(t %% iter_print == 0){
      print(sprintf(fmt = "%10s%3s%2s", "[",(t%/%iter_print)*10,"%]"))
    }
  }
  
  return(Xm_mcmc)
  
}


