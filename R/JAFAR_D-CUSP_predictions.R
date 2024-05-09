
# Missing Data Identification ----

get_NA_X <- function(X_m){
  
  M <- length(X_m)
  n <- unlist(sapply(X_m,nrow))[1]
  
  na_idx <- na_row_idx <- list()
  
  for(m in 1:M){
    
    if(!is.null(X_m[[m]])){
      
      # identifying indexes of NA
      Xm_na <- is.na(unname(as.matrix(X_m[[m]])))
      
      # Containers for output of imputed NA values
      na_idx[[m]]     <- apply(Xm_na,1,which)
      na_row_idx[[m]] <- c(1:n)[sapply(na_idx[[m]],length)>0]
      na_idx[[m]]     <- na_idx[[m]][na_row_idx[[m]]]
      
    }
  }
  
  return(list(na_row_idx=na_row_idx,na_idx=na_idx))  
}

# Copula Functions ----

order_index = function(x){
  ind = sort(x, index.return = T)$ix
  x[ind] = (1:length(x))/(length(x)+1)
  return(x)
}

order_index_na = function(x){
  x_na = is.na(x)
  x[x_na == F] = order_index(x[x_na == F])
  n = length(x) - sum(x_na)
  return(x*n/(n+1))
}

F_hat <- function(xNew,xTrain){
  x_obs <- xTrain[!is.na(xTrain)]
  n = length(x_obs)
  knots = c(-Inf,sort(x_obs))
  values = c(1/(n+1),c(1:n)/(n+1))
  
  sapply(xNew,function(x) max(values[knots<=x]))
}

Q_hat <- function(pNew,xTrain) {
  x_obs <- xTrain[!is.na(xTrain)]
  n <- length(x_obs)
  knots = c(sort(x_obs),max(x_obs))
  values = c(1:n)/(n+1)
  
  knots[sapply(pNew, function(p) sum(p > values) + 1)]
}

get_F_smooth <- function(vec){
  vec <- vec[!is.na(vec)]
  x_range  <- c(-100,sort(vec),100)
  # y_range  <- c(0,sort(order_index_na(vec)),1)
  y_range  <- c(0.001,sort(order_index_na(vec)),0.999)
  F_spline <- stats::splinefun(x=x_range, y=y_range, method = "hyman", ties= list("ordered", mean))
  return(F_spline)
}

get_Q_smooth <- function(vec){
  vec <- vec[!is.na(vec)]
  x_range  <- c(0,sort(order_index(vec)),1)
  y_range  <- c(min(vec),sort(vec),max(vec))
  Q_spline <- stats::splinefun(x=x_range, y=y_range, method = "hyman", ties= list("ordered", mean))
  return(Q_spline)
}

cdf_transform <- function(Z_m,Z_m_test=NULL,smoothed=F){
  preprocess_X_m <- list()
  for(m in 1:Data$M){
    # Train Set
    Z_m[[m]] = qnorm(apply(Z_m[[m]], 2, order_index_na))
    # Center and Scale Predictors
    preprocess_X_m[[m]] = caret::preProcess(Z_m[[m]], method = c("center", "scale"))
    Z_m[[m]] = as.matrix(predict(preprocess_X_m[[m]], Z_m[[m]]))
    # Test Set
    if(!is.null(Z_m_test)){
      for(j in 1:Data$p_m[m]){
        if(smoothed){
          F_smooth <- get_F_smooth(Z_m[[m]][,j])
          Z_m_test[[m]][,j] <- qnorm(F_smooth(Z_m_test[[m]][,j]))
        } else{
          Z_m_test[[m]][,j] <- qnorm(F_hat(Z_m_test[[m]][,j], Z_m[[m]][,j]))
        }
      }
      Z_m_test[[m]] = as.matrix(predict(preprocess_X_m, Z_m_test[[m]]))
    }
  }
  output <- list(Z_m=Z_m,preprocess_X_m=preprocess_X_m)
  if(!is.null(Z_m_test)){output$Z_m_test=Z_m_test}
  return(output)
}


# Induced Regression Coefficients ----

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

# Prediction of Linear Predictor via Sampling of Latent Factors ----


y_cond_pred_JAFAR_sampling <- function(Xpred,nPred,M,p_m,K,K_m,
                                       Theta,Lambda_m,Gamma_m,mu_y,s2_inv_y,mu_m,s2_inv_m,
                                       rescale_pred=rescale_pred){
  
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


# Linear Predictor via Exact Expression with missing data in features -----

y_cond_pred_JAFAR_NA <- function(Xpred,nPred,M,p_m,K,K_m,
                                 Theta,Lambda_m,Gamma_m,
                                 mu_y,s2_inv_y,mu_m,s2_inv_m,
                                 na_row_idx,na_idx,
                                 rescale_pred=rescale_pred){
  
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

# Inferred Linear Predictor for out-of-sample observations -----

y_pred_JAFAR <- function(Xpred,risMCMC,rescale_pred=FALSE){
  
  M = length(Xpred)
  nPred = unlist(sapply(Xpred,nrow))[1]
  p_m = sapply(Xpred,ncol)
  
  NA_in_X <- max(sapply(X_m,function(df) max(is.na(df))))
  
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
    
    Ga_m_t  <- lapply(1:M, function(m) risMCMC$Gamma_m[[m]][t,,1:risMCMC$K_Gm[t,m],drop=F])
    La_m_t  <- lapply(risMCMC$Lambda_m, function(df) df[t,,1:risMCMC$K[t],drop=F])
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





















