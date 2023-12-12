
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

# Induced Regression Coefficients ----

y_pred_coeff_JAFAR <- function(M,K,K_m,p_m,Theta,s2_inv_y,Lambda_m,Gamma_m,s2_inv_m,rescale_pred=FALSE){
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

# Response Predictions ----

y_pred_JAFAR_from_coeff <- function(Xpred,ris_MCMC){
  
  M <- length(Xpred)
  nPred <- nrow(Xpred[[1]])
  
  tMCMC <- dim(ris_MCMC$K_Gm)[1]
  
  mu_MC  <- matrix(NA,tMCMC,nPred)
  var_MC <- matrix(NA,tMCMC,nPred)
  
  for(t in 1:tMCMC){
    mu_MC[t,] <- rep(ris_MCMC$mu_y[t],nPred)
    for(m in 1:M){
      X_res <- Xpred[[m]]-tcrossprod(rep(1,nPred),ris_MCMC$mu_m[[m]][t,])
      mu_MC[t,] <- mu_MC[t,] + tcrossprod(ris_MCMC$pred_coeff_m[[m]][t,],X_res)
    }
  }
  
  var_MC <- rep(mean(ris_MCMC$pred_var),nPred) + apply(mu_MC,2,var)
  mu_MC <- colMeans(mu_MC)
  
  return(list(mean=mu_MC,var=var_MC))
}

y_cond_pred_JAFAR_NA <- function(Xpred,nPred,M,p_m,K,K_m,
                                   Theta,Lambda_m,Gamma_m,
                                   mu_y,s2_inv_y,mu_m,s2_inv_m,
                                   na_row_idx,na_idx,
                                   rescale_pred=rescale_pred){
  
  # browser()
  # todo <- 'check'
  # browser()
  
  mean_y <- rep(mu_y,nPred)
  var_y  <- rep(1/s2_inv_y,nPred)
  
  Ga_m <- La_m <- s2_Ga_m <- s2_La_m <- C_inv_m <- list()
  GaT_s2_Ga_m <- GaT_s2_La_m <- LaT_s2_La_m <- D_GaT_s2_La_m <- list()

  for(m in 1:M){
    
    if(!is.null(Xpred[[m]])){
    
      mar_std_m = rep(1,p_m[m])
      if(rescale_pred){
        mar_std_m = sqrt(1/s2_inv_m[[m]] + rowSums(Lambda_m[[m]]^2) + rowSums(Gamma_m[[m]]^2))
      }
      s2_m      <- s2_inv_m[[m]]*(mar_std_m^2)
      Ga_m[[m]] <- Gamma_m[[m]][,1:K_m[m],drop=F]/mar_std_m
      La_m[[m]] <- Lambda_m[[m]][,1:K,drop=F]/mar_std_m
      
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
          
          # browser()
          # todo <- 'check'
          # browser()
          
          # La_D_Xi <- La_D_Xi + colSums(s2_La_m[[m]]*(Xpred[[m]][i,]-mu_m[[m]])) -
          #   colSums(D_GaT_s2_La_m[[m]]*colSums(s2_Ga_m[[m]]*(Xpred[[m]][i,]-mu_m[[m]])))

          La_D_Xi <- La_D_Xi + crossprod(s2_La_m[[m]],Xpred[[m]][i,]-mu_m[[m]]) -
            crossprod(D_GaT_s2_La_m[[m]],crossprod(s2_Ga_m[[m]],Xpred[[m]][i,]-mu_m[[m]]))
          
          # browser()
          # todo <- 'check'
          # browser()
          
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
          
          # browser()
          # todo <- 'check'
          # browser()
          
          La_D_Xi <- La_D_Xi + crossprod(s2_La_m[[m]][-idx_m_i,,drop=F],
              Xpred[[m]][i,-idx_m_i]-mu_m[[m]][-idx_m_i]) -
            crossprod(D_GaT_s2_La_m_NA,crossprod(s2_Ga_m[[m]][-idx_m_i,,drop=F],
              Xpred[[m]][i,-idx_m_i]-mu_m[[m]][-idx_m_i]))
          
          # browser()
          # todo <- 'check'
          # browser()
      
        }
      }
    }
    
    # browser()
    # todo <- 'check'
    # browser()
    
    C_chol = chol(Ci_inv)
    Theta_C = backsolve(C_chol,forwardsolve(t(C_chol),Theta[1:K,drop=F]))
    
    var_y[i]  <- var_y[i] + sum(Theta_C*Theta[1:K,drop=F])
    mean_y[i] <- mean_y[i] + sum(Theta_C*La_D_Xi)
    
    # browser()
    # todo <- 'check'
    # browser()
  }
  
  return(list(pred_mean=mean_y, pred_var=var_y))
}

y_pred_JAFAR_NA <- function(Xpred,risMCMC,rescale_pred=FALSE){

  M = length(Xpred)
  nPred = unlist(sapply(Xpred,nrow))[1]
  p_m = sapply(Xpred,ncol)
  
  get_NA_pred <- get_NA_X(Xpred)
  na_row_idx  <- get_NA_pred$na_row_idx
  na_idx      <- get_NA_pred$na_idx
  
  tMCMC = dim(risMCMC$K_Gm)[1]
  iter_print <- tMCMC %/% 10 
  
  var_MC <- mu_MC  <- matrix(NA,tMCMC,nPred)

  print(" - Computing Response Predictions - ")
  print(sprintf(fmt = "%10s%3s%2s", "[",0,"%]"))
  
  for(t in 1:tMCMC){
    
    La_m_t <- lapply(risMCMC$Lambda_m, function(df) df[t,,])
    Ga_m_t <- lapply(risMCMC$Gamma_m, function(df) df[t,,])
    mu_m_t <- lapply(risMCMC$mu_m, function(df) df[t,])
    s2_m_t <- lapply(risMCMC$s2_inv_m, function(df) df[t,])
    
    pred_t <- y_cond_pred_JAFAR_NA(Xpred,nPred,M,p_m,
      risMCMC$K[t],risMCMC$K_Gm[t,],
      risMCMC$Theta[t,],La_m_t,Ga_m_t,
      risMCMC$mu_y[t],risMCMC$s2_inv[t],mu_m_t,s2_m_t,
      na_row_idx,na_idx,rescale_pred=rescale_pred)
    
    var_MC[t,] <- pred_t$pred_var
    mu_MC[t,]  <- pred_t$pred_mean
    
    if(t %% iter_print == 0){
      print(sprintf(fmt = "%10s%3s%2s", "[",(t%/%iter_print)*10,"%]"))
    }
  }
  
  var_MC <- colMeans(var_MC) + apply(mu_MC,2,var)
  mu_MC  <- colMeans(mu_MC)
  
  return(list(mean=mu_MC,var=var_MC))
  
}

# Omics Predictions ----

X_cond_samples_JAFAR_NA <- function(Xpred,nPred,M,p_m,K,K_m,Lambda_m,Gamma_m,mu_m,s2_inv_m,
                                    na_row_idx,na_idx,rescale_pred=rescale_pred,nSample=10){
  
  X_sample <- lapply(1:M, function(m) array(NA,c(nSample,nPred,p_m[m])))
  
  Ga_m <- La_m <- s2_m <- s2_Ga_m <- s2_La_m <- C_inv_m <- list()
  GaT_s2_Ga_m <- GaT_s2_La_m <- LaT_s2_La_m <- D_GaT_s2_La_m <- list()
  
  for(m in 1:M){
    
    mar_std_m = rep(1,p_m[m])
    if(rescale_pred){
      mar_std_m = sqrt(1/s2_inv_m[[m]] + rowSums(Lambda_m[[m]]^2) + rowSums(Gamma_m[[m]]^2))
    }
    s2_m[[m]] <- s2_inv_m[[m]]*(mar_std_m^2)
    Ga_m[[m]] <- Gamma_m[[m]][,1:K_m[m],drop=F]/mar_std_m
    La_m[[m]] <- Lambda_m[[m]][,1:K,drop=F]/mar_std_m
    
    s2_Ga_m[[m]] <- s2_m[[m]]*Ga_m[[m]]
    s2_La_m[[m]] <- s2_m[[m]]*La_m[[m]]
    
    GaT_s2_Ga_m[[m]] <- crossprod(s2_Ga_m[[m]],Ga_m[[m]])
    GaT_s2_La_m[[m]] <- crossprod(s2_Ga_m[[m]],La_m[[m]])
    LaT_s2_La_m[[m]] <- crossprod(s2_La_m[[m]],La_m[[m]])
    
    D_m_chol0 <- chol(diag(1.,K_m[m],K_m[m]) + GaT_s2_Ga_m[[m]])
    
    D_GaT_s2_La_m[[m]] <- backsolve(D_m_chol0,forwardsolve(t(D_m_chol0),GaT_s2_La_m[[m]]))
    
    C_inv_m[[m]] <- LaT_s2_La_m[[m]] - crossprod(GaT_s2_La_m[[m]],D_GaT_s2_La_m[[m]])
  }
  
  for(i in 1:nPred){
    
    for(mPred in 1:M){
      
      Ci_inv  <- diag(1,K,K)
      La_D_Xi <- rep(0.,K)
      
      for(m in c(1:M)[-mPred]){
        
        idx_i <- match(i,na_row_idx[[m]])
        
        if(is.na(idx_i)){
          Ci_inv  <- Ci_inv + C_inv_m[[m]]
          
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
      
      C_chol = chol(Ci_inv)
      C_La_D_Xi = backsolve(C_chol,forwardsolve(t(C_chol),La_D_Xi))
      
      X_sample[[mPred]][,i,] = matrix(mu_m[[mPred]] + as.vector(La_m[[mPred]]%*%C_La_D_Xi),nSample,p_m[mPred]) +
        t(La_m[[mPred]]%*%forwardsolve(t(C_chol),matrix(rnorm(K*nSample),K,nSample)) +
          Ga_m[[mPred]]%*%matrix(rnorm(K_m[mPred]*nSample),K_m[mPred],nSample) +
          s2_m[[mPred]]*matrix(rnorm(p_m[mPred]*nSample),p_m[mPred],nSample) )
      
    }
  }
  
  return(X_sample)
}

X_sample_JAFAR_NA <- function(Xpred,risMCMC,Xtrain=NULL,rescale_pred=FALSE,
                              copula=F,smoothed=F,nSample=10){
  
  nPred = unlist(sapply(Xpred,nrow))[1]
  p_m = sapply(risMCMC$s2_inv_m,ncol)
  M = ncol(risMCMC$K_Gm)
  
  Z_m <- Q_smooth <- list()
  if(copula){
    for(m in 1:M){
      if(smoothed){
        Q_smooth[[m]] <- lapply(1:p_m[m], function(j) get_Q_smooth(Xtrain[[m]][,j]))
        F_smooth      <- lapply(1:p_m[m], function(j) get_F_smooth(Xtrain[[m]][,j]))
        Z_m[[m]] <- qnorm(sapply(1:p_m[m], function(j) F_smooth[[j]](Xpred[[m]][,j]) ))
      } else {
        Z_m[[m]] <- qnorm(sapply(1:p_m[m], function(j) F_hat(Xpred[[m]][,j], Xtrain[[m]][,j])))
      }
    }
  } else {
    Z_m <- Xpred
  }
  
  get_NA_pred <- get_NA_X(Xpred)
  na_row_idx  <- get_NA_pred$na_row_idx
  na_idx      <- get_NA_pred$na_idx
  
  tMCMC = nrow(risMCMC$K_Gm)
  iter_print <- tMCMC %/% 10 
  
  Z_MC <- lapply(1:M, function(m) array(NA,c(tMCMC*nSample,nPred,p_m[m])))
  
  print(" - Computing Omics Imputations - ")
  print(sprintf(fmt = "%10s%3s%2s", "[",0,"%]"))
  
  for(t in 1:tMCMC){
    
    La_m_t <- lapply(risMCMC$Lambda_m, function(df) df[t,,])
    Ga_m_t <- lapply(risMCMC$Gamma_m, function(df) df[t,,])
    mu_m_t <- lapply(risMCMC$mu_m, function(df) df[t,])
    s2_m_t <- lapply(risMCMC$s2_inv_m, function(df) df[t,])
    
    Z_sample <- X_cond_samples_JAFAR_NA(Z_m,nPred,M,p_m,
                                        risMCMC$K[t],risMCMC$K_Gm[t,],
                                        La_m_t,Ga_m_t,mu_m_t,s2_m_t,
                                        na_row_idx,na_idx,nSample=nSample,
                                        rescale_pred=rescale_pred)
    
    for(m in 1:M){ Z_MC[[m]][(1+(t-1)*nSample):(t*nSample),,] <- Z_sample[[m]] }
    
    if(t %% iter_print == 0){
      print(sprintf(fmt = "%10s%3s%2s", "[",(t%/%iter_print)*10,"%]"))
    }
  }
  
  if(copula){
    for(m in 1:M){
      if(smoothed){
        for(j in 1:p_m[m]){ Z_MC[[m]][,,j] <- matrix(Q_smooth[[m]][[j]](Z_MC[[m]][,,j]),ncol=nPred) }
      } else {
        for(j in 1:p_m[m]){ Z_MC[[m]][,,j] <- Q_hat(pnorm(Z_MC[[m]][,,j]),Xtrain[[m]][,j]) }
      }
    }
  }
  
  return(Z_MC)
}

# Burn-In Removal Ex-Post ----

remove_BurnIn <- function(risMCMC,tBurn){
  
  t_idx = c(1:tBurn)

  risMCMC$K <- risMCMC$K[-t_idx]
  risMCMC$K_T_eff <- risMCMC$K_T_eff[-t_idx]
  risMCMC$active_T <- risMCMC$active_T[-t_idx,]
  risMCMC$K_Gm <- risMCMC$K_Gm[-t_idx,]
  risMCMC$K_Lm_eff <- risMCMC$K_Lm_eff[-t_idx,]
  risMCMC$K_Gm_eff <- risMCMC$K_Gm_eff[-t_idx,]
  risMCMC$active_Lm <- risMCMC$active_Lm[-t_idx,,]
  
  risMCMC$Theta <- risMCMC$Theta[-t_idx,]
  risMCMC$var_y <- risMCMC$var_y[-t_idx]
  risMCMC$s2_inv <- risMCMC$s2_inv[-t_idx]
  risMCMC$mu_y <- risMCMC$mu_y[-t_idx]
  risMCMC$pred_var <- risMCMC$pred_var[-t_idx]
  risMCMC$eta <- risMCMC$eta[-t_idx,,]
  
  risMCMC$phi_m <- lapply(risMCMC$phi_m, function(df) df[-t_idx,,])
  risMCMC$Lambda_m <- lapply(risMCMC$Lambda_m, function(df) df[-t_idx,,])
  risMCMC$Gamma_m <- lapply(risMCMC$Gamma_m, function(df) df[-t_idx,,])
  risMCMC$s2_inv_m <- lapply(risMCMC$s2_inv_m, function(df) df[-t_idx,])
  risMCMC$mu_m <- lapply(risMCMC$mu_m, function(df) df[-t_idx,])
  risMCMC$Marg_Var_m <- lapply(risMCMC$Marg_Var_m, function(df) df[-t_idx,])
  risMCMC$pred_coeff_m <- lapply(risMCMC$pred_coeff_m, function(df) df[-t_idx,])
  
  if('Xm_MC' %in% names(risMCMC)){
    risMCMC$Xm_MC <- lapply(risMCMC$Xm_MC, function(ll) lapply(ll,function(df) df[-t_idx,]))
  }
  
  return(risMCMC)  
}









