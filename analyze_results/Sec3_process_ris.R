
if(FALSE){
  rm(list = ls())
}  

if(FALSE){
  library(fields)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(multiview)
  library(aperm)
}
 
if(T){ 
  credible_intervals_mcmc <- function(m_samples, s_samples, Smc=100){
    T_val <- nrow(m_samples)
    N_val <- ncol(m_samples)
    if(length(s_samples)==1){s_samples=rep(s_samples,T_val)}
    
    means_matrix_long <- m_samples[, rep(1:N_val, each = Smc)]
    sds_matrix_long <- matrix(rep(s_samples, each = N_val * Smc), nrow = T_val)
    
    noise_matrix <- matrix(rnorm(T_val * N_val * Smc), nrow = T_val)
    
    all_samples <- means_matrix_long + sds_matrix_long * noise_matrix
    
    all_samples_reshaped <- matrix(all_samples, ncol = N_val)
    
    # which_quantiles=c(0.025, 0.975)
    which_quantiles=c(0.05, 0.95)
    
    credible_intervals <- apply(all_samples_reshaped, 2, quantile, probs = which_quantiles)
    
    t(credible_intervals)
  }
}

data_path = '~/Documents/GitHub/jafar/data/Sec3_Simulations/'
ris_pathJAFAR = '~/Downloads/jafar_ris_paper/sim_sec3_J/'
ris_pathJFR = '~/Downloads/jafar_ris_paper/sim_sec3_jfr/'
ris_pathBSFP = '~/Downloads/jafar_ris_paper/sim_sec3_bsfp/'
ris_pathBIP = '~/Downloads/jafar_ris_paper/sim_sec3_bip//'
ris_pathCOOP = '~/Documents/GitHub/jafar/ris/sim_sec3/'

# LOAD RESULTS -----------
if(TRUE){
  
  source('~/Documents/GitHub/jafar/new_source/Sec3_simulations_preproc.R') # data  
  
  # Response predicitons
  mse_train <- mse_test <- matrix(NA,nrow=0,ncol=6)
  colnames(mse_train) <- colnames(mse_test) <- c('n','r','method','mse_y','coverage90s','R_squared')

  # ess and nrun 
  ess_mcmc <- matrix(NA,nrow=0,ncol=8)
  colnames(ess_mcmc) <- c('n','r','method','nTot','nBurn','nThin','ess_s2','ess_mu')
  
  ess_views <- matrix(NA,nrow=0,ncol=6)
  colnames(ess_views) <- c('n','r','method','m','mean_ess','median_ess')
    
  # time 
  time_run <- matrix(NA,nrow=0,ncol=5)
  colnames(time_run) <- c('n','r','method','timeRun','timePred')
  
  # n. of factors
  n_fact <- matrix(NA,nrow=0,ncol=7)
  colnames(n_fact) <- c('n','r','method','K','K1','K2','K3')
  
  # Correlations 
  cor_fn2 <- matrix(NA,nrow=0,ncol=6)
  colnames(cor_fn2) <- c('n','r','method','m','exact','cor_fn2')
  
  n_values = c(50,100,150,200)
  s_values = c(1:100)
  # n_values = c(200)
  # s_values = c(1:3) # c(1:10)*10
  
  n0=200
  
  v_meth = c('JAFAR','JFR','BSFP','BIP','IntegL','CoopL')
  
  for(nn in n_values){
    for(ss in s_values){
      
      print(paste0(' --- Loading Results and Data for n=',nn,' s=',ss, ' --- '))
      
      dataFile <- paste0('Simulated_data_n',n0,'_s',ss,'.rds')
      Data <- readRDS(file.path(data_path, dataFile))
      
      if(TRUE){
        # Data ----
        Data <- data_scale_subset(Data,nn)
      }
      
      cor_true_m <- list()
      for(m in 1:Data$M){  
        cor_true_m[[m]] <- cov2cor(tcrossprod(Data$Lambda_m[[m]])+tcrossprod(Data$Gamma_m[[m]])+diag(Data$s2m[[m]]))
      }
      
      for(meth in 1:length(v_meth)){
        
        y_pred_train <- y_pred_test <- sig_y_mcmc <- NULL
        
        ## BIP ------
        if(v_meth[meth]=='BIP'){
          
          which_rep=0
          risFile <- paste0('Sec3_Simulated_data_n',nn,'_s',ss,'_y_bip_nMC6000_nBurn3000_rep',which_rep,'.rds')
          risBIP <- tryCatch(readRDS(file.path(ris_pathBIP, risFile)), 
                             error = function(ee) NULL, warning= function(ww) NULL)
          if (!is.null(risBIP)) {
            
            # Response predictions
            y_pred_train = matrix(risBIP$train_pred$ypredict,nrow=1)
            y_pred_test = matrix(risBIP$test_pred$ypredict,nrow=1)
            # sig_y_mcmc <- risBIP$runBIP$EstSig2[[Data$M+1]]
            
            # time 
            time_run <- rbind(time_run, c(nn,ss,meth,risBIP$time_run,risBIP$time_pred))
            
            # n. of factors
            n_fact <- rbind(n_fact, c(nn,ss,meth,ncol(risBIP$runBIP$EstU),rep(0,Data$M)))
            
            # Correlations 
            for(m in 1:Data$M){
              cor_fn2 <- rbind(cor_fn2, c(nn,ss,meth,m,0,sum((cor_true_m[[m]]-cov2cor(risBIP$cov_mean_m[[m]]))^2) / (Data$p_m[m])^2))
            }
            
            rm(risBIP)
          } else {
            print(paste0('BIP ris NOT available for n=',nn,' s=',ss))
          }
        }
        
        # BSFP ------
        if(v_meth[meth]=='BSFP'){
          
          which_rep=0
          risFile <- paste0('Sec3_Simulated_data_n',nn,'_s',ss,'_y_bsfp_nMC6000_nBurn3000_nThin10_rep',which_rep,'.rds')
          risBSFP <- tryCatch(readRDS(file.path(ris_pathBSFP, risFile)), 
                             error = function(ee) NULL, warning= function(ww) NULL)
          if (!is.null(risBSFP)) {
          
            # Response predictions
            y_pred_train = t(do.call(cbind,do.call(c,risBSFP$train_pred$EY.draw)))
            y_pred_train = y_pred_train[-c(1:floor(nrow(y_pred_train)/2)),]
            
            y_pred_test = t(do.call(cbind,do.call(c,risBSFP$test_pred$EY.draw)))
            y_pred_test = y_pred_test[-c(1:floor(nrow(y_pred_test)/2)),]
            
            sig_y_mcmc <- unlist(risBSFP$bsfp_mcmc$tau2.draw)
            sig_y_mcmc <- sqrt(tail(sig_y_mcmc,floor(length(sig_y_mcmc)/2)))
            
            # time 
            time_run <- rbind(time_run, c(nn,ss,meth,risBSFP$time_run,risBSFP$time_pred))
            
            # n. of factors
            n_fact <- rbind(n_fact, c(nn,ss,meth,risBSFP$bsfp_mcmc$ranks))
            
            # ess and nrun 
            s2_y = unlist(risBSFP$bsfp_mcmc$tau2.draw)[-c(1:300)]
            mu_y = matrix(unlist(risBSFP$bsfp_mcmc$beta.draw),ncol=600)[1,-c(1:300)]
            ess_mcmc <- rbind(ess_mcmc, c(nn,ss,meth,6000,3000,10,sns::ess(s2_y),sns::ess(mu_y)))
            
            ess_m_s2 = lapply(risBSFP$Marg_Var_m, function(mat) apply(mat,2,sns::ess))
            for(m in 1:Data$M){
              ess_views <- rbind(ess_views, c(nn,ss,meth,m,mean(ess_m_s2[[m]]),median(ess_m_s2[[m]])))
            }
            
            # Correlations 
            for(m in 1:Data$M){
              cor_fn2 <- rbind(cor_fn2, c(nn,ss,meth,m,0,sum((cor_true_m[[m]]-cov2cor(risBSFP$Cov_m_mean[[m]]))^2) / (Data$p_m[m])^2))
            }
            
            rm(risBSFP)
          } else {
            print(paste0('BSFP ris NOT available for n=',nn,' s=',ss))
          }
        }
        
        # JAFAR ------
        if(v_meth[meth]=='JAFAR'){
          
          which_rep=0 # 
          risFile <- paste0('Sec3_Simulated_data_n',nn,'_s',ss,'_y_i-cusp_nMC20000_nBurn15000_nThin10_rep',which_rep,'.rds')
          risJAFAR <- tryCatch(readRDS(file.path(ris_pathJAFAR, risFile)), 
                              error = function(ee) NULL, warning= function(ww) NULL)
          if (!is.null(risJAFAR)) {
            
            # Response predictions
            y_pred_train <- risJAFAR$y_JAFAR_train$mean
            y_pred_test  <- risJAFAR$y_JAFAR_test$mean
            
            sig_y_mcmc <- 1/sqrt(risJAFAR$ris_MCMC$s2_inv)
            
            # time 
            time_run <- rbind(time_run, c(nn,ss,meth,risJAFAR$time_run,risJAFAR$time_pred))
            
            # n. of factors
            n_fact <- rbind(n_fact, c(nn,ss,meth,ncol(risJAFAR$ris_MCMC$Theta),sapply(risJAFAR$ris_MCMC$Theta_m,ncol)))
            
            # ess and nrun
            ess_mcmc <- rbind(ess_mcmc, c(nn,ss,meth,20000,15000,10,sns::ess(1/risJAFAR$ris_MCMC$s2_inv),sns::ess(risJAFAR$ris_MCMC$mu_y)))
            
            ess_m_s2 = lapply(risJAFAR$Marg_Var_m, function(mat) apply(mat,2,sns::ess))
            for(m in 1:Data$M){
              ess_views <- rbind(ess_views, c(nn,ss,meth,m,mean(ess_m_s2[[m]]),median(ess_m_s2[[m]])))
            }
            
            # Correlations 
            for(m in 1:Data$M){
              cor_fn2 <- rbind(cor_fn2, c(nn,ss,meth,m,1,sum((cor_true_m[[m]]-cov2cor(risJAFAR$ris_MCMC$Cov_m_mean[[m]]))^2) / (Data$p_m[m])^2))
              cor_fn2 <- rbind(cor_fn2, c(nn,ss,meth,m,0,sum((cor_true_m[[m]]-cov2cor(risJAFAR$Cov_m_mean[[m]]))^2) / (Data$p_m[m])^2))
            }
            
            rm(risJAFAR)
          } else {
            print(paste0('JAFAR ris NOT available for n=',nn,' s=',ss))
          }
        }
        
        # JFR ------
        if(v_meth[meth]=='JFR'){
          
          which_rep=0 # 
          risFile <- paste0('Sec3_Simulated_data_n',nn,'_s',ss,'_y_jfr_nMC20000_nBurn15000_nThin10_rep',which_rep,'.rds')
          risJFR <- tryCatch(readRDS(file.path(ris_pathJFR, risFile)), 
                             error = function(ee) NULL, warning= function(ww) NULL)
          if (!is.null(risJFR)) {
            
            # Response predictions
            y_pred_train <- risJFR$y_JAFAR_train$mean
            y_pred_test  <- risJFR$y_JAFAR_test$mean
            
            sig_y_mcmc <- 1/sqrt(risJFR$ris_MCMC$s2_inv)
            
            # time 
            time_run <- rbind(time_run, c(nn,ss,meth,risJFR$time_run,risJFR$time_pred))
            
            # n. of factors
            n_fact <- rbind(n_fact, c(nn,ss,meth,ncol(risJFR$ris_MCMC$Theta),rep(0,Data$M)))
            
            # ess and nrun
            ess_mcmc <- rbind(ess_mcmc, c(nn,ss,meth,20000,15000,10,sns::ess(1/risJFR$ris_MCMC$s2_inv),sns::ess(risJFR$ris_MCMC$mu_y)))
            
            ess_m_s2 = lapply(risJFR$Marg_Var_m, function(mat) apply(mat,2,sns::ess))
            for(m in 1:Data$M){
              ess_views <- rbind(ess_views, c(nn,ss,meth,m,mean(ess_m_s2[[m]]),median(ess_m_s2[[m]])))
            }
            
            # Correlations 
            for(m in 1:Data$M){
              cor_fn2 <- rbind(cor_fn2, c(nn,ss,meth,m,1,sum((cor_true_m[[m]]-cov2cor(risJFR$ris_MCMC$Cov_m_mean[[m]]))^2) / (Data$p_m[m])^2))
              cor_fn2 <- rbind(cor_fn2, c(nn,ss,meth,m,0,sum((cor_true_m[[m]]-cov2cor(risJFR$Cov_m_mean[[m]]))^2) / (Data$p_m[m])^2))
            }
            
            rm(risJFR)
          } else {
            print(paste0('JFR ris NOT available for n=',nn,' s=',ss))
          }
        }
        
        # IntegL -----
        if(v_meth[meth]=='IntegL'){
          
          which_rep=0 # 
          risFile <- paste0('Sec3_Simulated_data_n',nn,'_s',ss,'_y_IntegL_rep',which_rep,'.rds')
          risIntegL <- tryCatch(readRDS(file.path(ris_pathCOOP, risFile)), 
                             error = function(ee) NULL, warning= function(ww) NULL)
          if (!is.null(risIntegL)) {
            
            # Response predictions
            y_pred_train <- t(risIntegL$train_pred)
            y_pred_test <- t(risIntegL$test_pred)
            sig_y_mcmc = 1
            
            # time 
            time_run <- rbind(time_run, c(nn,ss,meth,risIntegL$time_run,risIntegL$time_pred))
            
            rm(risIntegL)
          } else {
            print(paste0('Integ_Learn ris NOT available for n=',nn,' s=',ss))
          }
        }
        
        # CoopL -----
        if(v_meth[meth]=='CoopL'){
          
          which_rep=0
          risFile <- paste0('Sec3_Simulated_data_n',nn,'_s',ss,'_y_CoopL_rep',which_rep,'.rds')
          risCoopL <- tryCatch(readRDS(file.path(ris_pathCOOP, risFile)), 
                                error = function(ee) NULL, warning= function(ww) NULL)
          if (!is.null(risCoopL)) {
            
            # Response predictions
            y_pred_train <- matrix(risCoopL$y_pred_train,nrow=1)
            y_pred_test <- matrix(risCoopL$y_pred_test,nrow=1)
            
            # time 
            time_run <- rbind(time_run, c(nn,ss,meth,risCoopL$time_run,risCoopL$time_pred))
            
            rm(risCoopL)
          } else {
            print(paste0('CoopL ris NOT available for n=',nn,' s=',ss))
          }
        }
        
        # response train set ----
        if(!is.null(y_pred_train)){
          y_mean = colMeans(y_pred_train)
          if(v_meth[meth] %in% c('JAFAR','JFR','BSFP','IntegL')){
            train_uq <- credible_intervals_mcmc(y_pred_train,sig_y_mcmc, Smc=1)
            y_low <- train_uq[,1]
            y_up  <- train_uq[,2]
            
            mse_train <- rbind(mse_train, c(nn,ss,meth,mean((Data$yTrain-y_mean)^2),
                                            mean((Data$yTrain<=y_up)*(Data$yTrain>=y_low)),
                                            1-sum((Data$yTrain-y_mean)^2)/sum((Data$yTrain-mean(Data$yTrain))^2) ))
          } else if(v_meth[meth] %in% c('BIP','CoopL')){
            mse_train <- rbind(mse_train, c(nn,ss,meth,mean((Data$yTrain-y_mean)^2),
                                          NA, 1-sum((Data$yTrain-y_mean)^2)/sum((Data$yTrain-mean(Data$yTrain))^2) ))
          }
        }
        
        # response test set ----
        if(!is.null(y_pred_test)){
          y_mean = colMeans(y_pred_test)
          if(v_meth[meth] %in% c('JAFAR','JFR','BSFP','IntegL')){
            test_uq <- credible_intervals_mcmc(y_pred_test,sig_y_mcmc, Smc=1)
            y_low <- test_uq[,1]
            y_up  <- test_uq[,2]
            
            mse_test <- rbind(mse_test, c(nn,ss,meth,mean((Data$yTest-y_mean)^2),
                                            mean((Data$yTest<=y_up)*(Data$yTest>=y_low)),
                                            1-sum((Data$yTest-y_mean)^2)/sum((Data$yTest-mean(Data$yTest))^2) ))
          } else if(v_meth[meth] %in% c('BIP','CoopL')){
            mse_test <- rbind(mse_test, c(nn,ss,meth,mean((Data$yTest-y_mean)^2),
                                            NA, 1-sum((Data$yTest-y_mean)^2)/sum((Data$yTest-mean(Data$yTest))^2) ))
          }
        }
        
      }
      
      
    }
  }
  
  # save results ----
  ris_list = list(mse_train=mse_train, mse_test=mse_test, ess_mcmc=ess_mcmc,
                  ess_views=ess_views, time_run=time_run, n_fact=n_fact, cor_fn2=cor_fn2)
                  
  saveRDS(ris_list, file.path('~/Desktop/', 'Sec3_Simulations_results.rds'))
}



