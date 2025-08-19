
if(F){
  rm(list = ls())

  library(dplyr)
  library(ggplot2)
  library(scales)
  library(dplyr)
  library(tidyr)
}

credible_intervals_mcmc <- function(m_samples, s_samples, Smc=100) {
  T_val <- nrow(m_samples)
  N_val <- ncol(m_samples)
  if(length(s_samples)==1){s_samples=rep(s_samples,T_val)}
  
  means_matrix_long <- m_samples[, rep(1:N_val, each = Smc)]
  sds_matrix_long <- matrix(rep(s_samples, each = N_val * Smc), nrow = T_val)
  
  noise_matrix <- matrix(rnorm(T_val * N_val * Smc), nrow = T_val)
  
  all_samples <- means_matrix_long + sds_matrix_long * noise_matrix
  
  all_samples_reshaped <- matrix(all_samples, ncol = N_val)
  
  # credible_intervals <- apply(all_samples_reshaped, 2, quantile, probs = c(0.025, 0.975))
  credible_intervals <- apply(all_samples_reshaped, 2, quantile, probs = c(0.05, 0.95))
  
  t(credible_intervals)
}
  
ggplot_predictions_NEW <- function(yTrue, pred_means, pred_low, pred_up, pred_names, file_name, our_dir, plt_w=10,
                                  pt_sz = 0.7, rf_sz = 0.3, br_lw = 0.5, br_wd = 0.1, br_al = 0.5, tc_al = 0.5, sh_fct = 1) {
  
  # Extended color and shape palettes for up to 8 methods
  v_colors <- c('#FFB000',"#DC267F","#FE6100","#648FFF","#2DC4DC",'#2B5B43','#57B38C','#8D2FE5')
  v_pch <- c(21, 22, 23, 24, 25, 15, 17, 18)
  
  idx <- sort(yTrue, index.return = T)$ix
  
  n_meth <- length(pred_names)
  n_pred <- length(yTrue)
  
  Lup <- length(pred_up)
  Llow <- length(pred_low)
  
  # Updated shift logic for up to 8 methods
  if (n_meth > 8) {
    stop("Function developed for up to 8 methods")
  }
  
  # Calculate shifts to avoid hardcoding for each number of methods
  half_width <- 0.25 * (n_meth - 1)
  shift <- seq(-half_width, half_width, length.out = n_meth) * 10 / n_pred * sh_fct
  
  # Create an empty data frame to store all method predictions
  df_plot <- data.frame(
    x_column = numeric(),
    y_column = numeric(),
    meth = character(),
    y_low = numeric(),
    y_up = numeric()
  )
  
  # Loop through all methods to build the data frame
  for (l in 1:n_meth) {
    df_tmp <- data.frame(x_column = c(1:n_pred) + shift[l], y_column = pred_means[[l]][idx])
    df_tmp$meth <- pred_names[l]
    df_tmp$y_low <- NA
    df_tmp$y_up <- NA
    
    if (Lup >= l && Llow >= l) {
      if (!is.null(pred_up[[l]]) && !is.null(pred_low[[l]])) {
        df_tmp$y_low <- pred_low[[l]][idx]
        df_tmp$y_up <- pred_up[[l]][idx]
      }
    }
    df_plot <- rbind(df_plot, df_tmp)
  }
  
  # Reference data for true values
  df_ref <- data.frame(x_column = c(1:n_pred), y_column = yTrue[idx])
  df_ref$x_low <- df_ref$x_column + shift[1]
  df_ref$x_up <- df_ref$x_column + shift[n_meth]
  
  combined_plot <- ggplot() + theme_bw() +
    geom_point(data = df_plot, aes(x = x_column, y = y_column, color = factor(meth),
                                   shape = factor(meth), fill = factor(meth)), size = pt_sz) +
    geom_errorbar(data = df_plot, aes(x = x_column, y = y_column, ymin = y_low, ymax = y_up, color = factor(meth)),
                  width = 0, linewidth = br_lw, alpha = br_al) +
    
    geom_errorbar(data = df_plot, aes(x = x_column, y = y_low, ymin = y_low, ymax = y_low, color = factor(meth)),
                  width = br_wd, alpha = tc_al) +
    geom_errorbar(data = df_plot, aes(x = x_column, y = y_up, ymin = y_up, ymax = y_up, color = factor(meth)),
                  width = br_wd, alpha = tc_al) +
    
    geom_errorbar(data = df_ref, aes(x = x_column, y = y_column, ymin = y_column, ymax = y_column),
                  width = rf_sz, color = "black", alpha = 0.9) +
    
    theme(axis.text.x = element_blank(), legend.position = "bottom") +
    labs(x = "Observations", y = "Response") +
    scale_colour_manual(name = "Methods:", values = v_colors[1:n_meth]) +
    scale_shape_manual(name = "Methods:", values = v_pch[1:n_meth]) +
    scale_fill_manual(name = "Methods:", values = v_colors[1:n_meth]) +
    guides(colour = guide_legend(nrow = 1, override.aes = list(size = 1)))
  
  ggsave(paste(file.path(our_dir, file_name), '.pdf', sep = ''), combined_plot, height = 5, width = plt_w)
  
  print(combined_plot)
}
  
# Load data ----

#              BIP      coopL     IntegL     BSFP     JAFAR    JAFAR_T     JFR      JFR_T
new_col = c('#8D2FE5','#DC267F','#FE6100','#FFB000','#648FFF','#2DC4DC','#2B5B43','#57B38C')

n_split = 5

out_dir = '~/Desktop'

pow_temp <- readRDS('~/Downloads/jafar_ris_paper/StelzerEGA_cv-s_y_i-cusp_get-pow-temp.rds')
  
if(FALSE){
  
  ris_train_true <- ris_train_mean <- ris_train_low <- ris_train_up <- 
    ris_train_R2 <- ris_train_err2 <- ris_train_cover <- list()
  ris_test_true <- ris_test_mean <- ris_test_low <- ris_test_up <- 
    ris_test_R2 <- ris_test_err2 <- ris_test_cover <- list()
  
  for(ss in 1:n_split){
    
    print(paste0(" --- Loading Random split ", ss))
    
    # Data ----
      
    Data <- readRDS(paste0('~/Documents/GitHub/jafar/data/Sec4_StelzerEGA/StelzerEGA_cv-',ss,'_copula.rds'))
    
    train_true <- Data$yTrain
    test_true <- Data$yTest
    
    mean_y = Data$preprocess_y$mean
    std_y = Data$preprocess_y$std
    
    train_true = train_true*std_y + mean_y
    test_true = test_true*std_y + mean_y
    
    # CoopL, IntegL ----
    
    risCoopL  <- readRDS(paste0('~/Documents/GitHub/jafar/ris/StelzerEGA_cv-s',ss,'_y_CoopL_rep0.rds'))
    
    risIntegL <- readRDS(paste0('~/Documents/GitHub/jafar/ris/StelzerEGA_cv-s',ss,'_y_IntegL_rep0.rds'))
    
    sig_y_IntegL <- 1
    
    # jafar ----
    
    powT=0; 
    getR = colnames(pow_temp)[which(pow_temp[ss,]==powT)]
    risJAFAR <- readRDS(paste0('~/Downloads/jafar_ris_paper/StelzerEGA_cv-s',ss,'_y_i-cusp_nMC20000_nBurn15000_nThin10_rep',getR,'.rds'))
    if(risJAFAR$ris_MCMC$hyper_param$pow_tempering != powT){
      stop(paste0("Tempering power in risJAFAR not matching ",powT))
    }
    
    y_train_jafar <- risJAFAR$y_JAFAR_train$mean
    y_test_jafar <- risJAFAR$y_JAFAR_test$mean
    
    sig_y_jafar <- 1/sqrt(risJAFAR$ris_MCMC$s2_inv)
    
    rm(risJAFAR)
    
    # jafar_T ----
    
    powT=2/3; 
    getR = colnames(pow_temp)[which(pow_temp[ss,]==powT)]
    risJAFAR_T <- readRDS(paste0('~/Downloads/jafar_ris_paper/StelzerEGA_cv-s',ss,'_y_i-cusp_nMC20000_nBurn15000_nThin10_rep',getR,'.rds'))
    if(risJAFAR_T$ris_MCMC$hyper_param$pow_tempering != powT){
      stop(paste0("Tempering power in risJAFAR_T not matching ",powT))
    }
    
    y_train_jafar_T <- risJAFAR_T$y_JAFAR_train$mean
    y_test_jafar_T <- risJAFAR_T$y_JAFAR_test$mean
    
    sig_y_jafar_T <- 1/sqrt(risJAFAR_T$ris_MCMC$s2_inv)
    
    rm(risJAFAR_T)
    
    # jfr ----
    
    risJFR <- readRDS(paste0('~/Downloads/jafar_ris_paper/StelzerEGA_cv-s',ss,'_y_jfr_nMC20000_nBurn15000_nThin10_rep2.rds'))
    
    y_train_jfr <- risJFR$y_JAFAR_train$mean
    y_test_jfr <- risJFR$y_JAFAR_test$mean
    
    sig_y_jfr <- 1/sqrt(risJFR$ris_MCMC$s2_inv)
    
    rm(risJFR)
    
    # jfr_T ----
    
    risJFR_T <- readRDS(paste0('~/Downloads/jafar_ris_paper/StelzerEGA_cv-s',ss,'_y_jfr_nMC20000_nBurn15000_nThin10_rep0.rds'))
    
    y_train_jfr_T <- risJFR_T$y_JAFAR_train$mean
    y_test_jfr_T <- risJFR_T$y_JAFAR_test$mean
    
    sig_y_jfr_T <- 1/sqrt(risJFR_T$ris_MCMC$s2_inv)
    
    rm(risJFR_T)
    
    # bsfp ----
    
    risBSFP <- readRDS(paste0('~/Downloads/jafar_ris_paper/StelzerEGA_cv-s',ss,'_y_bsfp_nMC8000_nBurn4000_nThin10_rep0.rds'))
    
    y_train_bsfp = t(do.call(cbind,do.call(c,risBSFP$train_pred$EY.draw)))
    y_train_bsfp = y_train_bsfp[-c(1:floor(nrow(y_train_bsfp)/2)),]

    y_test_bsfp = t(do.call(cbind,do.call(c,risBSFP$test_pred$EY.draw)))
    y_test_bsfp = y_test_bsfp[-c(1:floor(nrow(y_test_bsfp)/2)),]
    
    sig_y_bsfp <- unlist(risBSFP$bsfp_mcmc$tau2.draw)
    sig_y_bsfp <- sqrt(tail(sig_y_bsfp,floor(length(sig_y_bsfp)/2)))
    
    rm(risBSFP)

    # Train Set Computations --------------------------------------------------------------------
    
    ## E[y] and 95% Intervals -----
    
    train_mean <- train_low <- train_up <- list() 
    
    train_mean[['jafar']] <- colMeans(y_train_jafar)
    train_uq <- credible_intervals_mcmc(y_train_jafar,sig_y_jafar)
    train_low[['jafar']]  <- train_uq[,1]
    train_up[['jafar']]   <- train_uq[,2]
    
    train_mean[['jafar_T']] <- colMeans(y_train_jafar_T)
    train_uq <- credible_intervals_mcmc(y_train_jafar_T,sig_y_jafar_T)
    train_low[['jafar_T']]  <- train_uq[,1]
    train_up[['jafar_T']]   <- train_uq[,2]
    
    train_mean[['jfr']] <- colMeans(y_train_jfr)
    train_uq <- credible_intervals_mcmc(y_train_jfr,sig_y_jfr)
    train_low[['jfr']]  <- train_uq[,1]
    train_up[['jfr']]   <- train_uq[,2]
    
    train_mean[['jfr_T']] <- colMeans(y_train_jfr_T)
    train_uq <- credible_intervals_mcmc(y_train_jfr_T,sig_y_jfr_T)
    train_low[['jfr_T']]  <- train_uq[,1]
    train_up[['jfr_T']]   <- train_uq[,2]
    
    train_mean[['bsfp']] <- colMeans(y_train_bsfp)
    train_uq <- credible_intervals_mcmc(y_train_bsfp,sig_y_bsfp)
    train_low[['bsfp']]  <- train_uq[,1]
    train_up[['bsfp']]   <- train_uq[,2]
    
    train_mean[['IntegL']] <- colMeans(t(risIntegL$train_pred))
    train_uq <- credible_intervals_mcmc(t(risIntegL$train_pred),1)
    train_low[['IntegL']]  <- train_uq[,1]
    train_up[['IntegL']]   <- train_uq[,2]
    
    train_mean[['CoopL']] <- risCoopL$y_pred_train
    
    method_names <- names(train_mean)
        
    # R2 
    train_R2 = sapply(train_mean,function(v) round(1-sum((Data$yTrain-v)^2)/sum((Data$yTrain-mean(Data$yTrain))^2),3))
    names(train_R2) <- method_names
    
    # Put Back Data to Original Scale
    
    train_mean = lapply(train_mean, function(v) v*std_y + mean_y)
    train_low = lapply(train_low, function(v) v*std_y + mean_y)
    train_up = lapply(train_up, function(v) v*std_y + mean_y)
    
    # Squared Errors
    
    train_err2 <- (train_true%*%t(rep(1,length(train_mean)))-do.call(cbind, train_mean))^2
    colnames(train_err2) <- names(train_mean)
    
    # Coverage of 95% Intervals
    
    train_cover <- (train_true%*%t(rep(1,length(train_mean)-1))<do.call(cbind, train_up))*
                   (train_true%*%t(rep(1,length(train_mean)-1))>do.call(cbind, train_low))
    colnames(train_cover) <- colnames(do.call(cbind, train_low))
    
    # Save Results
    
    ris_train_true[[ss]] <- train_true
    ris_train_mean[[ss]] <- train_mean
    ris_train_low[[ss]] <- train_low
    ris_train_up[[ss]] <- train_up
    ris_train_R2[[ss]] <- train_R2
    ris_train_err2[[ss]] <- train_err2
    ris_train_cover[[ss]] <- colMeans(train_cover)
    
    # Test Set Computations ---------------------------------------------------------------------
    
    ## E[y] and 95% Intervals -----
    
    test_mean <- test_low <- test_up <- list() 
    
    test_mean[['jafar']] <- colMeans(y_test_jafar)
    test_uq <- credible_intervals_mcmc(y_test_jafar,sig_y_jafar)
    test_low[['jafar']]  <- test_uq[,1]
    test_up[['jafar']]   <- test_uq[,2]
    
    test_mean[['jafar_T']] <- colMeans(y_test_jafar_T)
    test_uq <- credible_intervals_mcmc(y_test_jafar_T,sig_y_jafar_T)
    test_low[['jafar_T']]  <- test_uq[,1]
    test_up[['jafar_T']]   <- test_uq[,2]
    
    test_mean[['jfr']] <- colMeans(y_test_jfr)
    test_uq <- credible_intervals_mcmc(y_test_jfr,sig_y_jfr)
    test_low[['jfr']]  <- test_uq[,1]
    test_up[['jfr']]   <- test_uq[,2]
    
    test_mean[['jfr_T']] <- colMeans(y_test_jfr_T)
    test_uq <- credible_intervals_mcmc(y_test_jfr_T,sig_y_jfr_T)
    test_low[['jfr_T']]  <- test_uq[,1]
    test_up[['jfr_T']]   <- test_uq[,2]
    
    test_mean[['bsfp']] <- colMeans(y_test_bsfp)
    test_uq <- credible_intervals_mcmc(y_test_bsfp,sig_y_bsfp)
    test_low[['bsfp']]  <- test_uq[,1]
    test_up[['bsfp']]   <- test_uq[,2]
    
    test_mean[['IntegL']] <- colMeans(t(risIntegL$test_pred))
    test_uq <- credible_intervals_mcmc(t(risIntegL$test_pred),1)
    test_low[['IntegL']]  <- test_uq[,1]
    test_up[['IntegL']]   <- test_uq[,2]
    
    test_mean[['CoopL']] <- risCoopL$y_pred_test
    
    method_names <- names(test_mean)
    
    # R2 
    test_R2 = sapply(test_mean,function(v) round(1-sum((Data$yTest-v)^2)/sum((Data$yTest-mean(Data$yTest))^2),3))
    names(test_R2) <- method_names
    
    # Put Back Data to Original Scale
    
    test_mean = lapply(test_mean, function(v) v*std_y + mean_y)
    test_low = lapply(test_low, function(v) v*std_y + mean_y)
    test_up = lapply(test_up, function(v) v*std_y + mean_y)
    
    # Squared Errors
    
    test_err2 <- (test_true%*%t(rep(1,length(test_mean)))-do.call(cbind, test_mean))^2
    colnames(test_err2) <- names(test_mean)
    
    # Coverage of 95% Intervals
    
    test_cover <- (test_true%*%t(rep(1,length(test_mean)-1))<do.call(cbind, test_up))*
      (test_true%*%t(rep(1,length(test_mean)-1))>do.call(cbind, test_low))
    colnames(test_cover) <- colnames(do.call(cbind, test_low))
    
    # Save Results
    
    ris_test_true[[ss]] <- test_true
    ris_test_mean[[ss]] <- test_mean
    ris_test_low[[ss]] <- test_low
    ris_test_up[[ss]] <- test_up
    ris_test_R2[[ss]] <- test_R2
    ris_test_err2[[ss]] <- test_err2
    ris_test_cover[[ss]] <- colMeans(test_cover)
    
  }
  
  ris = list(train_true=ris_train_true,test_true=ris_test_true,
             train_mean=ris_train_mean,train_low=ris_train_low,
             train_up=ris_train_up,train_R2=ris_train_R2,train_err2=ris_train_err2,
             test_mean=ris_test_mean,test_low=ris_test_low,
             test_up=ris_test_up,test_R2=ris_test_R2,test_err2=ris_test_err2,
             train_cover=ris_train_cover,test_cover=ris_test_cover)
  
  saveRDS(ris,file.path(out_dir,paste0('StelzerEGA_predictions_ris.rds')))
   
} else {
  # ris <- readRDS(file.path(out_dir,paste0('StelzerEGA_predictions_ris.rds')))
  ris <- readRDS(file.path(out_dir,paste0('jafar_Sec4_predictions/StelzerEGA_predictions_ris.rds')))
}
  
if(F){
    
  for(ss in 1:n_split){
    
    print(paste0(" --- Loading Random split ", ss))
  
    # Train Set Plots --------------------------------------------------------------------
    
    ## predictions plot ----
    
    plot_title = paste0("s",ss,"_",'Train_set_Predictions')
    
    ggplot_predictions_NEW(ris$train_true[[ss]],ris$train_mean[[ss]],ris$train_low[[ss]],ris$train_up[[ss]],
                           method_names,plot_title,out_dir,
                           pt_sz=1.4,br_lw=1,br_wd=0.2,rf_sz=0.7,br_al=0.1,tc_al=0.5,sh_fct = 1.1,plt_w=14)
  
    ## SE boxplot -----------------------------------------------------------------
    
    plot_title = paste0("s",ss,"_",'Train_set_SquaredErrors')
    
    combined_plot <- reshape2::melt(ris$train_err2[[ss]]) %>%
      ggplot(aes(x=Var2, y=value, fill=Var2)) +
      geom_violin(width=1,color = scales::alpha('black',alpha=0.4),alpha=0.3) +
      geom_boxplot(width=0.3, color="black", alpha=0.8, outlier.size=0.7) +
      scale_fill_manual(values=c("#648FFF","#2DC4DC",'#2B5B43','#57B38C','#FFB000',"#FE6100","#DC267F")) +
      theme(aspect.ratio=1/2, legend.position="none") +
      labs(x = "", y = "Squared Errors") + 
      scale_y_continuous(trans='log10', labels = scales::label_scientific(),
                         limits = c(5*1e-5, 15), oob = scales::squish) +
      coord_fixed(ratio=0.25)
    
    ggsave(paste0(file.path(out_dir,plot_title),'_logScale','.pdf'),combined_plot,height=2.4,width=4.8)
    
    # Test Set Plots ---------------------------------------------------------------------
    
    ## predictions plot ----
    
    plot_title = paste0("s",ss,"_",'Test_set_Predictions')
    
    ggplot_predictions_NEW(ris$test_true[[ss]],ris$test_mean[[ss]],ris$test_low[[ss]],ris$test_up[[ss]],
                           method_names,plot_title,out_dir,
                           # pt_sz=0.8,rf_sz=0.3,br_lw=0.6,br_wd=1.1,br_al=0.3)
                           pt_sz=1.4,br_lw=2,br_wd=0.1,rf_sz=0.5,br_al=0.1,tc_al=0.5,sh_fct = 0.25,plt_w=7)
    
    ## SE boxplot -----
    
    plot_title = paste0("s",ss,"_",'Test_set_SquaredErrors')
    
    combined_plot <- reshape2::melt(ris$test_err2[[ss]]) %>%
      ggplot(aes(x=Var2, y=value, fill=Var2)) +
      geom_violin(width=1,color = scales::alpha('black',alpha=0.4),alpha=0.3) +
      geom_boxplot(width=0.3, color="black", alpha=0.8, outlier.size=0.7) +
      # scale_fill_manual(values=c("#DC267F","#FFB000","#785EF0","#004D40")) + 
      scale_fill_manual(values=c("#648FFF","#2DC4DC",'#2B5B43','#57B38C','#FFB000',"#FE6100","#DC267F")) + 
      theme(aspect.ratio=1/2, legend.position="none") +
      labs(x = "", y = "Squared Errors") + 
      scale_y_continuous(trans='log10', labels = scales::label_scientific(),
                         limits = c(5*1e-5, 15), oob = scales::squish) +
      coord_fixed(ratio=0.25)
    
    ggsave(paste0(file.path(out_dir,plot_title),'_logScale','.pdf'),combined_plot,height=2.4,width=4.8)
    
  }

  df_train<- bind_rows(lapply(1:n_split, function(k) {
    as.data.frame(ris$train_err2[[k]]) %>%
      pivot_longer(everything(), names_to = "Method", values_to = "errors") %>%
      mutate(split = paste("Split", k))
  }))
  
  combined_plot = ggplot(df_train, aes(x = split, y = errors, fill = Method)) +
    geom_boxplot(position = position_dodge(width = 0.8),
                 color="black", alpha=0.7, outlier.size=0.7,
                 outlier.colour = 'black',outlier.alpha = 0.3) +
    labs(title = NULL, x = NULL, y = "Squared Errors") +
    scale_fill_manual(values=c('#FFB000',"#DC267F","#FE6100","#648FFF","#2DC4DC",'#2B5B43','#57B38C')) + 
    scale_y_continuous(trans='log10', labels = scales::label_scientific(),
                       limits = c(5*1e-5, 15), oob = scales::squish) +
    theme(legend.position = "bottom") + guides(fill = guide_legend(nrow = 1))
  
  ggsave(paste(file.path('~/Desktop/', 'all_Train_set_SquaredErrors_logScale'), '.pdf', sep = ''),
         combined_plot, height = 4, width = 6)
  
  df_test <- bind_rows(lapply(1:n_split, function(k) {
    as.data.frame(ris$test_err2[[k]]) %>%
      pivot_longer(everything(), names_to = "Method", values_to = "errors") %>%
      mutate(split = paste("Split", k))
  }))
  
  combined_plot = ggplot(df_test, aes(x = split, y = errors, fill = Method)) +
    geom_boxplot(position = position_dodge(width = 0.8),
                 color="black", alpha=0.7, outlier.size=0.7,
                 outlier.colour = 'black',outlier.alpha = 0.3) +
    labs(title = NULL, x = NULL, y = "Squared Errors") +
    scale_fill_manual(values=c('#FFB000',"#DC267F","#FE6100","#648FFF","#2DC4DC",'#2B5B43','#57B38C')) + 
    scale_y_continuous(trans='log10', labels = scales::label_scientific(),
                       limits = c(5*1e-5, 15), oob = scales::squish) +
    theme(legend.position = "bottom") + guides(fill = guide_legend(nrow = 1))
  
  ggsave(paste(file.path('~/Desktop/', 'all_Test_set_SquaredErrors_logScale'), '.pdf', sep = ''),
         combined_plot, height = 4, width = 6)

}

#              BIP      coopL     IntegL     BSFP     JAFAR    JAFAR_T     JFR      JFR_T
new_col = c('#8D2FE5','#DC267F','#FE6100','#FFB000','#648FFF','#2DC4DC','#2B5B43','#57B38C')



if(F){
  # MSE
  colMeans(do.call(rbind,lapply(1:5, function(ss) colMeans((do.call(cbind,ris$train_mean[[ss]])-ris$train_true[[ss]][,1])^2))))
  colMeans(do.call(rbind,lapply(1:5, function(ss) colMeans((do.call(cbind,ris$test_mean[[ss]])-ris$test_true[[ss]][,1])^2))))

  # R2
  colMeans(do.call(rbind,ris$train_R2))
  colMeans(do.call(rbind,ris$test_R2))
   
  # 90% C.I. Emp Coverage
  colMeans(do.call(rbind,ris$train_cover))
  colMeans(do.call(rbind,ris$test_cover))
}

