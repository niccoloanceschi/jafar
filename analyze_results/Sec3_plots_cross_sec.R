

rm(list = ls())

library(fields)
library(dplyr)
library(tidyr)
library(ggplot2)
library(multiview)
library(knitr)
library(kableExtra)

save_cor_png <- function(cor_m,folder='~/Desktop/',filename='prova',n_colors=256,lab_col=NULL){
  
  if(is.null(lab_col)){
    lab_col=cor('X')
  }
  
  file_path <- file.path(folder, paste(filename,'.png',sep=""))
  png(file=file_path,height=5,width=5.6,res=300,pointsize=3.5,unit='cm')
  colors <- colorRampPalette(c("#AA4499", "white", "#117733"))(n_colors)
  ticks <- seq(-1, 1, by = 0.5)
  par(pty = "s", mar = c(1, 0, 1, 4) )  
  image(cor_m, useRaster = TRUE,asp=1,axes=F,col = colors,
        xlim=c(0,1),ylim=c(0,1),zlim=c(-1, 1))
  image.plot(zlim = c(-1, 1), col = colors, legend.only = TRUE, side = 4,
             axis.args = list(at = ticks, labels = TRUE,cex.axis = 0.9), legend.shrink = 0.5,
             legend.width = 0.9,legend.mar = 4.5)
  mtext(lab_col, side = 4, line = 1, cex = 1.3, las=1, at=0.8)
  dev.off()
}

# graphical paramters ----

mth_names = c('BIP',   'CoopL',  'IntegL', 'BSFP',   'JAFAR', 'JAFAR_T',  'JFR',    'JFR_T')
v_colors = c('#8D2FE5','#DC267F','#FE6100','#FFB000','#648FFF','#2DC4DC','#2B5B43','#57B38C')
v_pch <- c(8, 22, 23, 21, 24, 25, 15, 17)
v_pch <- c(8, 24, 25, 21, 22, 27, 23, 9)

v_meth = c('JAFAR','JFR','BSFP','BIP')


data_path = '~/Documents/GitHub/jafar/data/Sec3_Simulations/'
ris_pathJAFAR = '~/Downloads/jafar_ris_paper/sim_sec3_J/'
ris_pathJFR = '~/Downloads/jafar_ris_paper/sim_sec3_jfr/'
ris_pathBSFP = '~/Downloads/jafar_ris_paper/sim_sec3_bsfp/'
ris_pathBIP = '~/Downloads/jafar_ris_paper/sim_sec3_bip/'

# LOAD RESULTS -----------

source('~/Documents/GitHub/jafar/new_source/Sec3_simulations_preproc.R') # data  
source('~/Documents/GitHub/jafar/new_source/jafar_predictions.R') # data  

nn = 200
ss = 56

n0=200

print(paste0(' --- Loading Results and Data for n=',nn,' s=',ss, ' --- '))

dataFile <- paste0('Simulated_data_n',n0,'_s',ss,'.rds')
Data <- readRDS(file.path(data_path, dataFile))

if(TRUE){
  # Data ----
  Data <- data_scale_subset(Data,nn)
}
    
for(m in 1:Data$M){  
  cor_m <- cov2cor(tcrossprod(Data$Lambda_m[[m]])+tcrossprod(Data$Gamma_m[[m]])+diag(Data$s2m[[m]]))

  save_cor_png(cor_m,filename = paste("cor_true_m",m,"_",sep=''),
               lab_col = substitute(paste("cor(",X[v],")"), list(v=m)))
}

v_col <- c('#332288','#882255','#CC6677','#8E021F')

omics_color = c()
for(m in 1:Data$M){omics_color = c(omics_color,rep(v_col[m],Data$p_m[m]))}

K_0 = ncol(Data$eta)
Km_0 = sapply(Data$phi_m,ncol)

pred_true <- jafar_coeff_y_t(Data$X_m,Data$n,Data$M,Data$p_m,K_0,Km_0,
                             lapply(1:Data$M, function(m) (1/Data$preprocess_X_m[[m]]$std)*Data$Lambda_m[[m]]),
                             lapply(1:Data$M, function(m) (1/Data$preprocess_X_m[[m]]$std)*Data$Gamma_m[[m]]),
                             (1/Data$preprocess_y$std)*Data$Theta,
                             lapply(Km_0,function(kk) rep(0,kk)),
                             mean(Data$yTrain),
                             Data$preprocess_y$std^2/Data$s2y, 
                             lapply(Data$X_m, colMeans),
                             lapply(1:Data$M, function(m) Data$preprocess_X_m[[m]]$std^2/Data$s2m[[m]]),
                             rescale_pred=F)


coef_true = unlist(pred_true$coeff)

pdf(file='~/Desktop/coeff_true.pdf',height=3,width=5)
par(mar=c(3.2, 3.7, 1.5, 1.2), mgp=c(2,0.7,0), xpd=TRUE)
# plot(c(1:sum(Data$p_m)),coef_true,col = omics_color, ylim = 1.2*c(min(coef_true),max(coef_true)),
plot(c(1:sum(Data$p_m)),coef_true,col = omics_color, ylim =c(-0.075,0.06),
     xlab="",ylab='Regression Coefficients',main='',pch=20,cex=0.1,xaxt='n')
arrows(c(1:sum(Data$p_m)), pmax(rep(0,length(coef_true)),coef_true),
       c(1:sum(Data$p_m)), pmin(rep(0,length(coef_true)),coef_true),
       length=0, angle=90, code=3, col = omics_color)
legend("bottom", inset=c(-0.2,-0.25), legend=c("m = 1   ","m = 2   ","m = 3   "),
       ncol=3,pch=16, col=v_col, title=NULL,bty = "n")
dev.off()

# Inference -----

if(T){
      
  ## BIP ------
            
  which_rep=0
  risFile <- paste0('Sec3_Simulated_data_n',nn,'_s',ss,'_y_bip_nMC6000_nBurn3000_rep',which_rep,'.rds')
  risBIP <- tryCatch(readRDS(file.path(ris_pathBIP, risFile)), 
                     error = function(ee) NULL, warning= function(ww) NULL)
  
  if (!is.null(risBIP)) {
     
    for(m in 1:Data$M){
      cor_m <- cov2cor(risBIP$cov_mean_m[[m]])
      
      save_cor_png(cor_m,filename = paste("cor_bip_m",m,"_",sep=''),
                   lab_col = substitute(paste("cor(",X[v],")"), list(v=m)))
    }
    
    rm(risBIP)
  } else {
    stop(paste0('BIP ris NOT available for n=',nn,' s=',ss))
  }
          
  # BSFP ------
  
  which_rep=0
  risFile <- paste0('Sec3_Simulated_data_n',nn,'_s',ss,'_y_bsfp_nMC6000_nBurn3000_nThin10_rep',which_rep,'.rds')
  risBSFP <- tryCatch(readRDS(file.path(ris_pathBSFP, risFile)), 
                      error = function(ee) NULL, warning= function(ww) NULL)
  
  if (!is.null(risBSFP)) {
    
    for(m in 1:Data$M){
      cor_m <- cov2cor(risBSFP$Cov_m_mean[[m]])
      
      save_cor_png(cor_m,filename = paste("cor_bsfp_m",m,"_",sep=''),
                   lab_col = substitute(paste("cor(",X[v],")"), list(v=m)))
    }
    
    rm(risBSFP)
  } else {
    stop(paste0('BSFP ris NOT available for n=',nn,' s=',ss))
  }
          
  # JAFAR ------
  
  which_rep=0 # 
  risFile <- paste0('Sec3_Simulated_data_n',nn,'_s',ss,'_y_i-cusp_nMC20000_nBurn15000_nThin10_rep',which_rep,'.rds')
  risJAFAR <- tryCatch(readRDS(file.path(ris_pathJAFAR, risFile)), 
                       error = function(ee) NULL, warning= function(ww) NULL)
  if (!is.null(risJAFAR)) {
    
    y_JAFAR_train <- jafar_coeff_y(Data$X_m,risJAFAR$ris_MCMC,rescale_pred=T)
    
    coef_jafar = unlist(lapply(y_JAFAR_train$coeff,colMeans))
    
    pdf(file='~/Desktop/coeff_jafar.pdf',height=3,width=5)
    par(mar=c(3.2, 3.7, 1.5, 1.2), mgp=c(2,0.7,0), xpd=TRUE)
    # plot(c(1:sum(Data$p_m)),coef_true,col = omics_color, ylim = 1.2*c(min(coef_true),max(coef_true)),
    plot(c(1:sum(Data$p_m)),coef_jafar,col = omics_color, ylim =c(-0.075,0.06),
         xlab="",ylab='Regression Coefficients',main='',pch=20,cex=0.1,xaxt='n')
    arrows(c(1:sum(Data$p_m)), pmax(rep(0,length(coef_jafar)),coef_jafar),
           c(1:sum(Data$p_m)), pmin(rep(0,length(coef_jafar)),coef_jafar),
           length=0, angle=90, code=3, col = omics_color)
    legend("bottom", inset=c(-0.2,-0.25), legend=c("m = 1   ","m = 2   ","m = 3   "),
           ncol=3,pch=16, col=v_col, title=NULL,bty = "n")
    dev.off()
    
    for(m in 1:Data$M){
      cor_m <- cov2cor(risJAFAR$ris_MCMC$Cov_m_mean[[m]])
      
      save_cor_png(cor_m,filename = paste("cor_jafar_m",m,"_",sep=''),
                   lab_col = substitute(paste("cor(",X[v],")"), list(v=m)))
    }
      
    # rm(risJAFAR)
  } else {
    stop(paste0('JAFAR ris NOT available for n=',nn,' s=',ss))
  }
          
  # JFR ------
        
  which_rep=0 # 
  risFile <- paste0('Sec3_Simulated_data_n',nn,'_s',ss,'_y_jfr_nMC20000_nBurn15000_nThin10_rep',which_rep,'.rds')
  risJFR <- tryCatch(readRDS(file.path(ris_pathJFR, risFile)), 
                     error = function(ee) NULL, warning= function(ww) NULL)
  if (!is.null(risJFR)) {
    
    y_JFR_train <- jafar_coeff_y(Data$X_m,risJFR$ris_MCMC,rescale_pred=T)
    
    coef_jfr = unlist(lapply(y_JFR_train$coeff,colMeans))
    
    pdf(file='~/Desktop/coeff_jfr.pdf',height=3,width=5)
    par(mar=c(3.2, 3.7, 1.5, 1.2), mgp=c(2,0.7,0), xpd=TRUE)
    # plot(c(1:sum(Data$p_m)),coef_true,col = omics_color, ylim = 1.2*c(min(coef_true),max(coef_true)),
    plot(c(1:sum(Data$p_m)),coef_jfr,col = omics_color, ylim =c(-0.075,0.06),
         xlab="",ylab='Regression Coefficients',main='',pch=20,cex=0.1,xaxt='n')
    arrows(c(1:sum(Data$p_m)), pmax(rep(0,length(coef_jfr)),coef_jfr),
           c(1:sum(Data$p_m)), pmin(rep(0,length(coef_jfr)),coef_jfr),
           length=0, angle=90, code=3, col = omics_color)
    legend("bottom", inset=c(-0.2,-0.25), legend=c("m = 1   ","m = 2   ","m = 3   "),
           ncol=3,pch=16, col=v_col, title=NULL,bty = "n")
    dev.off()
    
    for(m in 1:Data$M){
      cor_m <- cov2cor(risJFR$ris_MCMC$Cov_m_mean[[m]])
      
      save_cor_png(cor_m,filename = paste("cor_jfr_m",m,"_",sep=''),
                   lab_col = substitute(paste("cor(",X[v],")"), list(v=m)))
    }
    
    # rm(risJFR)
  } else {
    stop(paste0('JFR ris NOT available for n=',nn,' s=',ss))
  }
}









