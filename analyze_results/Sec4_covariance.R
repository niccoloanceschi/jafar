
rm(list = ls())
# graphics.off()

# Functions --------------------------------------------------------------------

library(fields)
library(ggplot2)

save_cor_png <- function(cor_m,folder='~/Desktop/',filename='prova',n_colors=256,lab_col=NULL){
  
  if(is.null(lab_col)){
    lab_col=cor('X')
    # lab_col=substitute(paste("cor(",X[v],")"), list(v=mPlot))
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
  abline(a=0, b=1, col=adjustcolor("#117733", alpha.f = 0.5), lwd=0.3) 
  dev.off()
}

pow_temp <- readRDS('~/Downloads/jafar_ris_paper/StelzerEGA__y_i-cusp_pow_temp.rds')

rel_diff = array(NA,c(5,3,2))

for(ss in 1:5){
  
  # MCMC Results & Data ---------

  if(T){ 

    Data <- readRDS(paste0('~/Documents/GitHub/jafar/data/Sec4_StelzerEGA/StelzerEGA_cv-',ss,'_copula.rds'))
  
    powT=0; 
    getR = colnames(pow_temp)[which(pow_temp[ss,]==powT)]
    risJAFAR <- readRDS(paste0('~/Downloads/jafar_ris_paper/StelzerEGA_cv-s',ss,'_y_i-cusp_nMC20000_nBurn15000_nThin10_rep',getR,'.rds'))
    
    powT=2/3; 
    getR = colnames(pow_temp)[which(pow_temp[ss,]==powT)]
    risJAFAR_T <- readRDS(paste0('~/Downloads/jafar_ris_paper/StelzerEGA_cv-s',ss,'_y_i-cusp_nMC20000_nBurn15000_nThin10_rep',getR,'.rds'))
  }

  ## covariances ----

  for(mPlot in 1:Data$M){
    
    title_m = paste0("$cor(X_",mPlot,")$")
    
    # Empirical --
    cor_emp_m <- cor(matrix(matrix(c(Data$X_m[[mPlot]]),ncol=1),ncol=Data$p_m[mPlot]),use = 'pairwise.complete.obs')
    save_cor_png(cor_emp_m,filename = paste("s",ss,"_cor_Emp_m",mPlot,"_",sep=''),
                 lab_col = substitute(paste("cor(",X[v],")"), list(v=mPlot)))
    
    # JAFAR --
    cor_jar_m = cov2cor(risJAFAR$ris_MCMC$Cov_m_mean[[mPlot]])
    save_cor_png(cor_jar_m,filename = paste("s",ss,"_cor_JAFAR_m",mPlot,"_",sep=''),
                 lab_col = substitute(paste("cor(",X[v],")"), list(v=mPlot)))
    
    # JAFAR TEMPERING --
    cor_jar_T_m = cov2cor(risJAFAR_T$ris_MCMC$Cov_m_mean[[mPlot]])
    save_cor_png(cor_jar_T_m,filename = paste("s",ss,"_cor_JAFAR_T_m",mPlot,"_",sep=''),
                 lab_col = substitute(paste("cor(",X[v],")"), list(v=mPlot)))
    
    rel_diff[ss,mPlot,1] = sum((cor_jar_m-cor_emp_m)^2)/sum(cor_emp_m^2)
    rel_diff[ss,mPlot,2] = sum((cor_jar_T_m-cor_emp_m)^2)/sum(cor_emp_m^2)
  }
  
  ## n. of factors ----
  M=Data$M
  
  pdf(paste0('~/Desktop/s',ss,'_nfact_shared.pdf'),width=4.5,height=3)
  par(mar=c(4,4,1,1))
  plot(risJAFAR$ris_MCMC$K,type='l',ylim=c(0,20),cex.main = 0.8,
       ylab="n. of factors",xlab='mcmc iter',main='Shared Component');
  matplot(risJAFAR$ris_MCMC$K_Lm_eff,type='l',lty=1,col=1+1:M,add=T)
  matplot(risJAFAR$ris_MCMC$K_T_eff,type='l',lty=2,col=1,add=T)
  dev.off()
  
  pdf(paste0('~/Desktop/s',ss,'_nfact_specif.pdf'),width=4.5,height=3)
  par(mar=c(4,4,1,1))
  matplot(risJAFAR$ris_MCMC$K_Gm,type='l',lty=1,col=1+1:M,ylim=c(0,30),cex.main = 0.8,
          ylab="n. of factors",xlab='mcmc iter',main='Specific Components')
  matplot(risJAFAR$ris_MCMC$K_Tm_eff,type='l',lty=2,col=1+1:M,ylim=c(0,30),add=T)
  dev.off()
  
  pdf(paste0('~/Desktop/s',ss,'_nfact_shared_T.pdf'),width=4.5,height=3)
  par(mar=c(4,4,1,1))
  plot(risJAFAR_T$ris_MCMC$K,type='l',ylim=c(0,20),cex.main = 0.8,
       ylab="n. of factors",xlab='mcmc iter',main='Shared Component');
  matplot(risJAFAR_T$ris_MCMC$K_Lm_eff,type='l',lty=1,col=1+1:M,add=T)
  matplot(risJAFAR_T$ris_MCMC$K_T_eff,type='l',lty=3,col=1,add=T)
  dev.off()
  
  pdf(paste0('~/Desktop/s',ss,'_nfact_specif_T.pdf'),width=4.5,height=3)
  par(mar=c(4,4,1,1))
  matplot(risJAFAR_T$ris_MCMC$K_Gm,type='l',lty=1,col=1+1:M,ylim=c(0,30),cex.main = 0.8,
          ylab="n. of factors",xlab='mcmc iter',main='Specific Components')
  matplot(risJAFAR_T$ris_MCMC$K_Tm_eff,type='l',lty=3,col=1+1:M,ylim=c(0,30),add=T)
  dev.off()
  
}

dimnames(rel_diff) = list(paste0('s',c(1:5)),paste0('M',c(1:M)),c('jafar','jafar_T'))

saveRDS(rel_diff,'~/Desktop/StelzerEGA_jafar_rel_frob_nrom.RDS')

apply(rel_diff,c(2,3),median)
