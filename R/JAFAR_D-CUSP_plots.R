
col_omics <- c('#332288','#882255','#CC6677','#8E021F','#648FFF','#785EF0','#8606D8','#DC267F')
col_theta <- '#117733'

#' Plot the of view-wise number of active shared factors throught the MCMC chain evolution
#' 
#' @param risMCMC Output of the Gibbs Sampler for JAFAR under the D-CUSP prior
#' @param out_folder Directory where to save output
#' 
#' @export
#' 
plot_n_fact_shared <- function(risMCMC,out_folder='~/Desktop'){
  
  M <- length(risMCMC$Lambda_m)
  
  active_J <- array(NA,dim(risMCMC$active_Lm)+c(0,0,1))
  active_J[,,-1] <- risMCMC$active_Lm
  active_J[,,1] <- risMCMC$active_T
  y_max <- max(apply(active_J,c(1,3),sum))+1
  
  pdf(file=file.path(out_folder,'n_active_shared.pdf'),height=3,width=5)
  par(mar=c(3.7, 3.7, 1.2,3.8), mgp=c(2,0.7,0), xpd=TRUE)
  p1 <- plot(rowSums(risMCMC$active_T),col=col_theta,type='l',ylim=c(0,y_max),
             xlab="MCMC Iterations",ylab='N. of Active Columns',lwd=1.4)
  for(m in 1:M){
    p1 <- p1 + lines(rowSums(risMCMC$active_Lm[,,m]),col=col_omics[m],lwd=1.4)
  }
  my_labs <- c(TeX('$\\Theta$'), sapply(1:M, function(m) TeX(paste0('$\\Lambda_{',m,'}$'))))
  
  legend('right', legend=my_labs, pch=16, cex=0.8, ncol=1,inset=c(-0.17,0), col=c(col_theta,col_omics[1:M]))
  
  dev.off()
}

#' Plot the inferred number of view-specific factors throught the MCMC chain evolution
#' 
#' @param risMCMC Output of the Gibbs Sampler for JAFAR under the D-CUSP prior
#' @param out_folder Directory where to save output
#' 
#' @export
#' 
plot_n_fact_specific <- function(risMCMC,out_folder='~/Desktop'){

  M <- length(risMCMC$Lambda_m)
  
  y_max <- max(risMCMC$K_Gm_eff)+1
  
  pdf(file=file.path(out_folder,'n_active_specific.pdf'),height=3,width=5)
  par(mar=c(3.7, 3.7, 1.2,3.8), mgp=c(2,0.7,0), xpd=TRUE)
  p1 <- plot(risMCMC$K_Gm_eff[,1],col=col_omics[1],type='l',ylim=c(0,y_max),
             xlab="MCMC Iterations",ylab='N. of Active Columns',lwd=1.4)
  if(M>1){
    for(m in 2:M){
      p1 <- p1 + lines(risMCMC$K_Gm_eff[,m],col=col_omics[m],lwd=1.4)
    }
  }
  
  my_labs <- sapply(1:M, function(m) TeX(paste0('$\\Gamma_{',m,'}$')))
  
  legend('right', legend=my_labs, pch=16, cex=0.8, ncol=1,inset=c(-0.17,0), col=col_omics[1:M])
  
  dev.off()
}

#' Plot the induced regression coefficients in the linear predictor of y|X=x
#' 
#' @param coef_jafar samples of the induced regression coefficients (output of \code{coeff_JAFAR})
#' @param out_folder Directory where to save output
#' 
#' @export
#' 
plot_coeff_jafar <- function(coef_jafar,out_folder='~/Desktop'){
  
  M = length(coef_jafar)
  p_m = sapply(coef_jafar,nrow)
  
  coef_jafar = unlist(lapply(coef_jafar,rowMeans))
  
  omics_color = unlist(lapply(1:M, function(m) rep(col_omics[m],p_m[m])))
  
  pdf(file=file.path(out_folder,'jafar_reg_coeff.pdf'),height=3,width=5)
  par(mar=c(3.2, 3.7, 1.5, 1.2), mgp=c(2,0.7,0), xpd=TRUE)
  plot(c(1:sum(p_m)),coef_jafar,col = omics_color, ylim = 1.2*c(min(coef_jafar),max(coef_jafar)),
       xlab="",ylab='Regression Coefficients',main='',pch=20,cex=0.1,xaxt='n')
  arrows(c(1:sum(p_m)), pmax(rep(0,length(coef_jafar)),coef_jafar),
         c(1:sum(p_m)), pmin(rep(0,length(coef_jafar)),coef_jafar),
         length=0, angle=90, code=3, col = omics_color)
  legend("bottom", inset=c(-0.2,-0.25), legend=sapply(1:M, function(m) paste0("m = ",m,"   ")),
         ncol=3,pch=16, col=col_omics[1:M], title=NULL,bty = "n")
  dev.off()
}

#' Plot the empirical and inferred within-viewa correlation matrices
#' 
#' @param risMCMC Output of the Gibbs Sampler for JAFAR under the D-CUSP prior
#' @param Z_m Train set multi-view predictors
#' @param out_folder Directory where to save output
#'
#' @export
#'
plot_cor_jafar <- function(risMCMC,Z_m,out_folder='~/Desktop/'){
  
  M <- length(risMCMC$Cov_m_mean)
  p_m <- sapply(Z_m, ncol) 
  
  n_colors=256
  
  for(which_cor in c('jafar','emp')){
    for(m in 1:M){
      
      if(which_cor=='jafar'){
        cor_mat = cov2cor(risMCMC$Cov_m_mean[[m]])
      } else {
        cor_mat = cor(Z_m[[m]])
      }
        
      lab_col=substitute(paste("cor(",X[v],")"), list(v=m))
      
      file_path <- file.path(out_folder, paste0(which_cor,'_cor_m',m,'.png'))
      png(file=file_path,height=5,width=5.6,res=300,pointsize=3.5,unit='cm')
      colors <- colorRampPalette(c("#AA4499", "white", "#117733"))(n_colors)
      ticks <- seq(-1, 1, by = 0.5)
      par(pty = "s", mar = c(1, 0, 1, 4) )  
      image(cor_mat, useRaster = TRUE,asp=1,axes=F,col = colors,
            xlim=c(0,1),ylim=c(0,1),zlim=c(-1, 1))
      image.plot(zlim = c(-1, 1), col = colors, legend.only = TRUE, side = 4,
                 axis.args = list(at = ticks, labels = TRUE,cex.axis = 0.9), legend.shrink = 0.5,
                 legend.width = 0.9,legend.mar = 4.5)
      mtext(lab_col, side = 4, line = 1, cex = 1.3, las=1, at=0.8)
      dev.off()
    }
  }
}

#' Plot the postprocessed shared-component loadings 
#' 
#' @param ris_PostProc Postprocessed shared component latent variables (output of \code{postprocess_JAFAR})
#' @param out_folder Directory where to save output
#'
#' @export
#'
plot_loadings_jafar <- function(ris_PostProc,out_folder='~/Desktop/'){
  
  K0 <- ncol(ris_PostProc$Theta)
  
  for(m in 1:M){
    
    mat_plot <- colMeans(ris_PostProc$Lambda_m[[m]])
    mat_plot <- mat_plot[,rev(1:K0)]
    
    colLabel <- TeX(paste0('$\\Lambda_{',m,'}$'))
      
    xMin = min(mat_plot)
    xMax = max(mat_plot)
    
    melted_mat <- melt(mat_plot)
    
    # ratio.display <- 1.5
    # ratio.values <- (max(melted_mat$Var2)-min(melted_mat$Var2))/(max(melted_mat$Var1)-min(melted_mat$Var1))
    
    covPlot <- ggplot(data = melted_mat, aes(x=Var2, y=Var1, fill=value)) +
      geom_tile() + labs(fill = colLabel) +
      scale_fill_gradient2(low="#AA4499", mid="white", high="#117733", oob = squish,
                           limits = c(xMin, xMax)) + scale_x_reverse() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black")) +
      theme(axis.title=element_blank(), axis.line=element_blank(),
            axis.text=element_blank(), axis.ticks=element_blank()) +
      theme(legend.key.width = unit(0.2, "cm"), legend.text=element_text(size=13)) +
      theme(legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,0,0,0)) +
      # coord_fixed(ratio.values / ratio.display) +
      # coord_fixed() + 
      theme(legend.title=element_text(size=14),
            aspect.ratio=1.9,
            panel.spacing = unit(-10, "cm"),
            panel.background = element_rect(colour = "gray", fill="gray", linewidth=0.5))
    
    ggsave(file.path(out_folder, paste0('jafar_Lambda_m',m,'.png')), covPlot, height=7.5, width=4.4)
  }
}

#' Deterministic jitter for plotting out-of-sample binary responses
#' 
#' @param vPred Predicted values
#' @param range Jitter's range
#'
unif_jitter <- function(vPred,range=0.25){
  n  = length(vPred)
  n1 = sum(vPred)
  n0 = n-n1
  
  delta_0 = 2*range/max(n0,n1)
  
  y_jitter <- vPred
  y_jitter[vPred==0] <- seq(-delta_0*floor(n0/2-0.5),delta_0*floor(n0/2),by=delta_0)
  y_jitter[vPred==1] <- 1+seq(-delta_0*floor(n1/2),delta_0*floor(n1/2-0.5),by=delta_0)
  
  y_jitter
}

#' Plot out-of-sample predictions of the response's linear predictors
#' 
#' @param vPred Predicted values
#' @param yTrue True values
#' @param s2y_inv MCMC samples of the response noise precision
#' @param is_binary Is the response binary?
#' @param out_folder Directory where to save output
#' @param size Point size
#' @param linewidth Linewidth of credible intervals bars
#' @param width Width of credible intervals extrema
#' @param alpha Opacity of credible intervals bars
#' @param range Jitter's range for binary response
#' @param filename Custom file name (optional)
#'
#' @export
#' 
plot_pred_jafar <- function(vPred,yTrue,s2y_inv=NULL,is_binary=F,out_folder='~/Desktop/',
                                  size=0.7,linewidth=0.5,width=0.7,alpha=0.4,range=0.25,filename=NULL){
  
  ref_size=size
  
  v_colors <- c("#FFB000")
  v_pch <- c(21)
  
  n_pred = nrow(vPred)
  
  if(is_binary){vPred <- pnorm(vPred)}
  
  df_plot = data.frame(y_mean = rowMeans(vPred),
                       y_up = apply(vPred,1,quantile,0.975),
                       y_low = apply(vPred,1,quantile,0.025))
  
  idx <- sort(df_plot$y_mean,index.return=T)$ix
  df_plot = df_plot[idx,]
  
  if(is_binary){
    df_plot$x_column = unif_jitter(yTrue[idx], range=range)
    width = 0.5*range/max(sum(yTrue),n_pred-sum(yTrue))
  } else {
    df_plot$x_column = c(1:n_pred)
    
    sig_y_jafar = mean(1/sqrt(s2y_inv))
    df_plot$y_up = df_plot$y_up + 1.96*sig_y_jafar
    df_plot$y_low = df_plot$y_low - 1.96*sig_y_jafar
    
    df_ref <- data.frame(x_column = c(1:n_pred), y_column = yTrue[idx])
  }
  
  combined_plot <- ggplot() + theme_bw() +
    geom_point(data = df_plot, aes(x = x_column, y = y_mean), 
               color=v_colors, fill=v_colors,shape=v_pch,size=size) +
    geom_errorbar(data = df_plot, aes(x=x_column, y = y_mean, ymin = y_low, ymax = y_up),
                  color=v_colors, width=width, linewidth=linewidth, alpha=alpha) +
    theme(aspect.ratio=0.7, legend.position="bottom") +
    guides(colour = guide_legend(override.aes = list(size=1)))
  
  if(is_binary){
    combined_plot <- combined_plot + labs(y = "P[ y=1 | X=x ]", x = "True Response") +
      geom_hline(yintercept=0,col='gray',linewidth=.7) + geom_hline(yintercept=1,col='gray',linewidth=.7) +
      scale_x_continuous(breaks=c(0,1),labels=c(0,1)) +
      scale_y_continuous(breaks=seq(0,1,by=0.1),labels=seq(0,1,by=0.1), limits=c(0, 1)) 
  } else { 
    combined_plot <- combined_plot + labs(y = "y | X=x", x = NULL) +
      geom_point(data = df_ref, aes(x = x_column, y = y_column), color = "black", size=0.) + 
      scale_x_continuous(breaks = NULL, labels=NULL) + # theme(axis.text.x = element_blank()) +
      geom_errorbar(data = df_ref, aes(x=x_column, y = y_column, ymin = y_column, ymax = y_column),
                      width = ref_size, color = "black",alpha=0.8)
  }
  
  outfile = paste0('jafar_y_pred.pdf')
    if(!is.null(filename)){outfile=paste0('jafar_y_pred_',filename,'.pdf')}
  ggsave(file.path(out_folder,outfile), combined_plot, height=5*0.8, width=5)
}







