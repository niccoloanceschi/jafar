
omics_colors = c(c('#2DC4DC','#648FFF','#8D2FE5','#DC267F','#FE6100','#FFB000','#57B38C','#2B5B43'),
                 rev(c('#DDCC77','#88CCEE','#679DCC','#BEB0D6','#CC6677','#882255')))

#' Plot MCMC samples of the inferred number of factors
#' 
#' @param risMCMC Posterior samples, output of \code{gibbs_jafar} or \code{gibbs_jfr}
#' @param out_path Output path where the generated plot will be saved
#' @param out_name_shared Output file name for the shared component plot (default: "n_factors_shared")
#' @param out_name_specific Output file name for the specific components plot (default: "n_factor_specific")
#' 
#' @export
#' 
plot_n_factors <- function(risMCMC,out_path='~/Desktop',
                           out_name_shared='n_factors_shared',
                           out_name_specific='n_factor_specific'){
  
  M = length(risMCMC$s2_inv_m)
  
  x_vals = seq(risMCMC$hyper_param$tBurnIn%/%risMCMC$hyper_param$tThin+1,
               risMCMC$hyper_param$tMCMC%/%risMCMC$hyper_param$tThin, by=1)
  
  # shared component
  
  K_max = max(risMCMC$K)
  
  pdf(file=file.path(out_path, paste0(out_name_shared,".pdf")), height=3, width=5)
  par(mar=c(3.7, 3.7, 1.2, 3.8), mgp=c(2,0.7,0))
  plot(risMCMC$K, type="l", lty=1, col=1, ylim=c(0, K_max),
       ylab="N. of Shared Factors", xlab="MCMC Iterations", main=NULL, lwd=1.3)
  matplot(x_vals, risMCMC$K_Lm_eff, type="l", lty=1, col=omics_colors[1:M],
          lwd=1.3, add=TRUE)
  par(xpd=FALSE)
  abline(v=risMCMC$hyper_param$tBurnIn %/% risMCMC$hyper_param$tThin,
         col="gray", lty=3, lwd=1.5)
  my_labs <- sapply(1:M, function(m) latex2exp::TeX(paste0("$\\Lambda_{", m, "}$")))
  par(xpd=TRUE)
  legend("right", legend=my_labs, pch=16, cex=0.8, ncol=1, inset=c(-0.17,0), col=omics_colors[1:M])
  dev.off()
  
  # specific components
  
  is_jafar = grepl('jafar',risMCMC$hyper_param$model)
  
  if(is_jafar){
    
    Km_max = max(unlist(risMCMC$K_Gm))
    
    pdf(file=file.path(out_path, paste0(out_name_specific,".pdf")), height=3, width=5)
    par(mar=c(3.7, 3.7, 1.2, 3.8), mgp=c(2,0.7,0))
    matplot(y=risMCMC$K_Gm,type='l',lty=1, col=omics_colors[1:M] ,ylim=c(0,Km_max),
            ylab='N. of Specific Factors',xlab='MCMC Iterations',main=NULL, lwd=1.3)
    par(xpd=FALSE)
    abline(v=risMCMC$hyper_param$tBurnIn %/% risMCMC$hyper_param$tThin,
           col="gray", lty=3, lwd=1.5)
    my_labs <- sapply(1:M, function(m) latex2exp::TeX(paste0("$\\Gamma_{", m, "}$")))
    par(xpd=TRUE)
    legend("right", legend=my_labs, pch=16, cex=0.8, ncol=1, inset=c(-0.17,0), col=omics_colors[1:M])
    dev.off()
  }
  
}


#' Plot the empirical and inferred within-view correlation matrices
#' 
#' @param risMCMC Posterior samples, output of \code{gibbs_jafar} or \code{gibbs_jfr}
#' @param X_m Training set multi-view predictors (optional, default: NULL).
#'    If NULL, only inferred correlation matrices are visualized.
#'    If not NULL, the empirical correlation matrices are displayed besides the inferred ones
#' @param out_path Output path where the generated plot will be saved (default: "~/Desktop/")
#' @param out_name Output file name (default: "correlations")
#'
#' @export
#'
plot_correlations <- function(risMCMC,X_m=NULL,out_path='~/Desktop/',
                              out_name='correlations'){
  
  M <- length(risMCMC$Cov_m_mean)
  
  plot_emp = !is.null(X_m)
  
  width0=4
  colors <- colorRampPalette(c("#AA4499", "white", "#117733"))(256)
  ticks <- seq(-1, 1, by = 0.5)
  
  png(file.path(out_path,paste0(out_name,".png")),
      width = (width0+1)*(1+plot_emp), height = width0*M, units = "in", res = 300)
  
  par(mfcol = c(M,1+plot_emp), mar = c(2,2,2,5), pty = "s") 
  
  for (m in 1:M) {
    # Inferred correlation
    cor_mat <- cov2cor(risMCMC$Cov_m_mean[[m]])
    image(cor_mat, useRaster=TRUE, asp=1, axes=FALSE, col=colors, zlim=c(-1,1))
    fields::image.plot(zlim=c(-1,1), col=colors, legend.only=TRUE,
                       smallplot = c(0.84, 0.88, 0.2, 0.8),
                       axis.args=list(at=ticks, labels=TRUE, cex.axis=0.9))
    mtext(bquote("Inferred" ~ cor(X[.(m)])), side = 3, line = 0.1, cex = 1)
    
    # Empirical correlation
    if(plot_emp) {
      cor_mat <- cor(X_m[[m]])
      image(cor_mat, useRaster=TRUE, asp=1, axes=FALSE, col=colors, zlim=c(-1,1))
      fields::image.plot(zlim=c(-1,1), col=colors, legend.only=TRUE,
                         smallplot = c(0.84, 0.88, 0.2, 0.8),
                         axis.args=list(at=ticks, labels=TRUE, cex.axis=0.9))
      mtext(bquote("Empirical" ~ cor(X[.(m)])), side = 3, line = 0.1, cex = 1)
    }
  }
  dev.off()
  
}


credible_intervals_mcmc <- function(m_samples, s_samples, Smc=10) {
  T_val <- nrow(m_samples)
  N_val <- ncol(m_samples)
  if(length(s_samples)==1){s_samples=rep(s_samples,T_val)}
  
  means_matrix_long <- m_samples[, rep(1:N_val, each = Smc)]
  sds_matrix_long <- matrix(rep(s_samples, each = N_val * Smc), nrow = T_val)
  
  noise_matrix <- matrix(rnorm(T_val * N_val * Smc), nrow = T_val)
  
  all_samples <- means_matrix_long + sds_matrix_long * noise_matrix
  
  all_samples_reshaped <- matrix(all_samples, ncol = N_val)
  
  credible_intervals <- apply(all_samples_reshaped, 2, quantile, probs = c(0.05, 0.95))
  
  t(credible_intervals)
}

unif_jitter <- function(yPred,range=0.25){
  n  = length(yPred)
  n1 = sum(yPred)
  n0 = n-n1
  
  delta_0 = 2*range/max(n0,n1)
  
  y_jitter <- yPred
  y_jitter[yPred==0] <- seq(-delta_0*floor(n0/2-0.5),delta_0*floor(n0/2),by=delta_0)
  y_jitter[yPred==1] <- 1+seq(-delta_0*floor(n1/2),delta_0*floor(n1/2-0.5),by=delta_0)
  
  y_jitter
}

#' Plot response predictions against true values
#' 
#' @param yPred Response predictions, output of \code{predict_y} or \code{predict_y_raw}
#' @param yTrue True values of the responses
#' @param risMCMC Posterior samples, output of \code{gibbs_jafar} or \code{gibbs_jfr}
#' @param out_path Output path where the generated plot will be saved (default: "~/Desktop/")
#' @param out_name Output file name (default: "predictions")
#'
#' @importFrom ggplot2 ggplot aes geom_line theme
#'
#' @export
#'
plot_predictions <- function(yPred,yTrue,risMCMC,out_path='~/Desktop/',out_name='predictions'){
  
  size=0.7
  linewidth=0.5
  width=0.7
  alpha=0.4
  range=0.25
  
  v_colors <- c("#57B38C")
  v_pch <- c(21)
  
  yPred = yPred$mean
  
  n_pred = ncol(yPred)
  
  yBinary = ('y_MC' %in% names(risMCMC))
  
  sig_y_jafar = sqrt(1/risMCMC$s2_inv)
  
  if(yBinary){
    yPred <- pnorm(yPred)
    sig_y_jafar = 0
  }
  
  df_plot = data.frame(y_mean = colMeans(yPred),
                       y_up = colMeans(yPred),
                       y_low = colMeans(yPred))
  
  train_uq <- credible_intervals_mcmc(yPred,sig_y_jafar)
  df_plot$y_low <- train_uq[,1]
  df_plot$y_up <- train_uq[,2]
  
  idx <- sort(df_plot$y_mean,index.return=T)$ix
  df_plot = df_plot[idx,]
  
  if(yBinary){
    df_plot$x_column = unif_jitter(yTrue[idx], range=range)
    width = 0.5*range/max(sum(yTrue),n_pred-sum(yTrue))
  } else {
    df_plot$x_column = c(1:n_pred)
    df_ref <- data.frame(x_column = c(1:n_pred), y_column = yTrue[idx])
  }
  
  combined_plot <- ggplot() + theme_bw() +
    geom_point(data = df_plot, aes(x = x_column, y = y_mean), 
               color=v_colors, fill=v_colors,shape=v_pch,size=size) +
    geom_errorbar(data = df_plot, aes(x=x_column, y = y_mean, ymin = y_low, ymax = y_up),
                  color=v_colors, width=width, linewidth=linewidth, alpha=alpha) +
    theme(aspect.ratio=0.7, legend.position="bottom") +
    guides(colour = guide_legend(override.aes = list(size=1)))
  
  if(yBinary){
    combined_plot <- combined_plot + labs(y = NULL, title = "P[ y = 1 | X = x ]", x = "True Response") +
      geom_hline(yintercept=0,col='gray',linewidth=.7) + geom_hline(yintercept=1,col='gray',linewidth=.7) +
      scale_x_continuous(breaks=c(0,1),labels=c(0,1)) + theme(plot.title = element_text(hjust = 0.5)) +
      scale_y_continuous(breaks=seq(0,1,by=0.1),labels=seq(0,1,by=0.1), limits=c(0, 1)) 
  } else { 
    combined_plot <- combined_plot + labs(y = NULL, title = "y | X = x", x = "Sorted Units") +
      geom_point(data = df_ref, aes(x = x_column, y = y_column), color = "black", size=0.) + 
      scale_x_continuous(breaks = NULL, labels=NULL) + theme(plot.title = element_text(hjust = 0.5)) +
      geom_errorbar(data = df_ref, aes(x=x_column, y = y_column, ymin = y_column, ymax = y_column),
                    width = size, color = "black",alpha=0.8)
  }
  
  outfile = paste0(out_name,'.pdf')
  ggsave(file.path(out_path,outfile), combined_plot, height=5*0.8, width=5)
}



#' Plot induced regression coefficients for y|X=x
#' 
#' @param yPred Response predictions, output of \code{predict_y} or \code{predict_y_raw}
#' @param out_path Output path where the generated plot will be saved (default: "~/Desktop/")
#' @param out_name Output file name (default: "coefficients")
#'
#' @export
#'
plot_coefficients <- function(yPred,out_path='~/Desktop/',out_name='coefficients'){
  
  ris_coeff = unlist(lapply(yPred$coeff,colMeans))
  
  M = length(yPred$coeff)
  p_m = sapply(yPred$coeff,ncol)
  
  coeff_color = c()
  for(m in 1:M){coeff_color = c(coeff_color,rep(omics_colors[m],p_m[m]))}
  
  eps = 0.1
  y_max = max(abs(ris_coeff))
  coeff_alpha = coeff_color
  for(j in 1:sum(p_m)){
    coeff_alpha[j] = adjustcolor(coeff_alpha[j], alpha.f = eps + (1-eps)*abs(ris_coeff[j])/y_max)
  }
  
  pdf(file=file.path(output_dir,paste0(out_name,'.pdf')),height=3,width=5)
  par(mar=c(3.2, 3.7, 1.5, 1.2), mgp=c(2,0.7,0), xpd=TRUE)
  plot(c(1:sum(p_m)),ris_coeff,col = coeff_color, ylim = 1.2*c(min(ris_coeff),max(ris_coeff)),
       xlab="",ylab='Regression Coefficients',main='',pch=20,cex=0.1,xaxt='n')
  arrows(c(1:sum(p_m)), pmax(rep(0,length(ris_coeff)),ris_coeff),
         c(1:sum(p_m)), pmin(rep(0,length(ris_coeff)),ris_coeff),
         length=0, angle=90, code=3, col = coeff_alpha)
  legend("bottom", inset=c(-0.2,-0.25), legend=paste0('View ',c(1:M)),
         ncol=3,pch=16, col=omics_colors[1:M], title=NULL,bty = "n")
  dev.off()
  
}

#'
#' @importFrom ggplot2 ggplot aes geom_line theme
#'
visualizeMatrix <- function(Sigma, texLabel = NULL, clim = NULL) {
  
  if (is.null(texLabel)) texLabel <- "$(\\Sigma_m)_{ij}$"
  
  xMin <- min(Sigma)
  xMax <- max(Sigma)
  if (!is.null(clim)) {
    xMin <- clim[1]
    xMax <- clim[2]
  }
  
  melted_cov <- reshape2::melt(Sigma)
  
  matPlot <- ggplot(data = melted_cov, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    labs(fill = latex2exp::TeX(texLabel)) +
    scale_fill_gradient2(low = "#AA4499", mid = "white", high = "#117733",
                         oob = scales::squish,
                         limits = c(xMin, xMax)) +
    theme_void() +
    theme(panel.background = element_rect(colour = "gray", fill = "gray", linewidth = 0.5))
  
  matPlot <- matPlot + theme(legend.position = "none") +
    labs(title = latex2exp::TeX(texLabel)) +
    theme(plot.title = element_text(size = 10,hjust = 0.5))
  
  return(matPlot)
}

#'
#' @importFrom ggplot2 ggplot aes geom_line theme
#'
colorBar_strip <- function(clim = c(-1, 1), texLabel = NULL) {
  
  if (is.null(texLabel)) texLabel <- "$\\Sigma$"
  
  ncolors <- 512
  
  z <- seq(clim[1], clim[2], length.out = ncolors)
  df <- data.frame(x = z, y = 1, z = z)
  
  blim = floor(clim[2] / 0.05) * 0.05
  breaks <- round(seq(-blim, blim, length.out = 5),2)
  
  stripPlot <- ggplot(df, aes(x, y, fill = z)) +
    geom_raster() +
    scale_fill_gradient2(limits = clim, low = "#AA4499", mid = "white", high = "#117733",
                         oob = scales::squish) +
    scale_x_continuous(breaks = breaks) +
    # coord_fixed(ratio = 1/15) +
    theme_minimal(base_size = 12) +
    theme(
      axis.title.y = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_blank(),
      panel.grid   = element_blank(),
      legend.position = "none",
      axis.text.x  = element_text(size = 9) 
    ) +
    labs(title = latex2exp::TeX(texLabel)) +
    theme(plot.title = element_text(size = 10,hjust = 0.5))
  
  return(stripPlot)
}

#'
#' @importFrom ggplot2 ggplot aes geom_line theme
#'
white_spacer <- function() {
  white_plot <- ggplot() + theme_void() +
    theme(panel.background = element_rect(fill = "white", colour = NA))
  return(white_plot)
}


#' Plot posterior means of factor loadings. 
#' 
#' @description Rotational alignment must be performed in advanced through the function \code{multiviewMatchAlign}
#'
#' @param risMCMC Posterior samples, output of \code{gibbs_jafar} or \code{gibbs_jfr}
#' @param out_path Output path where the generated plot will be saved
#' @param out_name_shared Output file name for the shared component plot (default: "n_factors_shared")
#' @param out_name_specific Output file name for the specific components plot (default: "n_factor_specific")
#' 
#' @importFrom ggplot2 ggplot aes geom_line theme
#'
#' @export
#'
plot_loadings <- function(risMCMC,out_path='~/Desktop/',
                          out_name_shared='shared_loadings',
                          out_name_specific='specific_loadings'){
  
  is_supervised = ('Theta'%in%names(risMCMC))
  is_jafar = grepl('jafar',risMCMC$hyper_param$model)
  
  M <- length(risMCMC$mu_m)
  p_m <- sapply(risMCMC$Lambda_m, ncol)
  
  # Compute posterior means
  
  L_mean_m <- lapply(risMCMC$Lambda_m, function(mat) apply(mat,c(2,3),mean))
  if(is_jafar){G_mean_m <- lapply(risMCMC$Gamma_m, function(mat) apply(mat,c(2,3),mean))} 
  
  if(is_supervised){
    T_mean = apply(risMCMC$Theta,2,mean)
    if(is_jafar){T_mean_m <- lapply(risMCMC$Theta_m, function(mat) apply(mat,2,mean))}
  }
  
  # Get limits
  
  K0 = ncol(L_mean_m[[1]])
  if(is_jafar){K0_m = sapply(G_mean_m,ncol)}
  
  lim_load = c(-1,1)*max(1,ceiling(100*max(abs(unlist(c(L_mean_m))))+1)/100)
  if(is_supervised){lim_th = c(-1,1)*ceiling(100*max(abs(c(T_mean)))+1)/100}
  
  if(is_jafar){
    lim_load = c(-1,1)*max(1,ceiling(100*max(abs(unlist(c(L_mean_m,G_mean_m))))+1)/100)
    if(is_supervised){lim_th = c(-1,1)*ceiling(100*max(abs(c(T_mean,unlist(T_mean_m))))+1)/100}
  }
  
  # Re-order columns
  
  idx_shared = order(rowMeans(sapply(L_mean_m, function(mat) colMeans(mat^2))),decreasing = T)
  idx_shared_m <- list()
  
  for(m in 1:M){
    L_mean_m[[m]] = L_mean_m[[m]][,idx_shared]
    if(is_jafar){
      idx_shared_m[[m]] <- order(colMeans(G_mean_m[[m]]^2),decreasing = T)
      G_mean_m[[m]] = G_mean_m[[m]][,idx_shared_m[[m]]]
    }
  }
  
  if(is_supervised){
    T_mean = T_mean[idx_shared]
    if(is_jafar){for(m in 1:M){T_mean_m[[m]] = T_mean_m[[m]][idx_shared_m[[m]]]}}
  }
  
  # plotting parameters
  
  Nrow <- 2
  rel_heights = c(12.5, 2)
  tot_height = 5
  tot_width = 2*M
  
  if(is_supervised){
    Nrow <- 3
    rel_heights = c(2, 13, 2)
    tot_height = 6
  }
  
  # plot shared component
  
  pL <- lapply(1:M, function(m) visualizeMatrix(t(L_mean_m[[m]]),
                                                texLabel = paste0("$\\Lambda_", m, "$"), clim = lim_load) + 
                 theme(plot.margin = margin(3, 5, 3, 5)))
  
  p0 <- replicate(M, white_spacer(), simplify = FALSE)
  p0[[M]] <- colorBar_strip(clim = lim_load, texLabel = paste0("$\\Lambda_m$"))
  
  all_plots <- c(pL,p0)
  
  if(is_supervised){
    p0[[1]] <- colorBar_strip(clim = lim_th, texLabel = paste0("$\\theta$"))
    
    pT <- visualizeMatrix(matrix(T_mean,ncol=1),
                          texLabel = paste0("$\\theta$"), clim = lim_th) + 
      theme(plot.margin = margin(0, 5, 0, 5))
    
    pT0 <- replicate(M, white_spacer(), simplify = FALSE)
    pT0[[(M+1)%/%2]] <- pT
    
    all_plots <- c(pT0,pL,p0)
  }
  
  combined_plots <- cowplot::plot_grid(plotlist=all_plots,nrow=Nrow,ncol=M,align="hv",rel_heights=rel_heights)
  
  ggsave(file.path(out_path, paste0(out_name_shared,".pdf")),
         combined_plots, width = tot_width, height = tot_height)
  
  # plot specific components
  
  if(is_jafar){
    
    pG <- lapply(1:M, function(m) visualizeMatrix(t(G_mean_m[[m]]),
                                                  texLabel = paste0("$\\Gamma_", m, "$"), clim = lim_load) + 
                   theme(plot.margin = margin(3, 5, 3, 5)))
    
    p0 <- replicate(M, white_spacer(), simplify = FALSE)
    p0[[M]] <- colorBar_strip(clim = lim_load, texLabel = paste0("$\\Gamma_m$"))
    
    all_plots <- c(pG,p0)
    
    if(is_supervised){
      p0[[1]] <- colorBar_strip(clim = lim_th, texLabel = paste0("$\\theta_m$"))
      
      pT <- lapply(1:M, function(m) visualizeMatrix(matrix(T_mean_m[[m]],ncol=1),
                                                    texLabel = paste0("$\\theta_", m,"$"), clim = lim_th) + 
                     theme(plot.margin = margin(0, 5, 0, 5)))
      
      all_plots <- c(pT,pG,p0)
    }
    
    combined_plots <- cowplot::plot_grid(plotlist=all_plots,nrow=Nrow,ncol=M,align="hv",rel_heights=rel_heights)
    
    ggsave(file.path(out_path, paste0(out_name_specific,".pdf")),
           combined_plots, width = tot_width, height = tot_height)
    
  }
  
}
