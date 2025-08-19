
rm(list = ls())

latex_table_from_matrices <- function(mean_mat, sd_mat, caption = "Table", label = "tab:results") {
  if (!all(dim(mean_mat) == dim(sd_mat))) {
    stop("Both matrices must have the same dimensions.")
  }
  
  # Combine mean and sd into formatted entries
  formatted <- matrix(
    # paste0(sprintf("%.1f",round(mean_mat,1)), "$_{", sprintf("%.1f",round(sd_mat,1)), "}$"),
    paste0(sprintf("%.4f",round(mean_mat,4)), "$_{", sprintf("%.4f",round(sd_mat,4)), "}$"),
    nrow = nrow(mean_mat),
    ncol = ncol(mean_mat),
    dimnames = dimnames(mean_mat)
  )
  
  # Build LaTeX code
  header <- paste0(" & ", paste(colnames(formatted), collapse = " & "), " \\\\ \\midrule")
  rows <- apply(formatted, 1, function(r) paste0(names(r)[1], " & ", paste(r, collapse = " & "), " \\\\"))
  
  latex_code <- c(
    "\\begin{table}[ht]",
    "\\centering",
    "\\begin{tabular}{l", paste(rep("c", ncol(mean_mat)), collapse = ""), "}",
    "\\toprule",
    header,
    rows,
    "\\bottomrule",
    "\\end{tabular}",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"),
    "\\end{table}"
  )
  
  return(paste(latex_code, collapse = "\n"))
}


data_path = '~/Documents/GitHub/jafar/data/Sec2_Simulations/'
ris_pathNAIVE = '~/Downloads/jafar_ris_paper/sim_sec2_Jnaive/'
ris_pathJAFAR = '~/Downloads/jafar_ris_paper/sim_sec2_J/'
ris_pathJFR = '~/Downloads/jafar_ris_paper/sim_sec2_jfr/'

# cycle ----

n0=50

M=3

s_max = 100 # 
n_meth = 3

if(T){
  
  mm_tot = 0
  for(m in 1:M){
    for(mm in 1:M){
      if(mm>m){
        mm_tot =  mm_tot + 1
      }
    }
  }
  
  ris_K <- array(NA,c(s_max,n_meth,1+2*M+2))
  ris_frob_norm <- array(NA,c(s_max,n_meth,3*M))
  ris_frob_norm_intra <- array(NA,c(s_max,n_meth,mm_tot))
  
  for(ss in 1:s_max){
    
    ## load ----
    
    print(paste0(' --- Loading Results and Data for n=',n0,' s=',ss, ' --- '))
    
    dataFile <- paste0('Simulated_data_n',n0,'_s',ss,'.rds')
    Data <- readRDS(file.path(data_path, dataFile))
    
    filename = paste0('Sec2_Simulated_data_n50_s',ss,'_i-cusp-naive_nMC10000_nBurn5000_nThin10_rep0.rds')
    risNAIVE <- readRDS(file.path(ris_pathNAIVE,filename))
    
    filename = paste0('Sec2_Simulated_data_n50_s',ss,'_i-cusp_nMC10000_nBurn5000_nThin10_rep0.rds')
    risJAFAR <- readRDS(file.path(ris_pathJAFAR,filename))
    
    filename = paste0('Sec2_Simulated_data_n50_s',ss,'_jfr_nMC10000_nBurn5000_nThin10_rep0.rds')
    risJFR <- readRDS(file.path(ris_pathJFR,filename))
    
    # n fact ----
  
    K_jafar = round(c(mean(risJAFAR$ris_MCMC$K)-1,
                        colMeans(risJAFAR$ris_MCMC$K_Lm_eff),
                        colMeans(risJAFAR$ris_MCMC$K_Gm_eff)),1)
    
    K_naive = round(c(mean(risNAIVE$ris_MCMC$K)-1,
                        colMeans(risNAIVE$ris_MCMC$K_Lm_eff),
                        colMeans(risNAIVE$ris_MCMC$K_Gm_eff)),1)
    
    K_jfr = round(c(mean(risJFR$ris_MCMC$K)-1,
                            colMeans(risJFR$ris_MCMC$K_Lm_eff),
                            rep(0,M)),1)
    
    tab_K = rbind(K_jafar,K_naive,K_jfr) 
    colnames(tab_K) <- c('KJ',paste0('KJ-',1:M),paste0('KI-',1:M))
    
    tab_01 <- round(rbind(
    c(mean(apply(risJAFAR$ris_MCMC$active_Lm,1,function(mat) sum(rowSums(mat)==0)))-1,
      mean(apply(risJAFAR$ris_MCMC$active_Lm,1,function(mat) sum(rowSums(mat)==1)))) ,
    c(mean(apply(risNAIVE$ris_MCMC$active_Lm,1,function(mat) sum(rowSums(mat)==0)))-1,
      mean(apply(risNAIVE$ris_MCMC$active_Lm,1,function(mat) sum(rowSums(mat)==1)))) , 
    c(mean(apply(risJFR$ris_MCMC$active_Lm,1,function(mat) sum(rowSums(mat)==0)))-1,
      mean(apply(risJFR$ris_MCMC$active_Lm,1,function(mat) sum(rowSums(mat)==1)))) ),1)
    
    rownames(tab_01) <- c('jafar','naive','jfr')
    colnames(tab_01) <- c('All-zero-cols','One-only-cols')
    
    # corr tot ----
    
    cor_true_m <- cor_m_J <- cor_m_N <- cor_m_jfr <- list()
    
    frob_norm = matrix(0,nrow=3,ncol=Data$M)
      
    for(m in 1:Data$M){
      cor_true_m[[m]] <- (tcrossprod(Data$Lambda_m[[m]])+tcrossprod(Data$Gamma_m[[m]])+diag(Data$s2m[[m]]))
      # cor_true_m[[m]] <- cov2cor(cor_true_m[[m]])
      
      cor_m_J[[m]] <- (risJAFAR$ris_MCMC$Cov_m_mean[[m]])
      frob_norm[1,m] <- sum((cor_true_m[[m]]-cov2cor(cor_m_J[[m]]))^2) / (Data$p_m[m])^2
      
      cor_m_N[[m]] <- (risNAIVE$ris_MCMC$Cov_m_mean[[m]])
      frob_norm[2,m] <- sum((cor_true_m[[m]]-cov2cor(cor_m_N[[m]]))^2) / (Data$p_m[m])^2
      
      cor_m_jfr[[m]] <- (risJFR$ris_MCMC$Cov_m_mean[[m]])
      frob_norm[3,m] <- sum((cor_true_m[[m]]-cov2cor(cor_m_jfr[[m]]))^2) / (Data$p_m[m])^2
    
    }
  
    # corr ----
    
    nMC = length(risJAFAR$ris_MCMC$K)
    
    cor_L_true_m <- cor_L_m_J <- cor_L_m_N <- cor_L_m_jfr <- list()
    cor_G_true_m <- cor_G_m_J <- cor_G_m_N <- cor_G_m_jfr <- list()
    
    cor_LG_true_m <- cor_LG_m_J <- cor_LG_m_N <- cor_LG_m_jfr <- list()
    
    frob_norm_L = matrix(0,nrow=3,ncol=Data$M)
    frob_norm_G = matrix(0,nrow=3,ncol=Data$M)
    
    frob_norm_LG = matrix(0,nrow=3,ncol=mm_tot)
    
    for(m in 1:Data$M){
      
      cor_L_true_m[[m]] <- (tcrossprod(Data$Lambda_m[[m]]))
      cor_G_true_m[[m]] <- (tcrossprod(Data$Gamma_m[[m]]))
      
      for(mm in 1:Data$M){
        if(mm>m){
          m_idx = paste0(m,mm)
          cor_LG_true_m[[m_idx]] <- (tcrossprod(Data$Lambda_m[[m]],Data$Lambda_m[[mm]]))

          cor_LG_m_J[[m_idx]] <- cor_LG_m_N[[m_idx]] <- cor_LG_m_jfr[[m_idx]] <- matrix(0,Data$p_m[m],Data$p_m[mm])
        }
      }
      
      cor_L_m_J[[m]] <- cor_L_m_N[[m]] <- matrix(0,Data$p_m[m],Data$p_m[m])
      cor_G_m_J[[m]] <- cor_G_m_N[[m]] <- matrix(0,Data$p_m[m],Data$p_m[m])
    }
    
    for(t in 1:nMC){
    
      ## JAFAR ----
      
      s2_m <- La_m <- Ga_m <- list()
      
      for(m in 1:Data$M){
        
        s2_m[[m]] = risJAFAR$ris_MCMC$s2_inv_m[[m]][t,]
        La_m[[m]] = risJAFAR$ris_MCMC$Lambda_m[[m]][t,,]
        Ga_m[[m]] = risJAFAR$ris_MCMC$Gamma_m[[m]][t,,]
        
        mar_std_m = sqrt(1/s2_m[[m]] + rowSums(La_m[[m]]^2) + rowSums(Ga_m[[m]]^2))
        
        s2_m[[m]] = s2_m[[m]]*(mar_std_m^2)
        La_m[[m]] = La_m[[m]]/mar_std_m
        Ga_m[[m]] = Ga_m[[m]]/mar_std_m
        
        cor_L_m_J[[m]] <- cor_L_m_J[[m]] + tcrossprod(La_m[[m]]) / nMC
        cor_G_m_J[[m]] <- cor_G_m_J[[m]] + tcrossprod(Ga_m[[m]]) / nMC
      }
        
      for(m in 1:Data$M){
        for(mm in 1:Data$M){
          if(mm>m){
            m_idx = paste0(m,mm)
            cor_LG_m_J[[m_idx]] <- cor_LG_m_J[[m_idx]] + tcrossprod(La_m[[m]],La_m[[mm]]) / nMC
          }
        }
      }
        
      # NAIVE ----
      
      s2_m <- La_m <- Ga_m <- list()
      
      for(m in 1:Data$M){
        s2_m[[m]] = risNAIVE$ris_MCMC$s2_inv_m[[m]][t,]
        La_m[[m]] = risNAIVE$ris_MCMC$Lambda_m[[m]][t,,]
        Ga_m[[m]] = risNAIVE$ris_MCMC$Gamma_m[[m]][t,,]
        
        mar_std_m = sqrt(1/s2_m[[m]] + rowSums(La_m[[m]]^2) + rowSums(Ga_m[[m]]^2))
        
        s2_m[[m]] = s2_m[[m]]*(mar_std_m^2)
        La_m[[m]] = La_m[[m]]/mar_std_m
        Ga_m[[m]] = Ga_m[[m]]/mar_std_m
        
        cor_L_m_N[[m]] <- cor_L_m_N[[m]] + tcrossprod(La_m[[m]]) / nMC
        cor_G_m_N[[m]] <- cor_G_m_N[[m]] + tcrossprod(Ga_m[[m]]) / nMC
      }
      
      for(m in 1:Data$M){  
        for(mm in 1:Data$M){
          if(mm>m){
            m_idx = paste0(m,mm)
            cor_LG_m_N[[m_idx]] <- cor_LG_m_N[[m_idx]] + tcrossprod(La_m[[m]],La_m[[mm]]) / nMC
          }
        }
      }
        
      # JFR ----
      
      s2_m <- La_m <- Ga_m <- list()
        
      for(m in 1:Data$M){
        s2_m[[m]] = risJFR$ris_MCMC$s2_inv_m[[m]][t,]
        La_m[[m]] = risJFR$ris_MCMC$Lambda_m[[m]][t,,]
        
        mar_std_m = sqrt(1/s2_m[[m]] + rowSums(La_m[[m]]^2) )
        
        s2_m[[m]] = s2_m[[m]]*(mar_std_m^2)
        La_m[[m]] = La_m[[m]]/mar_std_m
      }
      
      for(m in 1:Data$M){  
        for(mm in 1:Data$M){
          if(mm>m){
            m_idx = paste0(m,mm)
            cor_LG_m_jfr[[m_idx]] <- cor_LG_m_jfr[[m_idx]] + tcrossprod(La_m[[m]],La_m[[mm]]) / nMC
          }
        }
        
      }
    }
    
    # dif ----
    
    m_count = 0
    for(m in 1:Data$M){
      
      frob_norm_L[1,m] <- sum((cor_L_true_m[[m]]-(cor_L_m_J[[m]]))^2) / (Data$p_m[m])^2  
      frob_norm_G[1,m] <- sum((cor_G_true_m[[m]]-(cor_G_m_J[[m]]))^2) / (Data$p_m[m])^2
      
      frob_norm_L[2,m] <- sum((cor_L_true_m[[m]]-(cor_L_m_N[[m]]))^2) / (Data$p_m[m])^2
      frob_norm_G[2,m] <- sum((cor_G_true_m[[m]]-(cor_G_m_N[[m]]))^2) / (Data$p_m[m])^2
      
      for(mm in 1:Data$M){
        if(mm>m){
          m_idx = paste0(m,mm)
          m_count =  m_count + 1

          frob_norm_LG[1,m_count] <- sum((cor_LG_true_m[[m_idx]]-(cor_LG_m_J[[m_idx]]))^2) / (Data$p_m[m]*Data$p_m[mm])
          frob_norm_LG[2,m_count] <- sum((cor_LG_true_m[[m_idx]]-(cor_LG_m_N[[m_idx]]))^2) / (Data$p_m[m]*Data$p_m[mm])
          frob_norm_LG[3,m_count] <- sum((cor_LG_true_m[[m_idx]]-(cor_LG_m_jfr[[m_idx]]))^2) / (Data$p_m[m]*Data$p_m[mm])
        }
      }
      
    }
    
    # store ----
    
    ris_K[ss,,] <- cbind(tab_K,tab_01)
    ris_frob_norm[ss,,] <- cbind(frob_norm,frob_norm_L,frob_norm_G)
    ris_frob_norm_intra[ss,,] <- frob_norm_LG
  }
  
  # save ----
  
  saveRDS(list(ris_K=ris_K,ris_frob_norm=ris_frob_norm,
               ris_frob_norm_intra=ris_frob_norm_intra),
          '~/Desktop/ris_Sec2_mod.rds')
} else {

  ris_all <- readRDS('~/Desktop/ris_Sec2_mod.rds')
  # ris_all <- readRDS('~/Desktop/ris_Sec2_all.rds')
  
  ris_K=ris_all[['ris_K']]
  ris_frob_norm=ris_all[['ris_frob_norm']]
  ris_frob_norm_intra=ris_all[['ris_frob_norm_intra']]
  
}

# analysis ----

## n fact ----

ris_K_mean <- round(apply(ris_K,c(2,3),mean),1)
ris_K_sd <- round(apply(ris_K,c(2,3),sd),1)

rownames(ris_K_mean) <- c('jafar','naive','jfr')
colnames(ris_K_mean) <- c('KJ',paste0('KJ-',1:M),paste0('KI-',1:M),'All-zero','One-only')

rownames(ris_K_sd) <- rownames(ris_K_mean)
colnames(ris_K_sd) <- colnames(ris_K_mean)

ris_K_mean[,c(1:4,8,9,5:7)]
ris_K_sd[,c(1:4,8,9,5:7)]

cat(latex_table_from_matrices(ris_K_mean[,c(5:7,1:4,8:9)],
                 ris_K_sd[,c(5:7,1:4,8:9)]))

## corr ----

ris_frob_norm_mean = round(apply(ris_frob_norm,c(2,3),mean),4)
ris_frob_norm_sd = round(apply(ris_frob_norm,c(2,3),sd),4)

colnames(ris_frob_norm_mean) <- c(paste0('cor_',1:M),paste0('cor_L_',1:M),paste0('cor_G_',1:M))

ris_frob_norm_intra_mean = round(apply(ris_frob_norm_intra,c(2,3),mean),4)
ris_frob_norm_intra_sd = round(apply(ris_frob_norm_intra,c(2,3),sd),4)

colnames(ris_frob_norm_intra_mean) <- c('cor_12','cor_13','cor_23')

cat(latex_table_from_matrices(
  cbind(ris_frob_norm_mean[,1:3],ris_frob_norm_intra_mean),
  cbind(ris_frob_norm_sd[,1:3],ris_frob_norm_intra_sd) ))

cat(latex_table_from_matrices(ris_frob_norm_mean[,4:9],ris_frob_norm_sd[,4:9]))
# others ----

if(F){
  # unique shared patterns (NO zero enforcement)
  which_comb0 <- c('1010','1010') 
  value_comb0 <- c(1,1)
  
  # unique shared patterns (zero enforcement)
  which_comb <- c('1100','0110')
  value_comb <- c(1,1)
  
  patterns_matrix <- do.call(rbind, lapply(c(which_comb0, which_comb), function(p) as.integer(strsplit(p, "")[[1]])))
  
  n_active_factors <- colSums(patterns_matrix * c(value_comb0, value_comb))[1:Data$M]
  
  K_true = c(ncol(Data$eta),n_active_factors,sapply(Data$phi,ncol))
}
  


