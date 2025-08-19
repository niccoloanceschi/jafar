
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

source('~/Documents/GitHub/jafar/new_source/Sec3_simulations_preproc.R') # data  
source('~/Documents/GitHub/jafar/new_source/jafar_predictions.R') # data  

# load ----

ss = 5
n0=50

print(paste0(' --- Loading Results and Data for n=',n0,' s=',ss, ' --- '))

data_path = '~/Documents/GitHub/jafar/data/Sec2_Simulations/'
dataFile <- paste0('Simulated_data_n',n0,'_s',ss,'.rds')
Data <- readRDS(file.path(data_path, dataFile))

ris_pathNAIVE = '~/Downloads/jafar_ris_paper/sim_sec2_Jnaive/'
ris_pathJAFAR = '~/Downloads/jafar_ris_paper/sim_sec2_J/'
ris_pathJFR = '~/Downloads/jafar_ris_paper/sim_sec2_jfr/'

filename = paste0('Sec2_Simulated_data_n50_s',ss,'_i-cusp-naive_nMC10000_nBurn5000_nThin10_rep0.rds')
risNAIVE <- readRDS(file.path(ris_pathNAIVE,filename))

filename = paste0('Sec2_Simulated_data_n50_s',ss,'_i-cusp_nMC10000_nBurn5000_nThin10_rep0.rds')
risJAFAR <- readRDS(file.path(ris_pathJAFAR,filename))

filename = paste0('Sec2_Simulated_data_n50_s',ss,'_jfr_nMC10000_nBurn5000_nThin10_rep0.rds')
risJFR <- readRDS(file.path(ris_pathJFR,filename))

M=Data$M

# true fact ----

# unique shared patterns (NO zero enforcement)
which_comb0 <- c('1010','1010') 
value_comb0 <- c(1,1)

# unique shared patterns (zero enforcement)
which_comb <- c('1100','0110')
value_comb <- c(1,1)

patterns_matrix <- do.call(rbind, lapply(c(which_comb0, which_comb), function(p) as.integer(strsplit(p, "")[[1]])))

n_active_factors <- colSums(patterns_matrix * c(value_comb0, value_comb))[1:Data$M]

K_true = c(ncol(Data$eta),n_active_factors,sapply(Data$phi,ncol))

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

tab_K = rbind(K_true,K_jafar,K_naive,K_jfr) 
colnames(tab_K) <- c('KJ',paste0('KJ-',1:M),paste0('KI-',1:M))

tab_K

tab_01 <- round(rbind(
  c(mean(apply(risJAFAR$ris_MCMC$active_Lm,1,function(mat) sum(rowSums(mat)==0)))-1,
    mean(apply(risJAFAR$ris_MCMC$active_Lm,1,function(mat) sum(rowSums(mat)==1)))) ,
  c(mean(apply(risNAIVE$ris_MCMC$active_Lm,1,function(mat) sum(rowSums(mat)==0)))-1,
    mean(apply(risNAIVE$ris_MCMC$active_Lm,1,function(mat) sum(rowSums(mat)==1)))) , 
  c(mean(apply(risJFR$ris_MCMC$active_Lm,1,function(mat) sum(rowSums(mat)==0)))-1,
    mean(apply(risJFR$ris_MCMC$active_Lm,1,function(mat) sum(rowSums(mat)==1)))) ),1)

rownames(tab_01) <- c('jafar','naive','jfr')
colnames(tab_01) <- c('All-zero-cols','One-only-cols')

tab_01

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

frob_norm


# corr ----

nMC = length(risJAFAR$ris_MCMC$K)

cor_L_true_m <- cor_L_m_J <- cor_L_m_N <- cor_L_m_jfr <- list()

frob_norm_L = matrix(0,nrow=2,ncol=Data$M)

cor_G_true_m <- cor_G_m_J <- cor_G_m_N <- cor_G_m_jfr <- list()

frob_norm_G = matrix(0,nrow=2,ncol=Data$M)

for(m in 1:Data$M){
  
  cor_L_true_m[[m]] <- (tcrossprod(Data$Lambda_m[[m]]))
  cor_G_true_m[[m]] <- (tcrossprod(Data$Gamma_m[[m]]))
  
  cor_L_m_J[[m]] <- cor_L_m_N[[m]] <- matrix(0,Data$p_m[m],Data$p_m[m])
  cor_G_m_J[[m]] <- cor_G_m_N[[m]] <- matrix(0,Data$p_m[m],Data$p_m[m])
  
  for(t in 1:nMC){
    
    # JAFAR
    
    s2_m = risJAFAR$ris_MCMC$s2_inv_m[[m]][t,]
    La_m = risJAFAR$ris_MCMC$Lambda_m[[m]][t,,]
    Ga_m = risJAFAR$ris_MCMC$Gamma_m[[m]][t,,]
    
    mar_std_m = sqrt(1/s2_m + rowSums(La_m^2) + rowSums(Ga_m^2))
    
    s2_m = s2_m*(mar_std_m^2)
    La_m = La_m/mar_std_m
    Ga_m = Ga_m/mar_std_m
    
    cor_L_m_J[[m]] <- cor_L_m_J[[m]] + tcrossprod(La_m) / nMC
    cor_G_m_J[[m]] <- cor_G_m_J[[m]] + tcrossprod(Ga_m) / nMC
    
    # NIAVE
    
    s2_m = risNAIVE$ris_MCMC$s2_inv_m[[m]][t,]
    La_m = risNAIVE$ris_MCMC$Lambda_m[[m]][t,,]
    Ga_m = risNAIVE$ris_MCMC$Gamma_m[[m]][t,,]
    
    mar_std_m = sqrt(1/s2_m + rowSums(La_m^2) + rowSums(Ga_m^2))
    
    s2_m = s2_m*(mar_std_m^2)
    La_m = La_m/mar_std_m
    Ga_m = Ga_m/mar_std_m
    
    cor_L_m_N[[m]] <- cor_L_m_N[[m]] + tcrossprod(La_m) / nMC
    cor_G_m_N[[m]] <- cor_G_m_N[[m]] + tcrossprod(Ga_m) / nMC
    
  }
  
  frob_norm_L[1,m] <- sum((cor_L_true_m[[m]]-(cor_L_m_J[[m]]))^2) / (Data$p_m[m])^2  
  frob_norm_G[1,m] <- sum((cor_G_true_m[[m]]-(cor_G_m_J[[m]]))^2) / (Data$p_m[m])^2
  
  frob_norm_L[2,m] <- sum((cor_L_true_m[[m]]-(cor_L_m_N[[m]]))^2) / (Data$p_m[m])^2
  frob_norm_G[2,m] <- sum((cor_G_true_m[[m]]-(cor_G_m_N[[m]]))^2) / (Data$p_m[m])^2
  
}

frob_norm_L
frob_norm_G






