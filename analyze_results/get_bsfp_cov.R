
rm(list = ls())

library(stringr)

run_cluster = T

if(run_cluster){
  folder_name = '/hpc/group/herringlab/na224/multiomics/ris/sim_sec3'
  if (length(commandArgs(trailingOnly = TRUE)) != 2) {
    stop("Usage: your_script.R <n> <s>")
  }
  nn <- as.integer(commandArgs(trailingOnly = TRUE)[2])
  ss <- as.integer(commandArgs(trailingOnly = TRUE)[1])
} else {
  folder_name = '/Users/nico/Downloads/'
  nn <- 50
  ss <- 10 
}

all_rds_files <- list.files(path = folder_name, pattern = "\\.rds$", full.names = FALSE, ignore.case = TRUE)

file_list <- all_rds_files[grepl("bsfp", all_rds_files, ignore.case = TRUE) &
                           grepl(paste0("_n",nn,"_"), all_rds_files, ignore.case = TRUE) &
                           grepl(paste0("_s",ss,"_"), all_rds_files, ignore.case = TRUE) ]

num_files <- length(file_list)

cat("Found", num_files, "files matching criteria.\n")
cat("Starting processing...\n")

for (i in seq_along(file_list)) { # Use seq_along to get index 'i'
  
  file_name <- file_list[i] # Get the current file name using the index
  
  cat(sprintf("File %d out of %d: %s\n", i, num_files, file_name))
  
  risBSFP <- readRDS(file.path(folder_name,file_name))
  
  elem_rem = c("Xm_train_bsfp", "Xm_test_bsfp", "tXm_train_bsfp", "tXm_test_bsfp")
  
  if(any(elem_rem %in% names(risBSFP))){
    
    nMC <- as.numeric(str_replace(str_extract(file_name, "nMC(\\d+)"), "nMC", ""))
    nBurn <- as.numeric(str_replace(str_extract(file_name, "nBurn(\\d+)"), "nBurn", ""))
    nThin <- as.numeric(str_replace(str_extract(file_name, "nThin(\\d+)"), "nThin", ""))
    
    t_MC = round(nMC/nThin)
    t_Burn = round(nBurn/nThin)
    
    iter_print <- (t_MC-t_Burn) %/% 10
    
    M = nrow(risBSFP$bsfp_mcmc$U.draw[[1]])
    p_m = sapply(1:M, function(m) nrow(risBSFP$bsfp_mcmc$U.draw[[1]][[m,1]]))
    
    Cov_m_mean <- lapply(1:M, function(m) matrix(0,p_m[m],p_m[m]))
    Marg_Var_m <- lapply(1:M, function(m) matrix(NA,t_MC-t_Burn,p_m[m]))
    
    for(t in 1:(t_MC-t_Burn)){
      for(m in 1:M){
        Cov_m_mean[[m]] <- Cov_m_mean[[m]] + cov(t(risBSFP$Xm_train_bsfp[[m]][t_Burn+t-1,,])) / (t_MC-t_Burn)
        Marg_Var_m[[m]][t,] <- apply(t(risBSFP$Xm_train_bsfp[[m]][t_Burn+t-1,,]),2,var)
      }
      if(t %% iter_print == 0){
        print(sprintf(fmt = "%10s%3s%2s", "[",(t%/%iter_print)*10,"%]"))
      }
    }
    
    for(m in 1:M){Cov_m_mean[[m]] <- Cov_m_mean[[m]]*(risBSFP$bsfp_mcmc$sigma.mat[m,1]^2)}
    
    for(el in elem_rem){risBSFP[[el]] = NULL}
    
    risBSFP$Cov_m_mean = Cov_m_mean
    risBSFP$Marg_Var_m = Marg_Var_m
    
    saveRDS(risBSFP, file.path(folder_name,file_name))
  } else {
    cat(' --- Skipping files, must have been already processed')
  }
}



