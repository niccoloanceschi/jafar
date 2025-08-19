
rm(list = ls())

run_cluster = F

if(run_cluster){
  # path <- '/hpc/group/dunsonlab/na224/multiomics/' # cluster simulations
  path <- '/hpc/group/herringlab/na224/multiomics/' # cluster simulations
} else{ 
  path <- '~/Documents/GitHub/jafar/' # local simulations
}
setwd(path)

# Loading Packages -------------------------------------------------------------

# library(...)

# Files and Folder -------------------------------------------------------------

## Select model and data ----

run_simulations = F
run_supervised = F
predict_views = F

whichCUSP = 'jfr' # 'i-cusp-naive' # 'jfr' # 'i-cusp-naive' # 
if(!whichCUSP%in%c('d-cusp','i-cusp','i-cusp-naive', 'jfr')){
  stop("whichCUSP must be on of 'd-cusp', 'i-cusp', 'i-cusp-naive' or 'jfr'")
}

## Source code ----

source_dir = "new_source"
if(!dir.exists(source_dir)){stop("Source folder not found", call.=FALSE)}

source_files = c("jafar_gibbs.R", "jfr_gibbs.R", "jafar_initialization.R",
                 "jafar_updates.R", "jafar_predictions.R",
                 "Sec3_simulations_preproc.R")

for(s_file in source_files){
  if(!file.exists(file.path(source_dir, s_file))){stop(paste0("Source file '",s_file,"' not found"))}
  source(file.path(source_dir, s_file))
}

source_files_cpp = c("jafar_updates.cpp", "jafar_updates_parallel.cpp")

for(s_file in source_files_cpp){
  if(!file.exists(file.path(source_dir, s_file))){stop(paste0("cpp source file '",s_file,"' not found"))}
  Rcpp::sourceCpp(file.path(source_dir, s_file))
}

## Output Folder ----

out_dir = 'ris/'

if(!dir.exists(out_dir)){dir.create(out_dir)}

## Data import  ----

print(' | Loading Data ')

if(run_simulations){
  
  # out_dir <- paste0(out_dir,'sim_sec3_J')
  out_dir <- paste0(out_dir,'sim_sec2_J')
  if(!dir.exists(out_dir)){dir.create(out_dir)}
  
  if(run_cluster){
    if (length(commandArgs(trailingOnly = TRUE)) != 2) {
      stop("Usage: your_script.R <n> <s>")
    }
    nn <- as.integer(commandArgs(trailingOnly = TRUE)[2])
    ss <- as.integer(commandArgs(trailingOnly = TRUE)[1])
  } else {
    nn <- 50 # 100 # 200 # 
    ss <- 2 # 10 # 5 # 
  }
  # data_file <- paste0('Simulated_data_n200_s',ss)
  data_file <- paste0('Simulated_data_n50_s',ss)
} else {
  if(run_cluster){
    if (length(commandArgs(trailingOnly = TRUE)) < 1) {
      stop("Usage: your_script.R <s>")
    }
    ss <- as.integer(commandArgs(trailingOnly = TRUE)[1])
  } else {
    ss <- 2
  }
  data_file <- 'StelzerEGA_cv-'
}

if(run_simulations){
  # data_dir = "data/Sec3_Simulations"
  data_dir = "data/Sec2_Simulations"
  data_path = file.path(data_dir, paste0(data_file,'.rds'))
} else {
  data_dir = "data/Sec4_StelzerEGA"
  data_path = file.path(data_dir, paste0(data_file,ss,'_copula.rds'))
}
if(!dir.exists(data_dir)){stop("Input data folder not found", call.=FALSE)}
if(!file.exists(data_path)){stop("Data file not found")}

Data <- readRDS(data_path)

if(run_simulations){
  if(grepl('Sec3',data_dir)){
    Data <- data_scale_subset(Data,nn)
    data_file <- paste0('Sec3_Simulated_data_n',nn,'_s',ss)
  } else {
    data_file <- paste0('Sec2_',data_file)
  }
} else {
  data_file <- paste0(data_file,'s',ss)
}

y_run=NULL
binary_y=F
if(run_supervised){
  y_run = Data$yTrain
  binary_y=as.logical(min(Data$yTrain%%1==0))
}

# TODO: remove from here
for(m in 1:Data$M){
  set.seed(123)
  # Data$X_m[[m]] = Data$X_m[[m]][,1:(m*100)]
  # Data$X_m[[m]] = Data$X_m[[m]][,c(floor(Data$p_m[m]/2)+1:(m*100))]
  Data$X_m[[m]] = Data$X_m[[m]][,sample(1:Data$p_m[m],m*100)]
  # Data$X_m[[m]] = Data$X_m[[m]][,c(1:(m*50),rev(1:Data$p_m[m])[1:(m*50)])]
  # Data$X_m[[m]] = Data$X_m[[m]][,rev(1:Data$p_m[m])[1:(m*100)]]
  # Data$X_m[[m]] = Data$X_m[[m]][,rev(1:Data$p_m[m])[100+1:(m*100)]]
}
Data$p_m = sapply(Data$X_m, ncol)
# TODO: remove until here

# MCMC Parameters --------------------------------------------------------------

K0 = 40 # 25 # 20 #
K0_m = c(20,20,20) # c(25,35,25) #
if(whichCUSP=='jfr'){
  # K0 = K0+sum(K0_m) #
  K0 = 50 # 50 #
  # K0 = 60 # 85 # 60 # 100 #
}

nMCMC = 20000 # 10000 # 
nBurnIn = 15000 # 5000 # 
nThin = 10 # 1 #  

hp_mcmc <- list(
  t0=-0.5, t1=-5e-4, t0_adapt=200, seed=123,
  # t0=-0.5, t1=-3e-4, t0_adapt=200, seed=123,
  a_sig=3, b_sig=1, prec0=1/4,
  a_theta=0.5, b_theta=0.1, var_spike_y=0.005,
  a_m=3, b_m=1,
  a_chi=0.5, b_chi=0.1, var_spike=0.005,
  alpha_L=5, alpha_G=5,
  # alpha_L=10, alpha_G=10,
  a_xi=3, b_xi=2
)

if(whichCUSP=='jfr'){
  hp_mcmc$alpha_L = hp_mcmc$alpha_L+Data$M*hp_mcmc$alpha_G 
}

# Gibbs Sampler ----------------------------------------------------------------
print(' | Running MCMC and computing predicitions ')

print(paste0('Which CUSPs: ',whichCUSP))

init_time <- proc.time()
if(whichCUSP %in% c('d-cusp','i-cusp','i-cusp-naive')){
  ris_MCMC <- gibbs_jafar(Data$X_m, y_run, yBinary=binary_y, K0=K0, K0_m=K0_m,
                          which_prior=whichCUSP,hyperparams = hp_mcmc,
                          nMCMC=nMCMC, nBurnIn=nBurnIn, nThin=nThin,
                          get_latent_vars=T, rescale_pred=T, pow_tempering=0)
} else {
  ris_MCMC <- gibbs_jfr(Data$X_m, y_run, yBinary=binary_y, K0=K0,
                        hyperparams = hp_mcmc,
                        nMCMC=nMCMC, nBurnIn=nBurnIn, nThin=nThin,
                        get_latent_vars=T, rescale_pred=T, pow_tempering=0)
}
end_time <- proc.time() - init_time
time_run = end_time["user.self"] + end_time["sys.self"]

ris_JAFAR <- list(ris_MCMC=ris_MCMC,time_run=time_run)

# Response Predictions ---------------------------------------------------------

if(run_supervised){
  init_time = proc.time()
  y_JAFAR_train <- jafar_pred_y(Data$X_m,ris_MCMC,rescale_pred=T)
  y_JAFAR_test <- jafar_pred_y(Data$X_m_test,ris_MCMC,rescale_pred=T)
  end_time <- proc.time() - init_time
  time_pred = end_time["user.self"] + end_time["sys.self"]
  
  ris_JAFAR <- c(ris_JAFAR,list(time_pred=time_pred,
                                y_JAFAR_train=y_JAFAR_train,
                                y_JAFAR_test=y_JAFAR_test))
  
  # if(!run_simulations){
  #   y_JAFAR_valid <- jafar_pred_y(Data$X_m_val,ris_MCMC,rescale_pred=T)
  #   ris_JAFAR$y_JAFAR_valid
  # }
}

# Omics Predictions ------------------------------------------------------------

if(predict_views){
  Xm_JAFAR_train <- jafar_pred_X(Data$X_m,ris_MCMC,rescale_pred=T)
  # Xm_JAFAR_test <- jafar_pred_X(Data$X_m_test,ris_MCMC,rescale_pred=T)
  
  nMCtest = round((nMCMC-nBurnIn)/nThin)
  
  Cov_m_mean <- lapply(1:Data$M, function(m) matrix(0,Data$p_m[m],Data$p_m[m]))
  Marg_Var_m <- lapply(1:Data$M, function(m) matrix(NA,nMCtest,Data$p_m[m]))
  
  for(t in 1:nMCtest){
    for(m in 1:Data$M){
      Cov_m_mean[[m]] <- Cov_m_mean[[m]] + cov(Xm_JAFAR_train[[m]][,,t]) / (nMCtest)
      Marg_Var_m[[m]][t,] <- apply(Xm_JAFAR_train[[m]][,,t],2,var)
    }
  }
  
  ris_JAFAR$Cov_m_mean <- Cov_m_mean
  ris_JAFAR$Marg_Var_m <- Marg_Var_m
}

# Saving Output ----------------------------------------------------------------
print(' | Saving Output ')

if(T){
  repCount <- 0
  if(run_supervised){data_file <- paste0(data_file,'_y')}
  fileName <- paste0(data_file,'_',whichCUSP,'_nMC',nMCMC,'_nBurn',nBurnIn,'_nThin',nThin)
  
  risFile  <- paste0(fileName,'_rep',repCount,'.rds')
  
  while(file.exists(file.path(out_dir, risFile))){
    repCount <- repCount + 1
    risFile <- paste(fileName,'_rep',repCount,'.rds',sep='')
  }
  
  saveRDS(ris_JAFAR,file.path(out_dir, risFile))
  
  print(risFile)
}





