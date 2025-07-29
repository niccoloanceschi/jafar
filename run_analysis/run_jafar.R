
rm(list = ls())

run_cluster = F

if(run_cluster){
  path <- '/hpc/group/dunsonlab/na224/multiomics/' # cluster simulations
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

whichCUSP = 'i-cusp' # 'jfr' # 
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
  
  out_dir <- paste0(out_dir,'sim_sec3')
  if(!dir.exists(out_dir)){dir.create(out_dir)}
  
  if(run_cluster){
    if (length(commandArgs(trailingOnly = TRUE)) != 2) {
      stop("Usage: your_script.R <n> <s>")
    }
    nn <- as.integer(commandArgs(trailingOnly = TRUE)[2])
    ss <- as.integer(commandArgs(trailingOnly = TRUE)[1])
  } else {
    nn <- 50 # 100 # 200 # 
    ss <- 10 # 5 # 
  }
  data_file <- paste0('Simulated_data_n200_s',ss)
} else {
  if(run_cluster){
    if (length(commandArgs(trailingOnly = TRUE)) < 1) {
      stop("Usage: your_script.R <s>")
    }
    ss <- as.integer(commandArgs(trailingOnly = TRUE)[1])
  } else {
    ss <- 1
  }
  data_file <- 'StelzerEGA_cv-'
}

if(run_simulations){
  data_dir = "data/Sec3_Simulations"
  data_path = file.path(data_dir, paste0(data_file,'.rds'))
} else {
  data_dir = "data/Sec4_StelzerEGA"
  data_path = file.path(data_dir, paste0(data_file,ss,'_copula.rds'))
}
if(!dir.exists(data_dir)){stop("Input data folder not found", call.=FALSE)}
if(!file.exists(data_path)){stop("Data file not found")}

data_path <- "~/Desktop/StelzerEGA_Data_cleaning/StelzerEGA_data_raw.rds" # TODO: remove
Data <- readRDS(data_path)

if(run_simulations){
  Data <- data_scale_subset(Data,nn)
  data_file <- paste0('Sec3_Simulated_data_n',nn,'_s',ss)
} else {
  data_file <- paste0(data_file,'s',ss)
}

y_run=NULL
if(run_supervised){
  y_run = Data$yTrain
}
binary_y=as.logical(min(Data$yTrain%%1==0))

# MCMC Parameters --------------------------------------------------------------

K0 = 20 # 40 #
K0_m = c(20,20,20) # rep(35,Data$M) #
if(whichCUSP=='jfr'){K0=K0*2}

nMCMC = 20000 # 20000 # 
nBurnIn = 15000 # 15000 #  
nThin = 10 # 1 #  

hp_mcmc <- list(
  t0=-1, t1=-1e-4, t0_adapt=200, seed=123,
  # a_sig=3, b_sig=1, prec0=1/4,
  a_sig=5, b_sig=1, prec0=1/4,
  a_theta=0.5, b_theta=0.1, var_spike_y=0.005,
  # a_theta=2, b_theta=0.2, var_spike_y=0.005,
  # a_m=3, b_m=1,
  a_m=5, b_m=1,
  a_chi=0.5, b_chi=0.1, var_spike=0.005,
  # a_chi=2, b_chi=0.2, var_spike=0.005,
  alpha=5, alpha_loc=5,
  # alpha=10, alpha_loc=10,
  a_xi=3, b_xi=2
)
if(whichCUSP %in% c('i-cusp','i-cusp-naive')){
  hp_mcmc$alpha = hp_mcmc$alpha / sqrt(Data$M)
  # hp_mcmc$alpha = hp_mcmc$alpha 
} else if(whichCUSP=='jfr'){
  hp_mcmc$alpha = sum(hp_mcmc$alpha+Data$M*hp_mcmc$alpha_loc) / sqrt(Data$M)
}

# Gibbs Sampler ----------------------------------------------------------------
print(' | Running MCMC and computing predicitions ')

print(paste0('Which CUSPs: ',whichCUSP))

init_time <- proc.time()
if(whichCUSP %in% c('d-cusp','i-cusp','i-cusp-naive')){
  ris_MCMC <- gibbs_jafar(Data$X_m, y_run, yBinary=binary_y, K0=K0, K0_m=K0_m,
                          which_prior=whichCUSP,hyperparams = hp_mcmc,
                          nMCMC=nMCMC, nBurnIn=nBurnIn, nThin=nThin,
                          get_latent_vars=T, rescale_pred=T,pow_rescale=1/2)
} else {
  ris_MCMC <- gibbs_jfr(Data$X_m, y_run, yBinary=binary_y, K0=K0,
                        hyperparams = hp_mcmc,
                        nMCMC=nMCMC, nBurnIn=nBurnIn, nThin=nThin,
                        get_latent_vars=T, rescale_pred=T,pow_rescale=2/3)
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
}

# Omics Predictions ------------------------------------------------------------

if(predict_views){
  Xm_JAFAR_train <- jafar_pred_X(Data$X_m,ris_MCMC,rescale_pred=T)
  Xm_JAFAR_test <- jafar_pred_X(Data$X_m_test,ris_MCMC,rescale_pred=T)
  
  ris_JAFAR$Xm_JAFAR_train <- Xm_JAFAR_train
  ris_JAFAR$Xm_JAFAR_test <- Xm_JAFAR_test
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





