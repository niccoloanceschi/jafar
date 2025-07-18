
rm(list = ls())

run_cluster = F

if(run_cluster){
  path <- '/hpc/group/dunsonlab/na224/multiomics/' # cluster simulations
} else{ 
  path <- '~/Documents/GitHub/jafar/' # local simulations
}
setwd(path)

# Loading Packages -------------------------------------------------------------

if(F){
  library(devtools)
  install_github('chekouo/BIPnet')
}
library(BIPnet)

# Files and Folder -------------------------------------------------------------

## Select model and data ----

run_simulations = T
run_supervised = T

## Source code ----

source_dir = "new_source"
if(!dir.exists(source_dir)){stop("Source folder not found", call.=FALSE)}

source_files = c("Sec3_simulations_preproc.R")

for(s_file in source_files){
  if(!file.exists(file.path(source_dir, s_file))){stop(paste0("Source file '",s_file,"' not found"))}
  source(file.path(source_dir, s_file))
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
    nn <- as.integer(commandArgs(trailingOnly = TRUE)[1])
    ss <- as.integer(commandArgs(trailingOnly = TRUE)[2])
  } else {
    nn <- 50 # 100 # 200 # 
    ss <- 10 # 5 # 
  }
  data_file <- paste0('Simulated_data_n200_s',ss)
} else {
  if(run_cluster){
    if (length(commandArgs(trailingOnly = TRUE)) != 1) {
      stop("Usage: your_script.R <s>")
    }
    ss <- as.integer(commandArgs(trailingOnly = TRUE)[1])
} else {
  ss <- 1
  }
  data_file <- 'StelzerEGA_cv-'
}

if(run_simulations){
  data_dir = "data/Sec3_Simulations/"
  data_path = file.path(data_dir, paste0(data_file,'.rds'))
} else {
  data_dir = "data/Sec4_StelzerEGA/"
  data_path = file.path(data_dir, paste0(data_file,ss,'_copula.rds'))
}
if(!dir.exists(data_dir)){stop("Input data folder not found", call.=FALSE)}
if(!file.exists(data_path)){stop("Data file not found")}
Data <- readRDS(data_path)

if(run_simulations){
  Data <- data_scale_subset(Data,nn)
  data_file <- paste0('Sec3_Simulated_data_n',nn,'_s',ss)
}

## Formatting Input Data -------------------------------------------------------

Xy = Data$X_m
var_types=rep(0,Data$M)

if(run_supervised){
  Xy[Data$M+1] = Data$yTrain
  var_types=c(var_types,1)
}

# MCMC Parameters --------------------------------------------------------------

nMCMC = 5000 # 
nBurnIn = nMCMC/2

Ktot = 15 # 20 # 

# Gibbs Sampler ----------------------------------------------------------------
print(' | Running MCMC and computing predicitions ')

timerun = NULL
tictoc::tic()
BIPrun=BIP(dataList=Xy,IndicVar=var_types,Method="BIP",
           nbrcomp=Ktot,sample=nMCMC,burnin=nBurnIn)
timerun = tictoc::toc()

# Response Predictions ---------------------------------------------------------

ris_BIP <- list(runBIP=BIPrun,time_run=time_run)

if(run_supervised){
  
  trainPred = BIPpredict(Data$X_m,Result=BIPrun,meth='BMA')
  testPred = BIPpredict(Data$X_m_test,Result=BIPrun,meth='BMA')
  
  ris_BIP <- c(ris_BIP,train_pred=trainPred,test_pred=testPred)
}

# Reconstructed Covariance -----------------------------------------------------

cov_mean_m <- list()
for(m in 1:Data$M){
  cov_estim_m[[m]] = diag(BIPrun$EstSig2[[m]],Data$p_m[m],Data$p_m[m])
  for(t in 1:length(BIPrun$EstLoadModel)){
    cov_estim_m[[m]] = tmp + BIPrun$PostGam[[t]]*crossprod(BIPrun$EstLoadModel[[t]][[m]])
  }
}

ris_BIP$cov_mean_m = cov_mean_m

# Saving Output ----------------------------------------------------------------
print(' | Saving Output ')

repCount <- 0
if(run_supervised){data_file <- paste0(data_file,'_y')}
fileName <- paste0(data_file,'_bip_nMC',nMCMC,'_nBurn',nBurnIn,'_nThin',nThin)

risFile  <- paste0(fileName,'_rep',repCount,'.rds')

while(file.exists(file.path(out_dir, risFile))){
  repCount <- repCount + 1
  risFile <- paste(fileName,'_rep',repCount,'.rds',sep='')
}

saveRDS(ris_BIP,file.path(out_dir, risFile))

print(risFile)





