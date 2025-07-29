
rm(list = ls())

run_cluster = F

if(run_cluster){
  path <- '/hpc/group/dunsonlab/na224/multiomics/' # cluster simulations
} else{ 
  path <- '~/Documents/GitHub/jafar/' # local simulations
}
setwd(path)

# Loading Packages -------------------------------------------------------------

if(FALSE){
  # devtools::install_github("sarahsamorodnitsky/BSFP")
  if('BSFP' %in% (.packages())){ detach("package:BSFP", unload = TRUE) }
  remove.packages("BSFP")
  # install.packages("~/Documents/PostDoc/Projects/Multi-Omics/Code/BSFP-master/", repos = NULL, type = "source")
  install.packages("~/Documents/PostDoc/Projects/Multi-Omics/Code/BSFP-Nico/", repos = NULL, type = "source")
}

library(BSFP)

# Files and Folder -------------------------------------------------------------

## Select model and data ----

run_simulations = F
run_supervised = T
predict_views = T

## Source code ----

source_dir = "new_source"
if(!dir.exists(source_dir)){stop("Source folder not found", call.=FALSE)}

source_files = c("bsfp_predict_oos.R", "Sec3_simulations_preproc.R")

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
} else {
  data_file <- paste0(data_file,'s',ss)
}

## Formatting Input Data -------------------------------------------------------

input_format_description <- "
Note that in the output of the bsfp_data function, the data are stored in a matrix of lists.
To access source 1, for example, we would run data.sim$data[[1,1]] and source 2 data.sim$data[[2,1]].
The response vector, data.sim$Y, is stored in a similar way and would be accessed using data.sim$Y[[1,1]]"

X.train <- X.test <- matrix(list(), nrow = Data$M, ncol = 1)
Y.train <- Y.test <- matrix(list(), nrow = 1, ncol = 1)

for(m in 1:Data$M) {
  X.train[[m,1]] <- t(Data$X_m[[m]]) 
  X.test[[m,1]]  <- t(Data$X_m_test[[m]])
}

Y.train[[1,1]] <- Data$yTrain
Y.test[[1,1]]  <- Data$yTest 

y_run=NULL
if(run_supervised){
  y_run = Y.train
}
binary_y=as.logical(min(Data$yTrain%%1==0))

# MCMC Parameters --------------------------------------------------------------

nMCMC = 20000 # 10000 # 
nBurnIn = 15000 # 6000 # 
nThin = 10 # 1 #  

nMCtest = floor(nMCMC/nThin)

# Gibbs Sampler ----------------------------------------------------------------
print(' | Running MCMC and computing predicitions ')

init_time <- proc.time()
bsfp_mcmc <- bsfp(data = X.train, Y = y_run, nsample = nMCMC, thinning = nThin, 
                  # ranks = rep(2,Data$M+1),
                  progress = T, save_init = T, save_structures = F,
                  save_predictive_model = T, save_loadings_scores = T,
                  save_imputations = F, save_last_sample = T)
end_time <- proc.time() - init_time
time_run = end_time["user.self"] + end_time["sys.self"]

ris_BSFP <- list(bsfp_mcmc=bsfp_mcmc,time_run=time_run)

# Response Predictions ---------------------------------------------------------

if(run_supervised){

  print(" | Response OOS Prediction")
  
  init_time <- proc.time()
  bsfp_train <- bsfp.predict.oos(bsfp.fit=bsfp_mcmc, test_data=X.train, nsample=nMCtest)
  bsfp_test <- bsfp.predict.oos(bsfp.fit=bsfp_mcmc, test_data=X.test, nsample=nMCtest)
  end_time <- proc.time() - init_time
  time_pred = end_time["user.self"] + end_time["sys.self"]
  
  ris_BSFP <- c(ris_BSFP,list(time_pred=time_pred,
                              train_pred=bsfp_train,
                              test_pred=bsfp_test))
}

# Omics Predictions ------------------------------------------------------------

if(predict_views){
  
  print(" | Views OOS Prediction")
  
  nBurnInTest = 0 # floor(nMCtest/2) # 
  
  Xm_train_bsfp <- lapply(1:Data$M, function(m) array(NA,c(nMCtest-nBurnInTest,Data$p_m[m],Data$n)))
  Xm_test_bsfp  <- lapply(1:Data$M, function(m) array(NA,c(nMCtest-nBurnInTest,Data$p_m[m],Data$nTest)))
  # double check
  tXm_train_bsfp <- Xm_train_bsfp
  tXm_test_bsfp  <- Xm_test_bsfp
    
  for(m in 1:Data$M){
    
    X.train.NA <- X.train
    X.train.NA[[m,1]] <- NA*X.train[[m,1]]
    
    X.test.NA <- X.test
    X.test.NA[[m,1]] <- NA*X.test[[m,1]]
    
    bsfp.X.train <- bsfp.predict.oos(bsfp.fit=bsfp_mcmc, test_data=X.train.NA, nsample=nMCtest)
    bsfp.X.test <- bsfp.predict.oos(bsfp.fit=bsfp_mcmc, test_data=X.test.NA, nsample=nMCtest)
    
    for(t in c((nBurnInTest+1):(nMCtest-1))){
      Xm_train_bsfp[[m]][t-nBurnInTest,,] <- matrix(bsfp.X.train$Xm.draw[[t]][[m,1]],ncol=Data$n)
      Xm_test_bsfp[[m]][t-nBurnInTest,,]  <- matrix(bsfp.X.test$Xm.draw[[t]][[m,1]],ncol=Data$nTest)
      # double check 
      tXm_train_bsfp[[m]][t-nBurnInTest,,] <- t(matrix(bsfp.X.train$Xm.draw[[t]][[m,1]],nrow=Data$n))
      tXm_test_bsfp[[m]][t-nBurnInTest,,]  <- t(matrix(bsfp.X.test$Xm.draw[[t]][[m,1]],nrow=Data$nTest))  
    }
  }
  
  ris_BSFP$Xm_train_bsfp <- Xm_train_bsfp
  ris_BSFP$Xm_test_bsfp <- Xm_test_bsfp
  # double check 
  ris_BSFP$tXm_train_bsfp <- tXm_train_bsfp
  ris_BSFP$tXm_test_bsfp <- tXm_test_bsfp
}

# Saving Output ----------------------------------------------------------------
print(' | Saving Output ')

repCount <- 0
if(run_supervised){data_file <- paste0(data_file,'_y')}
fileName <- paste0(data_file,'_bsfp_nMC',nMCMC,'_nBurn',nBurnIn,'_nThin',nThin)

risFile  <- paste0(fileName,'_rep',repCount,'.rds')

while(file.exists(file.path(out_dir, risFile))){
  repCount <- repCount + 1
  risFile <- paste(fileName,'_rep',repCount,'.rds',sep='')
}

saveRDS(ris_BSFP,file.path(out_dir, risFile))

print(risFile)





