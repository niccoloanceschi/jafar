
rm(list = ls())

run_cluster = F

if(run_cluster){
  path <- '/hpc/group/dunsonlab/na224/multiomics/' # cluster simulations
} else{ 
  path <- '~/Documents/GitHub/jafar/' # local simulations
}
setwd(path)

# Loading Packages -------------------------------------------------------------

library(multiview)

if(FALSE){ # Install `bartMachine`
  Sys.setenv(JAVA_HOME = '/Library/Internet Plug-Ins/JavaAppletPlugin.plugin/Contents/Home/')
  install.packages("rJava")
  install.packages("bartMachine")
}

if(FALSE){
  detach("package:IntegratedLearner", unload = TRUE)
  remove.packages("IntegratedLearner")
  install.packages("/Users/nico/Downloads/mcmcplots_0.4.tar.gz",repos = NULL,type = "source")
  install.packages("~/Documents/PostDoc/Projects/Multi-Omics/Code/IntegratedLearner-master/", repos = NULL, type = "source")
  devtools::install_github("himelmallick/IntegratedLearner")
}

Sys.setenv(JAVA_HOME = '/Library/Internet Plug-Ins/JavaAppletPlugin.plugin/Contents/Home/')
library(IntegratedLearner)
options(java.parameters = "-Xmx5g") # This is needed for running BART
library(SuperLearner)
library(bartMachine)
library(magrittr)

# Files and Folder -------------------------------------------------------------

## Select model and data ----

run_simulations = F

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

if(run_simulations){
  out_dir <- paste0(out_dir,'sim_sec3')
  if(!dir.exists(out_dir)){dir.create(out_dir)}
}

## Data import  ----

print(' | Loading Data ')

s_values = c(1:5)
# s_values = c(1:100)
n_values = c(50,100,150,200)

for(nn in n_values){
  for(ss in s_values){

    if(run_simulations){
      print(paste0(' --- Running Datset with n=',nn))
      data_file <- paste0('Simulated_data_n200_s',ss)
    } else {
      data_file <- 'StelzerEGA_cv-'
    }
    
    print(paste0(' --- Running Replica with s=',ss))
    
    if(run_simulations){
      data_dir = "data/Sec3_Simulations/"
      data_path = file.path(data_dir, paste0(data_file,'.rds'))
    } else {
      data_dir = "data/Sec4_StelzerEGA"
      data_path = file.path(data_dir, paste0(data_file,ss,'.rds'))
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
    
    binary_y=as.logical(min(Data$yTrain%%1==0))
    if(binary_y){
      family = stats::binomial(link = "probit")
    } else {
      family = gaussian()
    }
    
    # Fitting CoopL ----------------------------------------------------------------
    
    ## run model -----
    print(' | Running CoopL ')
    
    init_time <- proc.time() 
    coop_fit <- cv.multiview(x_list=Data$X_m, y=Data$yTrain, family=family,
                             rho=0.5, standardize=F, maxit=1e3, 
                             nfolds=length(Data$yTrain))
    end_time <- proc.time() - init_time
    time_run = end_time["user.self"] + end_time["sys.self"]
    
    ### Compute Response Predictions ----
    
    init_time <- proc.time() 
    yhat_coop_train = as.vector(predict(object = coop_fit,newx = Data$X_m, s="lambda.min"))
    yhat_coop_test = as.vector(predict(object = coop_fit,newx = Data$X_m_test, s="lambda.min"))
    
    if(binary_y){
      yhat_coop_train = pnorm(yhat_coop_train)
      yhat_coop_test = pnorm(yhat_coop_test)
    }
    end_time <- proc.time() - init_time
    time_pred = end_time["user.self"] + end_time["sys.self"]
    
    names(yhat_coop_train) <- rownames(Data$yTrain)
    names(yhat_coop_test) <- rownames(Data$yTest)
    
    risCoopL <- list(fit_coopL=coop_fit,time_run=time_run,time_pred=time_pred,
                     y_pred_train=yhat_coop_train,y_pred_test=yhat_coop_test)
    
    ## saving output ------
    
    repCount <- 0
    fileName <- paste0(data_file,'_y_','CoopL')
    
    risFile  <- paste0(fileName,'_rep',repCount,'.rds')
    
    while(file.exists(file.path(out_dir, risFile))){
      repCount <- repCount + 1
      risFile <- paste(fileName,'_rep',repCount,'.rds',sep='')
    }
    
    saveRDS(risCoopL,file.path(out_dir, risFile))
    
    print(risFile)
    
    # Fitting IntegL ---------------------------------------------------------------
    
    ## Format Data ----
    
    Data$yTrain <- matrix(Data$yTrain,ncol = 1)
    Data$yTest <- matrix(Data$yTest,ncol = 1)
    
    if(!('omics_types' %in% names(Data))){ Data$omics_types <- paste0(c(1:Data$M)) }
    if(!('omics_names' %in% names(Data))){
      Data$omics_names <- lapply(1:Data$M, function(m) paste0(m,'_',c(1:Data$p_m[m])))
    }
    if(is.null(rownames(Data$yTrain))){ rownames(Data$yTrain) <- paste0('n-',c(1:Data$n))}
    if(is.null(rownames(Data$yTest))){ rownames(Data$yTest) <- paste0('nTest-',c(1:Data$nTest))}
    
    # TRAIN SET
    
    # (i) feature_table
    feature_table <- t(matrix(unlist(Data$X_m),Data$n,sum(Data$p_m)))
    colnames(feature_table) <- rownames(Data$yTrain)
    rownames(feature_table) <- unlist(Data$omics_names)
    
    # (ii) sample_metadata
    sample_metadata <- data.frame(Y=unname(Data$yTrain),
                                  subjectID=sapply(strsplit(rownames(Data$yTrain), "_"), function(x) x[[1]]))
    rownames(sample_metadata) <- rownames(Data$yTrain)
    
    # (iii) feature_metadata
    feature_metadata <- data.frame(featureID=unlist(Data$omics_names),
                                   featureType=unlist(lapply(1:Data$M, function(m) rep(Data$omics_types[m],Data$p_m[m]))))
    rownames(feature_metadata) <- unlist(Data$omics_names)
    
    # TEST SET
      
    # (i) feature_table
    feature_table_test <- t(matrix(unlist(Data$X_m_test),Data$nTest,sum(Data$p_m)))
    colnames(feature_table_test) <- rownames(Data$yTest)
    rownames(feature_table_test) <- unlist(Data$omics_names)
    
    # (ii) sample_metadata
    sample_metadata_test <- data.frame(Y=unname(Data$yTest),
                                       subjectID=sapply(strsplit(rownames(Data$yTest), "_"), function(x) x[[1]]))
    rownames(sample_metadata_test) <- rownames(Data$yTest)
    
    ## run model ----
    print(' | Running IntegL ')
    
    init_time <- proc.time()
    fit_IntegL<-IntegratedLearner(feature_table = feature_table,
                           sample_metadata = sample_metadata, 
                           feature_metadata = feature_metadata,
                           base_learner = 'SL.BART',
                           meta_learner = 'SL.nnls.auc',  
                           family = family,
                           folds = 5,verbose = F)
    end_time <- proc.time() - init_time
    time_run = end_time["user.self"] + end_time["sys.self"]
    
    ### Compute Response Predictions ----
    
    init_time <- proc.time()
    
    weights <- fit_IntegL$weights
    
    dataX <- fit_IntegL$X_train_layers
    
    new_dataX <- lapply(Data$X_m_test,data.frame)
    for(m in 1:Data$M){colnames(new_dataX[[m]]) <- Data$omics_names[[m]]}
    names(new_dataX) <- Data$omics_types
    
    post.samples.train <- vector("list", length(weights))
    names(post.samples.train) <- names(dataX)
    
    post.samples.test <- vector("list", length(weights))
    names(post.samples.test) <- names(new_dataX)
    
    for(i in seq_along(post.samples.train)){
      post.samples.train[[i]] <- bart_machine_get_posterior(fit_IntegL$model_fits$model_layers[[i]],dataX[[i]])$y_hat_posterior_samples
    }
    
    for(i in seq_along(post.samples.test)){
      post.samples.test[[i]] <- bart_machine_get_posterior(fit_IntegL$model_fits$model_layers[[i]],new_dataX[[i]])$y_hat_posterior_samples
    }
    
    weighted.post.samples.train <-Reduce('+', Map('*', post.samples.train, weights))
    rownames(weighted.post.samples.train) <- rownames(dataX[[1]])
    
    weighted.post.samples.test <-Reduce('+', Map('*', post.samples.test, weights))
    rownames(weighted.post.samples.test) <- rownames(new_dataX[[1]])
    
    end_time <- proc.time() - init_time
    time_pred = end_time["user.self"] + end_time["sys.self"]
    
    risIntegLearn <- list(model_fit=fit_IntegL,time_run=time_run,
                          train_pred=weighted.post.samples.train,
                          test_pred=weighted.post.samples.test,
                          time_pred=time_pred)
    
    ## saving output ------
    
    repCount <- 0
    fileName <- paste0(data_file,'_y_','IntegL')
    
    risFile  <- paste0(fileName,'_rep',repCount,'.rds')
    
    while(file.exists(file.path(out_dir, risFile))){
      repCount <- repCount + 1
      risFile <- paste(fileName,'_rep',repCount,'.rds',sep='')
    }
    
    saveRDS(risIntegLearn,file.path(out_dir, risFile))
    
    print(risFile)
  }
  if(!run_simulations){stop()}
}


