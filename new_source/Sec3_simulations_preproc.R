
data_scale_subset <- function(Data,nn=NULL){
  
  if(is.null(nn)){nn=Data$n}
  
  Data$eta    <- Data$eta[1:nn,]
  Data$X_m    <- lapply(Data$X_m,function(df) df[1:nn,])
  Data$phi_m  <- lapply(Data$phi_m,function(df) df[1:nn,])
  Data$n      <- nn
  
  Data$preprocess_X_m <- list()
  for(m in 1:Data$M){
    # Center and Scale Omics Data  
    colnames(Data$X_m[[m]]) = seq(ncol(Data$X_m[[m]]))
    Data$preprocess_X_m[[m]] = caret::preProcess(Data$X_m[[m]], method = c("center", "scale"))
    Data$X_m[[m]] = as.matrix(predict(Data$preprocess_X_m[[m]], Data$X_m[[m]]))
    if('X_m_test' %in% names(Data)){
      colnames(Data$X_m_test[[m]]) = seq(ncol(Data$X_m_test[[m]]))
      Data$X_m_test[[m]] = as.matrix(predict(Data$preprocess_X_m[[m]], Data$X_m_test[[m]]))
    }
  }
  
  if('yTrain' %in% names(Data)){
    Data$yTrain <- Data$yTrain[1:nn]
    Data$yTrain <- matrix(Data$yTrain,ncol=1)
    # Center and Scale
    colnames(Data$yTrain) <- c("response")
    Data$preprocess_y = caret::preProcess(Data$yTrain, method = c("center", "scale"))
    Data$yTrain = as.vector(predict(Data$preprocess_y, Data$yTrain))
    if('yTest' %in% names(Data)){
      Data$yTest <- matrix(Data$yTest,ncol=1)
      colnames(Data$yTest)  <- c("response")
      Data$yTest = as.vector(predict(Data$preprocess_y, Data$yTest))
    }
  }
  
  return(Data)
}