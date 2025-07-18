
data_scale_subset <- function(Data,nn=NULL){
  
  if(is.null(nn)){nn=Data$n}
  
  Data$yTrain <- Data$yTrain[1:nn]
  Data$eta    <- Data$eta[1:nn,]
  Data$X_m    <- lapply(Data$X_m,function(df) df[1:nn,])
  Data$phi_m  <- lapply(Data$phi_m,function(df) df[1:nn,])
  Data$n      <- nn
  
  # Center and Scale
  
  Data$yTrain <- matrix(Data$yTrain,ncol=1)
  Data$yTest <- matrix(Data$yTest,ncol=1)
  
  colnames(Data$yTest)  <- c("response")
  colnames(Data$yTrain) <- c("response")
  
  Data$preprocess_y = caret::preProcess(Data$yTrain, method = c("center", "scale"))
  Data$yTrain = as.vector(predict(Data$preprocess_y, Data$yTrain))
  Data$yTest = as.vector(predict(Data$preprocess_y, Data$yTest))
  
  Data$preprocess_X_m <- list()
  for(m in 1:Data$M){
    colnames(Data$X_m[[m]]) = seq(ncol(Data$X_m[[m]]))
    colnames(Data$X_m_test[[m]]) = seq(ncol(Data$X_m_test[[m]]))
    
    # Center and Scale Omics Data  
    Data$preprocess_X_m[[m]] = caret::preProcess(Data$X_m[[m]], method = c("center", "scale"))
    Data$X_m[[m]] = as.matrix(predict(Data$preprocess_X_m[[m]], Data$X_m[[m]]))
    Data$X_m_test[[m]] = as.matrix(predict(Data$preprocess_X_m[[m]], Data$X_m_test[[m]]))
  }
  
  return(Data)
}