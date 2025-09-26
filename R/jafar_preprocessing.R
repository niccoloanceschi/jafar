
is.scalar <- function(x) is.atomic(x) && length(x) == 1L

get_NA_X <- function(X_m){
  
  M <- length(X_m)
  n <- unlist(sapply(X_m,nrow))[1]
  
  na_idx <- na_row_idx <- list()
  
  for(m in 1:M){
    
    if(!is.null(X_m[[m]])){
      
      # identifying indexes of NA
      Xm_na <- is.na(unname(as.matrix(X_m[[m]])))
      
      # Containers for output of imputed NA values
      na_idx[[m]]     <- apply(Xm_na,1,which)
      na_row_idx[[m]] <- c(1:n)[sapply(na_idx[[m]],length)>0]
      na_idx[[m]]     <- na_idx[[m]][na_row_idx[[m]]]
    }
  }
  
  return(list(na_row_idx=na_row_idx,na_idx=na_idx))  
}

order_index = function(x){
  ind = sort(x, index.return = T)$ix
  x[ind] = (1:length(x))/(length(x)+1)
  return(x)
}

order_index_na = function(x){
  x_na = is.na(x)
  x[x_na == F] = order_index(x[x_na == F])
  n = length(x) - sum(x_na)
  return(x*n/(n+1))
}

F_hat <- function(xNew,xTrain){
  x_obs <- xTrain[!is.na(xTrain)]
  n = length(x_obs)
  knots = c(-Inf,sort(x_obs))
  values = c(1/(n+1),c(1:n)/(n+1))
  
  sapply(xNew,function(x) max(values[knots<=x]))
}

Q_hat <- function(pNew,xTrain) {
  x_obs <- xTrain[!is.na(xTrain)]
  n <- length(x_obs)
  knots = c(sort(x_obs),max(x_obs))
  values = c(1:n)/(n+1)
  
  knots[sapply(pNew, function(p) sum(p > values) + 1)]
}

cdf_transform <- function(X_m,X_m_test=NULL){
  p_m <- sapply(X_m,ncol)
  for(m in 1:length(X_m)){
    # Train Set
    X_m[[m]] = qnorm(apply(X_m[[m]], 2, order_index_na))
    # Test Set
    if(!is.null(X_m_test)){
      for(j in 1:p_m[m]){
        X_m_test[[m]][,j] <- qnorm(F_hat(X_m_test[[m]][,j], X_m[[m]][,j]))
      }
    }
  }
  output <- list(X_m=X_m)
  if(!is.null(X_m_test)){output$X_m_test=X_m_test}
  return(output)
}

#' predictors preprocess: center & rescale + cdf transform (optional)
#' 
#' @param X_m Train set predictors
#' @param X_m_test Test set predictors
#' @param copula Apply cdf transformation
#' 
#' @return List of preprocessed features and rescaling factors
#'
#' @export
#' 
preprocess_X <- function(X_m,X_m_test=NULL,copula=F){
  
  preprocess_X_m <- list()

  n <- sapply(X_m,nrow)[1]
  if(!is.null(X_m_test)){nTest <- sapply(X_m_test,nrow)[1]}
  
  M <- length(X_m)
  p_m <- sapply(X_m,ncol)
  
  names_col <- lapply(X_m,colnames)
  if(any(sapply(names_col,is.null))){
    names_col <- lapply(p_m, seq)
    for(m in 1:M){
      colnames(X_m[[m]]) <- names_col
      if(!is.null(X_m_test)){colnames(X_m_test[[m]]) <- names_col}
    }
  }
  
  names_train <- rownames(X_m[[1]])
  if(is.null(names_train)){
    names_train <- seq(n)
    for(m in 1:M){rownames(X_m[[m]]) <- names_train}
  }
  
  if(!is.null(X_m_test)){
    names_test <- rownames(X_m_test[[1]])
    if(is.null(names_test)){
      names_test <- seq(nTest)
      for(m in 1:M){rownames(X_m_test[[m]]) <- names_test}
    }
  }
  
  if(copula){
    output <- cdf_transform(X_m,X_m_test=X_m_test)
    X_m = output$X_m
    if(!is.null(X_m_test)){X_m_test = output$X_m_test}
  }
  
  for(m in 1:M){
    # Center and Scale Predictors
    preprocess_X_m[[m]] = caret::preProcess(X_m[[m]], method = c("center", "scale"))
    
    X_m[[m]] = as.matrix(predict(preprocess_X_m[[m]], X_m[[m]]))
    rownames(X_m[[m]]) <- names_train
    
    if(!is.null(X_m_test)){
      X_m_test[[m]] = as.matrix(predict(preprocess_X_m[[m]], X_m_test[[m]]))
      rownames(X_m_test[[m]]) <- names_test
    }
  }
  
  output <- list(X_m=X_m,preprocess_X_m=preprocess_X_m)
  if(!is.null(X_m_test)){output$X_m_test=X_m_test}
  return(output)
}

#' response preprocess: center & rescale
#' 
#' @param yTrain Train set responses
#' @param yTest Test set responses
#' @return List of preprocessed responses and rescaling factors
#' 
#' @export
#' 
preprocess_y <- function(yTrain,yTest=NULL){
  
  yTrain <- as.matrix(yTrain,ncol=1)
  colnames(yTrain) <- c("response")
  
  names_train <- rownames(yTrain)
  if(is.null(names_train)){names_train <- seq(nrow(yTrain))}
  
  # Center and Scale Responses
  preprocess_y = caret::preProcess(yTrain, method = c("center", "scale"))
  
  yTrain = as.matrix(predict(preprocess_y, yTrain))
  rownames(yTrain) <- names_train
  
  output <- list(yTrain=yTrain,preprocess_y=preprocess_y)
  
  if(!is.null(yTest)){
    yTest <- as.matrix(yTest,ncol=1)
    colnames(yTest)  <- c("response")
    
    names_test <- rownames(yTest)
    if(is.null(names_test)){names_test <- seq(nrow(yTest))}
    
    yTest = as.matrix(predict(preprocess_y, yTest))
    rownames(yTest)  <- names_test
    
    output$yTest=yTest
  }
  
  return(output)
}

#' predictors preprocess: reorder features via hierarchical clustering for better visualization
#' 
#' @param X_m Train set predictors
#' @param X_m_test Test set predictors
#' @param K0_HC Reference number of clusters for hierarchical clustering (default: 15)
#' 
#' @return List of preprocessed features and rescaling factors
#'
#' @export
#' 
features_reorder_HC <- function(X_m,X_m_test=NULL,K0_HC=15){
  
  M <- length(X_m)
  p_m <- sapply(X_m,ncol)
  
  idx_sort_HC_m <- list()
  
  for(m in 1:M){
    idxHC_m  <- hclust(as.dist(1-cor(X_m[[m]],use='pairwise.complete.obs')),method = 'average')
    hHC_m    <- idxHC_m$height[p_m[m]-K0_HC]
    cutree_m <- cutree(tree = idxHC_m, h=hHC_m)
    
    idx_sort_HC_m[[m]] <- sort(cutree_m,index.return=T)$ix
    
    X_m[[m]] <- X_m[[m]][,idx_sort_HC_m[[m]]]
    if(!is.null(X_m_test[[m]])){X_m_test[[m]] <- X_m_test[[m]][,idx_sort_HC_m[[m]]]}
  }
  
  output <- list(X_m=X_m,idx_sort_HC_m=idx_sort_HC_m)
  if(!is.null(X_m_test)){output$X_m_test=X_m_test}
  return(output)
}












