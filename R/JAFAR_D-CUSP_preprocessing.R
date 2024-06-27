
#' Check in any given variable is a scalar 
#' 
#' @param x input variable
#' @return Boolean value
#' 
is.scalar <- function(x) is.atomic(x) && length(x) == 1L

#' Identification of missing entries in multi-view predictors
#' 
#' @param X_m Multi-view predictors
#' @return View-wise list if indices of meassing features
#'
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

#' Empirical cdf transform
#' 
#' @param x Predictor values
#'
order_index = function(x){
  ind = sort(x, index.return = T)$ix
  x[ind] = (1:length(x))/(length(x)+1)
  return(x)
}

#' Empirical cdf transform with missing data
#' 
#' @param x Predictor values
#'
order_index_na = function(x){
  x_na = is.na(x)
  x[x_na == F] = order_index(x[x_na == F])
  n = length(x) - sum(x_na)
  return(x*n/(n+1))
}

#' Point-wise evaluation empirical cdf function
#' 
#' @param xNew location at which the cdf function is to be evaluated
#' @param xTrain data to construct the empirical cdf
#'
F_hat <- function(xNew,xTrain){
  x_obs <- xTrain[!is.na(xTrain)]
  n = length(x_obs)
  knots = c(-Inf,sort(x_obs))
  values = c(1/(n+1),c(1:n)/(n+1))
  
  sapply(xNew,function(x) max(values[knots<=x]))
}

#' Point-wise evaluation empirical quantile function
#' 
#' @param pNew probability at which the quantile function is to be evaluated
#' @param xTrain data to construct the empirical cdf
#'
Q_hat <- function(pNew,xTrain) {
  x_obs <- xTrain[!is.na(xTrain)]
  n <- length(x_obs)
  knots = c(sort(x_obs),max(x_obs))
  values = c(1:n)/(n+1)
  
  knots[sapply(pNew, function(p) sum(p > values) + 1)]
}

#' Smoothed version of the empirical cdf
#' 
#' @param vec Predictor values
#'
get_F_smooth <- function(vec){
  vec <- vec[!is.na(vec)]
  x_range  <- c(-100,sort(vec),100)
  # y_range  <- c(0,sort(order_index_na(vec)),1)
  y_range  <- c(0.001,sort(order_index_na(vec)),0.999)
  F_spline <- splinefun(x=x_range, y=y_range, method = "hyman", ties= list("ordered", mean))
  return(F_spline)
}

#' Smoothed version of the quantile cdf
#' 
#' @param vec Predictor values
#'
get_Q_smooth <- function(vec){
  vec <- vec[!is.na(vec)]
  x_range  <- c(0,sort(order_index(vec)),1)
  y_range  <- c(min(vec),sort(vec),max(vec))
  Q_spline <- splinefun(x=x_range, y=y_range, method = "hyman", ties= list("ordered", mean))
  return(Q_spline)
}

#' apply cdf transform to multi-view predictors
#' 
#' @param Z_m Train set predictors
#' @param Z_m_test Test set predictors
#' @param smoothed Use smoothed cdf
#' @return List of transformed features
#'
cdf_transform <- function(Z_m,Z_m_test=NULL,smoothed=F){
  p_m <- sapply(Z_m,ncol)
  for(m in 1:length(Z_m)){
    # Train Set
    Z_m[[m]] = qnorm(apply(Z_m[[m]], 2, order_index_na))
    # Test Set
    if(!is.null(Z_m_test)){
      for(j in 1:p_m[m]){
        if(smoothed){
          F_smooth <- get_F_smooth(Z_m[[m]][,j])
          Z_m_test[[m]][,j] <- qnorm(F_smooth(Z_m_test[[m]][,j]))
        } else{
          Z_m_test[[m]][,j] <- qnorm(F_hat(Z_m_test[[m]][,j], Z_m[[m]][,j]))
        }
      }
    }
  }
  output <- list(Z_m=Z_m)
  if(!is.null(Z_m_test)){output$Z_m_test=Z_m_test}
  return(output)
}

#' JAFAR multi-view predictors preprocess: cdf transform + center & rescale
#' 
#' @param Z_m Train set predictors
#' @param Z_m_test Test set predictors
#' @param copula Apply cdf transformation
#' @param smoothed Use smoothed cdf
#' @return List of preprocessed features and rescaling factors
#'
#' @export
#' 
jafar_preprocess_X <- function(Z_m,Z_m_test=NULL,copula=F,smoothed=F){
  
  preprocess_X_m <- list()

  n <- sapply(Z_m,nrow)[1]
  if(!is.null(Z_m_test)){nTest <- sapply(Z_m_test,nrow)[1]}
  
  M <- length(Z_m)
  p_m <- sapply(Z_m,ncol)
  
  names_col <- lapply(Z_m,colnames)
  if(any(sapply(names_col,is.null))){
    names_col <- lapply(p_m, seq)
    for(m in 1:M){
      colnames(Z_m[[m]]) <- names_col
      if(!is.null(Z_m_test)){colnames(Z_m_test[[m]]) <- names_col}
    }
  }
  
  names_train <- rownames(Z_m[[1]])
  if(is.null(names_train)){
    names_train <- seq(n)
    for(m in 1:M){rownames(Z_m[[m]]) <- names_train}
  }
  
  if(!is.null(Z_m_test)){
    names_test <- rownames(Z_m_test[[1]])
    if(is.null(names_test)){
      names_test <- seq(nTest)
      for(m in 1:M){rownames(Z_m_test[[m]]) <- names_test}
    }
  }
  
  if(copula){
    output <- cdf_transform(Z_m,Z_m_test=NULL,smoothed=F)
    Z_m = output$Z_m
    if(!is.null(Z_m_test)){Z_m_test = output$Z_m_test}
  }
  
  for(m in 1:M){
    # Center and Scale Predictors
    preprocess_X_m[[m]] = preProcess(Z_m[[m]], method = c("center", "scale"))
    
    Z_m[[m]] = as.matrix(predict(preprocess_X_m[[m]], Z_m[[m]]))
    rownames(Z_m[[m]]) <- names_train
    
    if(!is.null(Z_m_test)){
      Z_m_test[[m]] = as.matrix(predict(preprocess_X_m[[m]], Z_m_test[[m]]))
      rownames(Z_m_test[[m]]) <- names_test
    }
  }
  
  output <- list(Z_m=Z_m,preprocess_X_m=preprocess_X_m)
  if(!is.null(Z_m_test)){output$Z_m_test=Z_m_test}
  return(output)
}

#' JAFAR response preprocess: center & rescale
#' 
#' @param yTrain Train set responses
#' @param yTest Test set responses
#' @return List of preprocessed responses and rescaling factors
#' 
#' @export
#' 
jafar_preprocess_y <- function(yTrain,yTest=NULL){
  
  yTrain <- as.matrix(yTrain,ncol=1)
  colnames(yTrain) <- c("response")
  
  names_train <- rownames(yTrain)
  if(is.null(names_train)){names_train <- seq(nrow(yTrain))}
  
  # Center and Scale Responses
  preprocess_y = preProcess(yTrain, method = c("center", "scale"))
  
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












