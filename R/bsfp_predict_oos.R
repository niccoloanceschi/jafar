
data.rearrange=function(data,rmt=F,sigma=NULL){
  out=NULL
  p=nrow(data)
  q=ncol(data)
  
  m.vec=rep(NA,p)
  n.vec= ncol(data[[1,1]]) # do.call(c, lapply(data[1,], ncol))
  
  if (is.null(sigma)) sigma=matrix(1,p,q)
  
  for (i in 1:p){
    dimm=do.call(cbind, lapply(data[i,],dim))
    m1=unique(dimm[1,])
    if (length(m1)==1 ){m.vec[i]=m1 }
    else{ stop("the number of rows do not match.") }
    if (!all(dimm[2,], n.vec)){ stop("the number of columns do not match")}
    
    for (j in 1:q){
      if (rmt) sigma[i,j]=sigma.rmt(data[[i,j]])
      data[[i,j]]=data[[i,j]]/sigma[i,j]
    }
    
    out=rbind(out,do.call(cbind,data[i,]))
  }
  
  return(list(out=out, nrows=m.vec, ncols=n.vec, sigma.mat=sigma))
}

#' Out-of-sample prediction for BSFP
#'
#' Modified version of the function \code{bsfp.predict} from the GitHub repo BSFP for out-of-sample predictions.
#'
#' @param bsfp.fit Results from fitting \code{bsfp} on training data.
#' @param test_data Matrix-list dataset of held-out test data.
#' @param response_type Continuous or binary response. Must be one of 'continuous' (deafult) or 'binary'.
#' @param model_params May be \code{NULL} if \code{model_params=NULL} in \code{bsfp} fit.
#' Otherwise, specify as \code{(error_vars, joint_vars, indiv_vars, beta_vars, response_vars)}.
#' @param nsample Integer specifing number of Gibbs sampling iterations
#' @param progress Boolean determining if progress of the sampler be displayed
#' @param starting_values List of starting values for \eqn{\mathbf{V}, \mathbf{U}_s, \mathbf{W}_s, \mathbf{V}_s} for \eqn{s=1,\dots, q}.
#' If \code{NULL}, initialize from prior.
#'
#' @details Generate new scores for held-out test data based on a
#' training fit of BSFP. Uses the estimated ranks and joint and individual loadings. Cannot
#' be used if missing values are present in test data.
#'
#' @return Returns a list with the following parameters:
#' \item{test_data}{Test data provided by user}
#' \item{EY.draw}{List of posterior samples for the E(Y|X), i.e. \eqn{\beta_0 + \mathbf{V}\boldsymbol{\beta}_{joint} + \sum_{s=1}^q \mathbf{V}_s \boldsymbol{\beta}_s} for each Gibbs sampling iteration.}
#' \item{V.draw}{List of posterior samples for joint scores, \eqn{\mathbf{V}}}
#' \item{U.train}{List of posterior samples for joint loadings for each source, \eqn{\mathbf{U}_s} for \eqn{s=1,\dots,q} given by the training BSFP fit}
#' \item{W.train}{List of posterior samples for individual loadings for each source,  \eqn{\mathbf{W}_s} for \eqn{s=1,\dots,q} given by the training BSFP fit}
#' \item{Vs.draw}{List of posterior samples for individual scores for each source, \eqn{\mathbf{V}_s} for \eqn{s=1,\dots,q}}
#' \item{ranks}{Vector with the estimated joint and individual ranks. \code{ranks[1]} is the estimated joint rank. \code{ranks[2:(q+1)]} correspond to the individual ranks for each source.}
#' \item{tau2.train}{List of posterior samples for the response variance if the response was continuous given by training BSFP fit}
#' \item{beta.train}{List of posterior samples for the regression coefficients used in the predictive model given by training BSFP fit}
#' \item{Xm.draw}{List of posterior samples for missing predictors imputations}
#' 
#' @export
#' 
bsfp.predict.oos <- function(bsfp.fit, test_data, response_type='continuous', model_params = NULL, nsample, progress = TRUE, starting_values = NULL) {
  
  # -------------------------------------------------------------------------- #
  # Extracting the dimensions -----
  # -------------------------------------------------------------------------- #
  
  q <- nrow(test_data) # Number of sources
  p.vec <- apply(test_data, 1, function(source) nrow(source[[1]])) # Number of features per source
  p <- sum(p.vec) # Total number of features
  n <- ncol(test_data[[1,1]]) # Number of subjects
  
  # -------------------------------------------------------------------------- #
  # Determining type of input data -----
  # -------------------------------------------------------------------------- #
  
  # Was the data input as a list?
  if (!("matrix" %in% class(test_data))) {
    
    # Save the number of sources
    q <- length(test_data)
    
    # Initialize new data matrix
    new_data <- matrix(list(), nrow = q, ncol = 1)
    
    # Add in the sources
    for (s in 1:q) {
      new_data[[s,1]] <- test_data[[s]]
    }
    
    # Rename the new version of the data
    test_data <- new_data
  }
  
  # -------------------------------------------------------------------------- #
  # Scaling the data
  # -------------------------------------------------------------------------- #
  for (s in 1:q) {
    test_data[[s,1]] <- test_data[[s,1]]/bsfp.fit$sigma.mat[s,]
  }
  
  # -------------------------------------------------------------------------- #
  # Extracting the model parameters -----
  # -------------------------------------------------------------------------- #
  
  if (is.null(model_params)) {
    model_params <- bsfp.fit$model_params
  }
  
  error_vars <- model_params$error_vars # Error variances
  sigma2_joint <- joint_var <- model_params$joint_var # Variance of joint structure
  sigma2_indiv <- indiv_vars <- model_params$indiv_vars # Variances of individual structure
  beta_vars <- model_params$beta_vars # Variances on betas
  response_vars <- model_params$response_vars; shape <- response_vars[1]; rate <- response_vars[2] # Hyperparameters of variance of response
  
  # -------------------------------------------------------------------------- #
  # Is there a response vector? -----
  # -------------------------------------------------------------------------- #
  
  # response_given <- !is.null(Y_test[[1,1]])
  response_given <- !is.null(bsfp.fit$beta.draw)
  
  # If so, what kind of response is it?
  if(!(response_type %in% c('binary','continuous'))){stop("response_type must be one of 'binary' and 'continuous'")}
  # if (response_given) {
  #   Y_test <- matrix(unlist(Y_test))
  #   
  #   response_type <- if (all(unique(Y_test) %in% c(0, 1, NA))) "binary" else "continuous"
  # }
  
  # -------------------------------------------------------------------------- #
  # Check for missingness in data -----
  # -------------------------------------------------------------------------- #
  
  # Check for missingness
  missingness_in_data <- any(sapply(test_data[,1], function(source) any(is.na(source))))
  
  # Which entries are missing?
  missing_obs <- lapply(test_data[,1], function(source) which(is.na(source)))
  
  # -------------------------------------------------------------------------- #
  # Obtaining the ranks -----
  # -------------------------------------------------------------------------- #
  
  # Saving the ranks from the training data fit
  ranks <- bsfp.fit$ranks
  
  r <- ranks[1]
  r.vec <- ranks[-1]
  
  r_total <- r + sum(r.vec)
  n_beta <- 1 + r_total
  
  # If a response is given, set up the variance matrix for the prior of the betas using the ranks
  if (response_given) {
    Sigma_beta <- matrix(0, nrow = n_beta, ncol = n_beta)
    beta_vars <- c(beta_vars[1], rep(beta_vars[-1], c(r, r.vec)))
    diag(Sigma_beta) <- beta_vars
  }
  
  # -------------------------------------------------------------------------- #
  # Save the loadings, betas, and outcome variance estimated on the training data -----
  # -------------------------------------------------------------------------- #
  
  # Save the joint loadings
  U.train <- bsfp.fit$U.draw
  
  # Save the individual loadings
  W.train <- bsfp.fit$W.draw
  
  # Save the betas
  beta.train <- bsfp.fit$beta.draw
  
  # Save the response variance
  tau2.train <- unlist(bsfp.fit$tau2.draw)
  
  # -------------------------------------------------------------------------- #
  # Storing the posterior samples -----
  # -------------------------------------------------------------------------- #
  
  V.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1))
  Vs.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = q))
  
  if (!response_given) {
    # Z.draw <- VStar.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1))
    VStar.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1))
  }

  if (response_given) {
    # Z.draw <- VStar.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1))
    VStar.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1))
  }
  
  if (missingness_in_data) {
    Xm.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = q, ncol = 1))
  }
  
  # -------------------------------------------------------------------------- #
  # Initialize V, Vs -----
  # -------------------------------------------------------------------------- #
  
  # If no starting values were provided, initialize from priors
  if (is.null(starting_values)) {
    V0 <- matrix(list(), nrow = 1, ncol = 1)
    if (r > 0) {
      V0[[1,1]] <- matrix(rnorm(n*r, mean = 0, sd = sqrt(sigma2_joint)), nrow = n, ncol = r)
    }
    if (r == 0) {
      V0[[1,1]] <- matrix(0, nrow = n, ncol = 1)
    }
    
    Vs0 <- matrix(list(), nrow = 1, ncol = q)
    
    for (s in 1:q) {
      
      # Initialize W and V
      if (r.vec[s] > 0) {
        Vs0[[1,s]] <- matrix(rnorm(n*r.vec[s], mean = 0, sd = sqrt(sigma2_indiv[s])), nrow = n, ncol = r.vec[s])
        
      }
      if (r.vec[s] == 0) {
        Vs0[[1,s]] <- matrix(0, nrow = n, ncol = 1)
        
      }
      
    }
  }
  
  # If starting values were provided, use them as initial values
  if (!is.null(starting_values)) {
    V0 <- starting_values$V
    Vs0 <- starting_values$Vs
  }
  
  if (response_given) {
    # Combining the scores together
    V0.star <- matrix(list(), nrow = 1, ncol = 1)
    if (r > 0) V0.star[[1,1]] <- V0[[1,1]] else V0.star[[1,1]] <- matrix(nrow = n, ncol = r)
    
    Vs0.star <- Vs0
    for (s in 1:q) {
      if (r.vec[s] > 0) Vs0.star[[1,s]] <- Vs0[[1,s]] else Vs0.star[[1,s]] <- matrix(nrow = n, ncol = r.vec[s])
    }
    
    VStar0 <- cbind(1, do.call(cbind, V0.star), do.call(cbind, Vs0.star))
    # Z0 <- matrix(rnorm(n, mean = VStar0 %*% beta.train[[1]][[1,1]], sd = 1))
    
  }
  
  # -------------------------------------------------------------------------- #
  # Storing the initial values -----
  # -------------------------------------------------------------------------- #
  
  V.draw[[1]] <- V0
  Vs.draw[[1]] <- Vs0
  
  if (response_given) {
    # Z.draw[[1]][[1,1]] <- Z0
    VStar.draw[[1]][[1,1]] <- VStar0
  }
  
  if (missingness_in_data) {
  
    Xm0 <- matrix(list(), ncol = 1, nrow = q)
    for (s in 1:q) {
      Xm0[[s,1]] <- matrix(0, nrow = length(missing_obs[[s]]), ncol = 1)
    }
    
    Xm.iter <- Xm0
    Xm.draw[[1]] <- Xm0
  }
  
  # -------------------------------------------------------------------------- #
  # Computing the inverses that don't change from iteration to iteration -----
  # -------------------------------------------------------------------------- #
  
  if (TRUE) {
  # if (!response_given) {
    # Error variance for X.
    SigmaVInv <- diag(rep(1/error_vars, p.vec))
  }
  
  # if (response_given) {
  #   if (response_type == "binary") {
  #     # For V - Combined precisions between data and Z
  #     SigmaVInv <- diag(c(rep(1/error_vars, p.vec), 1))
  #     
  #     # For Vs
  #     SigmaVsInv <- matrix(list(), nrow = q, ncol = q)
  #     
  #     for (s in 1:q) {
  #       SigmaVsInv[[s,s]] <- diag(c(rep(1/error_vars[s], p.vec[s]), 1))
  #     }
  #   }
  #   
  #   # For beta - Combined precisions between intercept and all betas
  #   SigmaBetaInv <- solve(Sigma_beta)
  # }
  
  # -------------------------------------------------------------------------- #
  # Start Gibbs sampling! -----
  # -------------------------------------------------------------------------- #
  
  for (iter in 1:(nsample-1)) {
    if (progress) svMisc::progress(iter/((nsample-1)/100))
    
    # -------------------------------------------------------------------------- #
    # Storing the current values of the parameters -----
    # -------------------------------------------------------------------------- #
    
    V.iter <- V.draw[[iter]]
    U.iter <- U.train[[iter]]
    Vs.iter <- Vs.draw[[iter]]
    W.iter <- W.train[[iter]]
    
    if (response_given) {
      # The current values of the betas
      beta.iter <- beta.train[[iter]][[1,1]]
      
      # Creating a matrix of the joint and individual effects
      beta_indiv.iter <- matrix(list(), nrow = q, ncol = 1)
      
      # Breaking beta down into the intercept,
      beta_intercept.iter <- beta.iter[1,, drop = FALSE]
      
      # Joint effect
      if (r != 0) beta_joint.iter <- beta.iter[2:(r+1),, drop = FALSE] else beta_joint.iter <- matrix(0)
      
      # Individual effects
      if (sum(r.vec) > 0) beta_indiv.iter.temp <- beta.iter[(r+2):n_beta,, drop = FALSE]
      
      for (s in 1:q) {
        # If there is no individual effect
        if (r.vec[s] == 0) beta_indiv.iter[[s, 1]] <- matrix(0)
        
        # If there is an individual effect
        if (r.vec[s] != 0) {
          if (s == 1) beta_indiv.iter[[s, 1]] <- beta_indiv.iter.temp[1:r.vec[s],, drop = FALSE]
          if (s != 1) beta_indiv.iter[[s, 1]] <- beta_indiv.iter.temp[(r.vec[s-1]+1):(r.vec[s-1] + r.vec[s]),, drop = FALSE]
        }
      }
      
      # if (response_type == "binary") {
      #   Z.iter <- Z.draw[[iter]][[1,1]]
      # }
      
      if (response_type == "continuous") {
        tau2.iter <- tau2.train[iter]
      }
      
      # Y_complete <- Y_test
      
    }
    
    X_complete <- test_data
    if (missingness_in_data) {
      # Fill in the completed matrices with the imputed values
      for (s in 1:q) {
        X_complete[[s,1]][missing_obs[[s]]] <- Xm.iter[[s,1]]
      }
    }
    
    # ------------------------------------------------------------------------- #
    # Computing the inverse that changes with tau2 -----
    # ------------------------------------------------------------------------- #
    
    # if (response_given) {
    #   if (response_type == "continuous") {
    #     # For V - Combined error variances between X1, X2, and Y
    #     SigmaVInv <- diag(c(rep(1/error_vars, p.vec), 1/tau2.iter))
    #     
    #     # For Vs
    #     SigmaVsInv <- matrix(list(), nrow = q, ncol = q)
    #     
    #     for (s in 1:q) {
    #       SigmaVsInv[[s,s]] <- diag(c(rep(1/error_vars[s], p.vec[s]), 1/tau2.iter))
    #     }
    #   }
    # }
    
    # ------------------------------------------------------------------------- #
    # Posterior sample for V -----
    # ------------------------------------------------------------------------- #
    
    if (r > 0) {
      if (TRUE) {
      # if (!response_given) {
        # Concatenating Ui's together
        U.iter.combined <- do.call(rbind, U.iter)
        tU_Sigma <- crossprod(U.iter.combined, SigmaVInv)
        
        # Computing the crossprod: t(U.iter) %*% solve(Sigma) %*% U.iter
        tU_Sigma_U <- crossprod(t(tU_Sigma), U.iter.combined)
        
        # The combined centered Xis with the latent response vector
        X.iter <- do.call(rbind, X_complete) - data.rearrange(W.iter)$out %*% do.call(rbind, lapply(Vs.iter, t))
        Bv <- solve(tU_Sigma_U + (1/sigma2_joint) * diag(r))
        
        V.draw[[iter+1]][[1,1]] <- t(matrix(sapply(1:n, function(i) {
          bv <-  tU_Sigma %*% X.iter[,i]
          
          Vi <- MASS::mvrnorm(1, mu = Bv %*% bv, Sigma = Bv)
          Vi
        }), nrow = r))
      }
      
      # if (response_given) {
      #   # Concatenating Ui's together
      #   U.iter.combined <- rbind(do.call(rbind, U.iter), t(beta_joint.iter))
      #   
      #   # Computing the crossprod: t(U.iter) %*% solve(Sigma) %*% U.iter
      #   tU_Sigma <- crossprod(U.iter.combined, SigmaVInv)
      #   tU_Sigma_U <- crossprod(t(tU_Sigma), U.iter.combined)
      #   
      #   Bv <- solve(tU_Sigma_U + (1/sigma2_joint) * diag(r))
      #   
      #   if (response_type == "binary") {
      #     # The combined centered Xis with the latent response vector
      #     X.iter <- rbind(do.call(rbind, X_complete) - data.rearrange(W.iter)$out %*% do.call(rbind, lapply(Vs.iter, t)),
      #                     t(Z.iter - c(beta_intercept.iter) - do.call(cbind, Vs.iter) %*% do.call(rbind, beta_indiv.iter)))
      #   }
      #   
      #   if (response_type == "continuous") {
      #     # The combined centered Xis with the latent response vector
      #     X.iter <- rbind(do.call(rbind, X_complete) - data.rearrange(W.iter)$out %*% do.call(rbind, lapply(Vs.iter, t)),
      #                     t(Y_complete - c(beta_intercept.iter) - do.call(cbind, Vs.iter) %*% do.call(rbind, beta_indiv.iter)))
      #   }
      #   
      #   V.draw[[iter+1]][[1,1]] <- t(matrix(sapply(1:n, function(i) {
      #     bv <- tU_Sigma %*% X.iter[,i]
      #     
      #     Vi <- MASS::mvrnorm(1, mu = Bv %*% bv, Sigma = Bv)
      #     Vi
      #   }), nrow = r))
      # }
      
    }
    
    if (r == 0) {
      V.draw[[iter+1]][[1,1]] <- matrix(0, nrow = n, ncol = 1)
    }
    
    # Updating the value of V
    V.iter <- V.draw[[iter+1]]
    
    # ------------------------------------------------------------------------- #
    # Save the training U, joint loadings -----
    # ------------------------------------------------------------------------- #
    
    U.iter <- U.train[[iter+1]]
    
    # ------------------------------------------------------------------------- #
    # Posterior sample for Vs, s=1,...,q -----
    # ------------------------------------------------------------------------- #
    
    if (TRUE) {
    # if (!response_given) {
      for (s in 1:q) {
        if (r.vec[s] > 0) {
          Xs.iter <- X_complete[[s,1]] - U.iter[[s,1]] %*% t(V.iter[[1,1]])
          Bvs <- solve((1/error_vars[s]) * t(W.iter[[s,s]]) %*% W.iter[[s,s]] + (1/indiv_vars[s]) * diag(r.vec[s]))
          
          Vs.draw[[iter+1]][[1,s]] <- t(matrix(sapply(1:n, function(i) {
            bvs <- (1/error_vars[s]) * t(W.iter[[s,s]]) %*% Xs.iter[, i]
            
            Vsi <- MASS::mvrnorm(1, mu = Bvs %*% bvs, Sigma = Bvs)
            Vsi
          }), nrow = r.vec[s]))
        }
        
        if (r.vec[s] == 0) {
          Vs.draw[[iter+1]][[1,s]] <- matrix(0, nrow = n, ncol = 1)
        }
      }
    }
    
    # if (response_given) {
    #   for (s in 1:q) {
    #     if (r.vec[s] > 0) {
    #       # Combined Ws and beta
    #       W.iter.combined <- rbind(W.iter[[s,s]], t(beta_indiv.iter[[s,1]]))
    #       
    #       tW_Sigma <- crossprod(W.iter.combined, SigmaVsInv[[s,s]])
    #       tW_Sigma_W <- crossprod(t(tW_Sigma), W.iter.combined)
    #       
    #       Bvs <- solve(tW_Sigma_W + (1/indiv_vars[s]) * diag(r.vec[s]))
    #       
    #       if (response_type == "binary") {
    #         # Combined centered Xs and Z
    #         Xs.iter <- rbind(X_complete[[s,1]] - U.iter[[s,1]] %*% t(V.iter[[1,1]]),
    #                          t(Z.iter - c(beta_intercept.iter) - V.iter[[1,1]] %*% beta_joint.iter -
    #                              do.call(cbind, Vs.iter[1, !(1:q %in% s)]) %*% do.call(rbind, beta_indiv.iter[!(1:q %in% s), 1])))
    #       }
    #       
    #       if (response_type == "continuous") {
    #         # Combined centered Xs and Y
    #         Xs.iter <- rbind(X_complete[[s,1]] - U.iter[[s,1]] %*% t(V.iter[[1,1]]),
    #                          t(Y_complete - c(beta_intercept.iter) - V.iter[[1,1]] %*% beta_joint.iter -
    #                              do.call(cbind, Vs.iter[1, !(1:q %in% s)]) %*% do.call(rbind, beta_indiv.iter[!(1:q %in% s), 1])))
    #       }
    #       
    #       Vs.draw[[iter+1]][[1,s]] <- t(matrix(sapply(1:n, function(i) {
    #         bvs <- tW_Sigma %*% Xs.iter[, i]
    #         
    #         Vsi <- MASS::mvrnorm(1, mu = Bvs %*% bvs, Sigma = Bvs)
    #         Vsi
    #       }), nrow = r.vec[s]))
    #     }
    #     
    #     if (r.vec[s] == 0) {
    #       Vs.draw[[iter+1]][[1,s]] <- matrix(0, nrow = n, ncol = 1)
    #     }
    #   }
    # }
    
    # Update the current value of V
    Vs.iter <- Vs.draw[[iter+1]]
    
    # Combine current values of V and V.
    V.iter.star.joint <- V.iter
    if (r == 0) {
      V.iter.star.joint[[1,1]] <- matrix(nrow = n, ncol = r)
    }
    
    Vs.iter.star <- Vs.iter
    for (s in 1:q) {
      if (r.vec[s] == 0) {
        Vs.iter.star[[1,s]] <- matrix(nrow = n, ncol = r.vec[s])
      }
    }
    
    VStar.iter <- cbind(1, do.call(cbind, V.iter.star.joint), do.call(cbind, Vs.iter.star))
    
    # Save the current VStar
    VStar.draw[[iter+1]][[1,1]] <- VStar.iter
    
    # ------------------------------------------------------------------------- #
    # Save training sample for W, individual loadings -----
    # ------------------------------------------------------------------------- #
    
    # Update the current value of W
    W.iter <- W.train[[iter+1]]
    
    # ------------------------------------------------------------------------- #
    # Posterior sample for beta -----
    # ------------------------------------------------------------------------- #
    
    if (response_given) {
      
      # Update the current value of beta
      beta.iter <- beta.train[[iter+1]][[1,1]]
      
      # Creating a matrix of the joint and individual effects
      beta_indiv.iter <- matrix(list(), ncol = 1, nrow = q)
      
      # Breaking beta down into the intercept
      beta_intercept.iter <- beta.iter[1,, drop = FALSE]
      
      # Joint effect
      if (r != 0) beta_joint.iter <- beta.iter[2:(r+1),, drop = FALSE] else beta_joint.iter <- matrix(0)
      
      # Individual effects
      if (sum(r.vec) > 0) beta_indiv.iter.temp <- beta.iter[(r+2):n_beta,, drop = FALSE]
      
      for (s in 1:q) {
        # If there is no individual effect
        if (r.vec[s] == 0) beta_indiv.iter[[s, 1]] <- matrix(0)
        
        # If there is an individual effect
        if (r.vec[s] != 0) {
          if (s == 1) beta_indiv.iter[[s, 1]] <- beta_indiv.iter.temp[1:r.vec[s],, drop = FALSE]
          if (s != 1) beta_indiv.iter[[s, 1]] <- beta_indiv.iter.temp[(r.vec[s-1]+1):(r.vec[s-1] + r.vec[s]),, drop = FALSE]
        }
      }
    }
    
    # ------------------------------------------------------------------------- #
    # Save the training tau2 -----
    # ------------------------------------------------------------------------- #
    
    if (response_given) {
      if (response_type == "continuous") {
        # Update the current value of tau2
        tau2.iter <- tau2.train[iter+1]
      }
    }
    
    # ------------------------------------------------------------------------- #
    # Posterior sample for latent continuous response Z -----
    # ------------------------------------------------------------------------- #
    
    # if (response_given) {
    #   if (response_type == "binary") {
    #     Z.draw[[iter+1]][[1,1]] <- matrix(sapply(1:n, function(i) {
    #       if (Y_complete[i,] == 1) {
    #         truncnorm::rtruncnorm(1, a = 0, mean = (VStar.iter %*% beta.iter)[i,], sd = 1)
    #       } else {
    #         truncnorm::rtruncnorm(1, b = 0, mean = (VStar.iter %*% beta.iter)[i,], sd = 1)
    #       }
    #     }), ncol = 1)
    #   }
    # }
    
    if (missingness_in_data) {
      for (s in 1:q) {
        Es <-  matrix(rnorm(p.vec[s]*n, 0, sqrt(error_vars[s])), nrow = p.vec[s], ncol = n)
        # Xm.draw[[iter+1]][[s,1]] <- matrix((U.iter[[s,1]] %*% t(V.iter[[1,1]]) + W.iter[[s,s]] %*% t(Vs.iter[[1,s]]) + Es)[missing_obs[[s]]])
        Xm.iter[[s,1]] <- matrix((U.iter[[s,1]] %*% t(V.iter[[1,1]]) + W.iter[[s,s]] %*% t(Vs.iter[[1,s]]) + Es)[missing_obs[[s]]])
      }
      
      Xm.draw[[iter]] <- Xm.iter
    }
  }
  
  # -------------------------------------------------------------------------- #
  # Calculating the joint and individual structure, scaled to the data -----
  # -------------------------------------------------------------------------- #
  
  # Calculating the structure for Y at each Gibbs sampling iteration
  EY.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1))
  
  if (missingness_in_data) {    # Scaling the imputed data back to the original data scale 
    for (iter in 1:nsample) {
      for (s in 1:q) {
        Xm.draw[[iter]][[s,1]] <- Xm.draw[[iter]][[s,1]] * bsfp.fit$sigma.mat[s,1]
      }
    }
  }
    
  for (iter in 1:nsample) {
    # Calculate the structure for Y
    if (response_given) {
      if (response_type == "continuous") {
        EY.draw[[iter]][[1,1]] <- VStar.draw[[iter]][[1,1]] %*% beta.train[[iter]][[1,1]]
      }
      
      if (response_type == "binary") {
        EY.draw[[iter]][[1,1]] <- pnorm(VStar.draw[[iter]][[1,1]] %*% beta.train[[iter]][[1,1]])
      }
    }
  }
  
  # Return
  ris_pred <- list(test_data = test_data, # Returning the test data
                   # Y_test = Y_test, # Return the test response vector
                   EY.draw = EY.draw, # Underlying structure
                   V.draw = V.draw, U.train = U.train, W.train = W.train, Vs.draw = Vs.draw, # Components of the structure
                   VStar.draw = VStar.draw, # Components that predict Y,
                   ranks = c(r, r.vec), # Ranks
                   tau2.train = tau2.train, beta.train = beta.train) # Regression parameters
  
  if(missingness_in_data){ris_pred$Xm.draw = Xm.draw}
  
  ris_pred
  
}
