
#' Multi-view varimax routine
#' 
#' @param lambda One single sample of the concatenate shared-loadings
#' @param p_views Number of features per view
#' @param normalize Normalization of loading columns (varimax internal parameter)
#' @param eps Threshold on trace
#' @param maxIter Maximum number of iterations
#'
multivew_varimax <- function(L, p_views, normalize = FALSE, eps = 1e-05, maxIter=1e3) {
  
  M <- length(p_views)
  
  idx_0 <- 1+c(0,cumsum(p_views)[-M])
  idx_F <- cumsum(p_views)
  
  pInv <- (1/p_views)[unlist(sapply(1:M,function(m) rep(m,p_views[m])))]
  
  nc <- ncol(L)
  if (nc < 2) 
    return(L)
  if (normalize) {
    sc <- sqrt(drop(apply(L, 1L, function(L) sum(L^2))))
    L <- L/sc
  }
  
  R <- diag(nc)
  dQ <- rep(0,maxIter)
  
  for (i in 2:maxIter) {
    
    LR  <- L%*%R
    LR2 <- LR^2
    
    colSumsLR2 <- do.call(rbind,sapply(1:M,function(m)
      matrix(colSums(LR2[idx_0[m]:idx_F[m],]),p_views[m],nc,byrow=T)))
    
    Q <- t(pInv*L) %*% (LR*LR2 - (pInv*LR) * colSumsLR2)
    
    svdQ <- La.svd(Q)
    R <- svdQ$u%*%svdQ$vt
    dQ[i] <- sum(svdQ$d)
    
    if (dQ[i] < dQ[i-1] * (1 + eps)){break}
  }
  LR <- L %*% R
  if (normalize) 
    LR <- LR * sc
  dimnames(LR) <- dimnames(L)
  class(LR) <- "loadings"
  if(i==maxIter){print("Warning: max iterations reached")}
  list(loadings = LR, rotmat = R, dQ=dQ[2:i])
}

#' Multi-view MatchAlign routine
#' 
#' @param lambda MCMC samples of the concatenate shared-loadings
#' @param theta MCMC samples of the shared latent factors
#' @param risMCMC MCMC samples of the response loadings
#' @param p_views Number of features per view
#' @param normalize Normalization of loading columns (varimax internal parameter)
#'
multiview_MatchAlign = function(lambda, eta, theta, p_views, normalize=F){
  
  nMCMC <- length(lambda)
  iter_print <- nMCMC%/%10
  
  print('Finding Optimal Rotations')
  
  vari <- list()
  print(sprintf(fmt = "%10s%3s%2s", "[",0,"%]"))
  for(t in 1:nMCMC){
    vari[[t]] <- multivew_varimax(lambda[[t]], p_views=p_views, normalize=normalize)
    if(t %% iter_print == 0){
      print(sprintf(fmt = "%10s%3s%2s", "[",(t%/%iter_print)*10,"%]"))
    }
  }
  
  print('Rotating Variables')
  
  loads = lapply(vari, `[[`, 1)
  rots = lapply(vari, `[[`, 2)
  rotfact = mapply(`%*%`, eta, rots, SIMPLIFY = FALSE)
  rottheta = mapply(`%*%`, theta, rots, SIMPLIFY = FALSE)
  
  print('Finding Pivot')
  
  norms = sapply(loads, norm, "2")
  piv = loads[order(norms)][[round(length(lambda)/2)]]
  
  print('Matching Variables Found')
  
  matches = lapply(loads, msfOUT, piv)
  
  lamout = mapply(aplr, loads, matches, SIMPLIFY = FALSE)
  etaout = mapply(aplr, rotfact, matches, SIMPLIFY = FALSE)
  thetaout = mapply(aplr, rottheta, matches, SIMPLIFY = FALSE)
  
  print('Postprocessing Complete')
  
  return(list(lambda = lamout, eta = etaout, Theta=thetaout))
}

#' Postprocessing of shared-component latent variable via multi-view MatchAlign 
#' 
#' @param risMCMC Output of the Gibbs Sampler for JAFAR under the D-CUSP prior
#' @param normalize_col Normalization of loading columns (varimax internal parameter)
#' @param nBurnIn Number of (extra) MCMC iterations to skip
#'
#' @export
#'
postprocess_JAFAR <- function(risMCMC,normalize_col=F,nBurnIn=0){
  
  M <- length(risMCMC$Lambda_m)
  n <- dim(risMCMC$eta)[2]
  p_m <- sapply(risMCMC$Lambda_m, ncol)
  
  idx_p0 <- 1+c(0,cumsum(p_m)[-M])
  idx_pF <- cumsum(p_m)
  
  Kmax <- max(risMCMC$K)
  
  nMC  <- length(risMCMC$K)
  nEff <- nMC-nBurnIn
  
  # Convert to list to abide MatchAling input format
  lambdaSamps <- etaSamps <- thetaSamps <- list()
  
  for(t in 1:nMC){
    
    if(t>nBurnIn){
      
      lambdaSamps[[t-nBurnIn]] <- matrix(0,sum(p_m),Kmax)
      for(m in 1:M){
        lambdaSamps[[t-nBurnIn]][idx_p0[m]:idx_pF[m],] <- risMCMC$Lambda_m[[m]][t,,1:Kmax]
      }
      
      etaSamps[[t-nBurnIn]] <- risMCMC$eta[t,,1:Kmax]
      
      thetaSamps[[t-nBurnIn]] <- risMCMC$Theta[t,1:Kmax]
    }
  }
  
  # Run Multiview MatchAlign
  postprocess_data <- multiview_MatchAlign(lambdaSamps, etaSamps, thetaSamps, p_views=p_m, normalize=normalize_col)
  
  # Convert back to arrays
  Theta    <- array(0.,c(nEff,Kmax))
  eta      <- array(0.,c(nEff,n,Kmax))
  Lambda_m <- lapply(1:M, function(m) array(0.,c(nEff,p_m[m],Kmax)))
  
  for(t in 1:(nEff)){
    Theta[t,] <- postprocess_data$Theta[[t]]
    eta[t,,]  <- postprocess_data$eta[[t]]
    for(m in 1:M){
      Lambda_m[[m]][t,,] <- postprocess_data$lambda[[t]][idx_p0[m]:idx_pF[m],]
    }
  }
  
  return(list(Lambda_m=Lambda_m,Theta=Theta,eta=eta))
}