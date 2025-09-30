
regular_varimax <- function(L, normalize = FALSE, eps = 1e-05, maxIter=1e3) {
  
  nc <- ncol(L)
  p  <- nrow(L)
  
  if (nc < 2) 
    return(L)
  if (normalize) {
    sc <- sqrt(drop(apply(L, 1L, function(L) sum(L^2))))
    L <- L/sc
  }
  
  R <- diag(nc)
  dQ <- rep(0,maxIter)
  
  for (i in 2:maxIter) {
    
    LR  <- L %*% R
    LR2 <- LR^2
    
    # Q <- t(L/p) %*% (LR2*LR - LR %*% diag(drop(rep(1, p) %*% LR^2))/p) 
    Q <- t(L/p) %*% (LR2*LR - t(colSums(LR2)*t(LR))/p)
    
    svdQ <- La.svd(Q)
    R <- svdQ$u %*% svdQ$vt
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


multiviewAlignment = function(lambda, eta=NULL, theta=NULL, p_views=NULL, normalize=F, eps=1e-5){
  
  nMCMC <- length(lambda)
  iter_print <- nMCMC%/%10
  
  print('Finding Optimal Rotations')
  
  vari <- list()
  print(sprintf(fmt = "%10s%3s%2s", "[",0,"%]"))
  for(t in 1:nMCMC){
    if(is.null(p_views)){
      vari[[t]] <- regular_varimax(lambda[[t]],normalize=normalize,eps=eps)
    } else {
      vari[[t]] <- multivew_varimax(lambda[[t]], p_views=p_views, normalize=normalize,eps=eps)
    }
    if(t %% iter_print == 0){
      print(sprintf(fmt = "%10s%3s%2s", "[",(t%/%iter_print)*10,"%]"))
    }
  }
  
  print('Rotating Variables')
  
  loads = lapply(vari, `[[`, 1)
  rots = lapply(vari, `[[`, 2)
  if(!is.null(eta)){rotfact <- mapply(`%*%`, eta, rots, SIMPLIFY = FALSE)}
  if(!is.null(theta)){rottheta <- mapply(`%*%`, theta, rots, SIMPLIFY = FALSE)}
  
  print('Finding Pivot')
  
  norms = sapply(loads, norm, "2")
  piv = loads[order(norms)][[round(length(lambda)/2)]]
  
  print('Matching Pivot')
  
  matches = lapply(loads, msfOUT, piv)
  
  loads = mapply(aplr, loads, matches, SIMPLIFY = FALSE)
  if(!is.null(eta)){rotfact <- mapply(aplr, rotfact, matches, SIMPLIFY = FALSE)}
  if(!is.null(theta)){rottheta <- mapply(aplr, rottheta, matches, SIMPLIFY = FALSE)}
  
  print('Postprocessing Complete')
  
  output <- list(lambda = loads)
  if(!is.null(eta)){output$eta <- rotfact}
  if(!is.null(theta)){output$theta <- rottheta}
  
  return(output)
}

#' Rotational alignment of latent factors and loading matrices
#' 
#' @description Post-processing routine to solve rotational ambiguity across MCMC samples of latent variables.
#' The alignment is performed using multi-view \code{MatchAlign} on the shared component and regular \code{MatchAlign} on the specific ones.
#'
#' @references 
#' Anceschi N., Ferrari F., Dunson D. B., & Mallick H. (2025). 
#' *Bayesian Joint Additive Factor Models for Multiview Learning*.
#' <https://arxiv.org/abs/2406.00778>
#' 
#' Poworoznek E., Anceschi N., Ferrari F., & Dunson D. B. (2025).
#' *Efficiently Resolving Rotational Ambiguity in Bayesian Matrix Sampling with Matching*.
#' Bayesian Analysis, 1–22. <https://doi.org/10.1214/25-BA1544>
#'
#' @param risMCMC Posterior samples, as returned by \code{\link{gibbs_jafar}} or \code{\link{gibbs_jfr}}.
#'
#' @return A modified version of the input \code{risMCMC}, with latent factors, loading matrices, 
#'   and response loadings (if supervised) rotated according to multi-view \code{MatchAlign}.
#'
#' @importFrom abind abind
#'
#' @export
#' 
multiviewMatchAlign <- function(risMCMC){
  
  if(!'eta'%in%names(risMCMC)){stop("Rotational alignment requires samples of latent quantities. 
                                   Re-run the Gibbs sampler with 'get_latent_vars=TRUE'")}
  
  M <- length(risMCMC$mu_m)
  p_m <- sapply(risMCMC$Lambda_m, ncol)
  
  nMC <- dim(risMCMC$eta)[1]
  K <- dim(risMCMC$eta)[3]
  
  idx_p0 <- 1+c(0,cumsum(p_m)[-M])
  idx_pF <- cumsum(p_m)
  
  is_supervised = ('Theta'%in%names(risMCMC))
  
  print("Rotating Shared Component")
  
  # Reshape to lists of matrices
  lambda_lists <- lapply(risMCMC$Lambda_m, function(x) apply(x, 1, identity, simplify = FALSE))
  lambda_lists <- lapply(1:nMC, function(t) do.call(rbind, lapply(lambda_lists, `[[`, t)))
  eta_list <- apply(risMCMC$eta, 1, identity, simplify = FALSE)
  theta_list = NULL
  if(is_supervised){theta_list <- apply(risMCMC$Theta, 1, identity, simplify = FALSE)}
  
  # Align and rotate
  LaEtaTh_R <- multiviewAlignment(lambda_lists,eta_list,theta_list,p_views=p_m,normalize=F)
  rm(lambda_lists,eta_list,theta_list)
  
  # Reshape back to arrays
  lambda_array <- unname(aperm(do.call(abind, c(LaEtaTh_R$lambda, along = 3)),c(3,1,2)))
  risMCMC$Lambda_m <- lapply(1:M, function(m) lambda_array[,idx_p0[m]:idx_pF[m],])
  risMCMC$eta = unname(aperm(do.call(abind, c(LaEtaTh_R$eta, along = 3)),c(3,1,2)))
  if(is_supervised){risMCMC$Theta = unname(do.call(abind, c(LaEtaTh_R$theta, along = 1)))}
  rm(lambda_array, LaEtaTh_R)
  
  if(grepl('jafar',risMCMC$hyper_param$model)){
    
    K_Gm <- sapply(risMCMC$phi_m, function(aa) dim(aa)[3])
    
    for(m in 1:M){
      
      print(paste0("Rotating Specific Component - m=",m))  
      
      # Reshape to lists of matrices
      Gamma_m_list <- apply(risMCMC$Gamma_m[[m]], 1, identity, simplify = FALSE)
      phi_m_list <- apply(risMCMC$phi_m[[m]], 1, identity, simplify = FALSE)
      theta_m_list <- NULL
      if(is_supervised){
        theta_m_list <- apply(risMCMC$Theta_m[[m]], 1, identity, simplify = FALSE)
      }
      
      # Align and rotate
      GaPhiTh_m_R <- multiviewAlignment(Gamma_m_list, phi_m_list,theta_m_list,normalize=F)
      
      # Reshape back to arrays
      risMCMC$Gamma_m[[m]] = unname(aperm(do.call(abind, c(GaPhiTh_m_R$lambda, along = 3)),c(3,1,2)))
      risMCMC$phi_m[[m]] = unname(aperm(do.call(abind, c(GaPhiTh_m_R$eta, along = 3)),c(3,1,2)))
      if(is_supervised){risMCMC$Theta_m[[m]] = unname(do.call(abind, c(GaPhiTh_m_R$theta, along = 1)))}
    }
  }
  
  return(risMCMC)
}
