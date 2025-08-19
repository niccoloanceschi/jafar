
rm(list = ls())

# Paths ---- 

run_cluster = T

if(run_cluster){
  path <- '/hpc/group/herringlab/na224/multiomics/' # cluster simulations
  setwd(path)
  source_dir <- 'new_source/'
  if(T){
    ris_dir_jafar = 'ris/sec4_i-cusp/'
    ris_dir_jfr = 'ris/sec4_jfr/'
  } else {
    ris_dir_jafar = 'ris/sec4_i-cusp_linear_resc/'
    ris_dir_jfr = 'ris/sec4_jfr_linear_resc/'
  }
  if(!dir.exists(ris_dir_jafar)){stop("jafar input folder not found", call.=FALSE)}
  if(!dir.exists(ris_dir_jfr)){stop("jfr input folder not found", call.=FALSE)}
} else{ 
  source_dir = '~/Documents/PostDoc/Projects/Multi-Omics/Code/infinitefactor/src/'
  ris_dir = '~/Downloads/jafar_ris_paper/'
}

# MatchAlign Functions ----

library(abind)

Rcpp::sourceCpp(file.path(source_dir,'msf.cpp'))

joint_varimax <- function(L, normalize = TRUE, eps = 1e-05, maxIter=1e3) {
  
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


multiviewRot = function(lambda, eta, theta=NULL, p_views=NULL, normalize=F, eps=1e-5){
  
  nMCMC <- length(lambda)
  iter_print <- nMCMC%/%10
  
  print('Finding Optimal Rotations')
  
  vari <- list()
  print(sprintf(fmt = "%10s%3s%2s", "[",0,"%]"))
  for(t in 1:nMCMC){
    if(is.null(p_views)){
      vari[[t]] <- joint_varimax(lambda[[t]],normalize=normalize,eps=eps)
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
  rotfact = mapply(`%*%`, eta, rots, SIMPLIFY = FALSE)
  if(!is.null(theta)){
    rottheta = mapply(`%*%`, theta, rots, SIMPLIFY = FALSE)
  }
  
  print('Finding Pivot')
  
  norms = sapply(loads, norm, "2")
  piv = loads[order(norms)][[round(length(lambda)/2)]]
  
  print('Matching Pivot')
  
  matches = lapply(loads, msfOUT, piv)
  
  lamout = mapply(aplr, loads, matches, SIMPLIFY = FALSE)
  etaout = mapply(aplr, rotfact, matches, SIMPLIFY = FALSE)
  if(!is.null(theta)){
    thetaout = mapply(aplr, rottheta, matches, SIMPLIFY = FALSE)
  } else {
    thetaout = NULL
  }
  
  print('Postprocessing Complete')
  
  return(list(lambda = lamout, eta = etaout, theta=thetaout))
}

MultiviewMatchAlign <- function(ris_MCMC){
  
  M <- length(ris_MCMC$mu_m)
  p_m <- sapply(ris_MCMC$Lambda_m, ncol)
  nMC <- length(ris_MCMC$K)
  
  K <- dim(ris_MCMC$eta)[3]
  
  idx_p0 <- 1+c(0,cumsum(p_m)[-M])
  idx_pF <- cumsum(p_m)
  
  print("Rotating Shared Component")
  
  # Reshape to lists of matrices
  lambda_lists <- lapply(ris_MCMC$Lambda_m, function(x) apply(x, 1, identity, simplify = FALSE))
  lambda_lists <- lapply(1:nMC, function(t) do.call(rbind, lapply(lambda_lists, `[[`, t)))
  eta_list <- apply(ris_MCMC$eta, 1, identity, simplify = FALSE)
  theta_list <- apply(ris_MCMC$Theta, 1, identity, simplify = FALSE)
  
  # Align and rotate
  LaEtaTh_R <- multiviewRot(lambda_lists,eta_list,theta_list,p_views=p_m,normalize=F)
  rm(lambda_lists,eta_list,theta_list)
  
  # Reshape back to arrays
  lambda_array <- unname(aperm(do.call(abind, c(LaEtaTh_R$lambda, along = 3)),c(3,1,2)))
  ris_MCMC$Lambda_m <- lapply(1:M, function(m) lambda_array[,idx_p0[m]:idx_pF[m],rev(1:K)])
  ris_MCMC$eta = unname(aperm(do.call(abind, c(LaEtaTh_R$eta, along = 3)),c(3,1,2)))
  ris_MCMC$Theta = unname(do.call(abind, c(LaEtaTh_R$theta, along = 1)))
  rm(lambda_array, LaEtaTh_R)
  
  # Revert factor order
  ris_MCMC$eta = ris_MCMC$eta[,,rev(1:K)]
  ris_MCMC$Theta = ris_MCMC$Theta[,rev(1:K)]

  if(grepl('jafar',ris_MCMC$hyper_param$model)){
    
    K_Gm <- sapply(ris_MCMC$phi_m, function(aa) dim(aa)[3])
    
    for(m in 1:M){
      
      print(paste0("Rotating Specific Component - m=",m))  
      
      # Reshape to lists of matrices
      Gamma_m_list <- apply(ris_MCMC$Gamma_m[[m]], 1, identity, simplify = FALSE)
      phi_m_list <- apply(ris_MCMC$phi_m[[m]], 1, identity, simplify = FALSE)
      theta_m_list <- apply(ris_MCMC$Theta_m[[m]], 1, identity, simplify = FALSE)
      
      # Align and rotate
      GaPhiTh_m_R <- multiviewRot(Gamma_m_list, phi_m_list,theta_m_list,normalize=F)
      
      # Reshape back to arrays
      ris_MCMC$Gamma_m[[m]] = unname(aperm(do.call(abind, c(GaPhiTh_m_R$lambda, along = 3)),c(3,1,2)))
      ris_MCMC$phi_m[[m]] = unname(aperm(do.call(abind, c(GaPhiTh_m_R$eta, along = 3)),c(3,1,2)))
      ris_MCMC$Theta_m[[m]] = unname(do.call(abind, c(GaPhiTh_m_R$theta, along = 1)))
      
      # Revert factor order
      ris_MCMC$Gamma_m[[m]] = ris_MCMC$Gamma_m[[m]][,,rev(1:K_Gm[m])]
      ris_MCMC$phi_m[[m]] = ris_MCMC$phi_m[[m]][,,rev(1:K_Gm[m])]
      ris_MCMC$Theta_m[[m]] = ris_MCMC$Theta_m[[m]][,rev(1:K_Gm[m])]
    }
  }
  
  return(ris_MCMC)
}


# Load Ris ----

if(run_cluster){
  if (length(commandArgs(trailingOnly = TRUE)) != 2) {
    stop("Usage: your_script.R <s> <r>")
  }
  ss <- as.integer(commandArgs(trailingOnly = TRUE)[1])
} else {
  ss <- 1
}

## jafar ----

if(T){
  
  pattern <- paste0("StelzerEGA_cv-s", ss, "_y_i-cusp_nMC20000_nBurn15000_nThin10_rep(\\d+)\\.rds")
  
  matching_files <- list.files(path = ris_dir_jafar, pattern = pattern, full.names = F)
  
  if(length(matching_files) >= 1) {
    
    for(filename in matching_files){
      
      print(filename)
      
      risJAFAR <- readRDS(file.path(ris_dir_jafar,filename))
      
      ris_rot = MultiviewMatchAlign(risJAFAR$ris_MCMC)
      
      saveRDS(ris_rot,file.path(ris_dir_jafar,paste0('rotated_',filename)))
      
      rm(risJAFAR,ris_rot)
      
    }
  } else {
    print(' --- No results found for jafar')
  }
}
  
## jfr ----

if(T){
  
  pattern <- paste0("StelzerEGA_cv-s", ss, "_y_jfr_nMC20000_nBurn15000_nThin10_rep(\\d+)\\.rds")
  
  matching_files <- list.files(path = ris_dir_jfr, pattern = pattern, full.names = F)
  
  if(length(matching_files) >= 1) {
    
    for(filename in matching_files){
      
      print(filename)
      
      risJFR <- readRDS(file.path(ris_dir_jfr,filename))
      
      ris_rot = MultiviewMatchAlign(risJFR$ris_MCMC)
      
      saveRDS(ris_rot,file.path(ris_dir_jfr,paste0('rotated_',filename)))
      
      rm(risJFR,ris_rot)
    }
  } else {
    print(' --- No results found for jfr')
  }
}

