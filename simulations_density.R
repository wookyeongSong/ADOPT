######################################################################################################
##### Title: ADOPT: Additive Optimal Transport Regression ############################################
##### Description: Simulations on Distribution-valued responses with Wasserstein metric ##############
##### Section: Section 5.1.  #########################################################################
######################################################################################################
library(frechet)
library(fdadensity)
library(fdapace)
library(matrixcalc)
library(parallel)
library(MASS)

##################################################
######### Required Functions (utils) #############
##################################################

## Comparison method (Additive functional regression (ADR), Han et al., 2020)
AddDensReg <- function (Ly, X, x = NULL, hu = NULL, hx = NULL, dSup = NULL) {
  
  n <- nrow(X)
  
  if (is.list(Ly) == FALSE) {
    
    return (message('The response input should be a list of random samples.'))
    
  }
  
  if (is.null(ncol(X))) {
    
    return (message('The design matrix must be multi-dimensional corresponding to additive component.'))
    
  }
  
  if (length(Ly) != nrow(X)) {
    
    return (message('The sample sizes of responses and regressors are different.'))
    
  }
  
  if (is.null(x) == TRUE) {
    
    message('The evaluation grid will be replaced by the observed design matrix.')
    x <- X
    
  } else {
    
    # if (sum(apply(x, 2, diff) < 0) > 0) {
    #   
    #   message('The evaluation grid must be in increasing order. Sorted by increasing for each column.')
    #   x <- apply(x, 2, sort)
    #   
    # }
    
    if (is.null(ncol(x))) {
      
      return (message('The evaluation grid must be multi-dimensional corresponding to additive component.'))
      
    } else if (ncol(X) != ncol(x)) {
      
      return (message('The lengths of columns between X and x are different.'))
      
    }
    
  }
  
  if (is.null(hu)) {
    
    hu <- 0.05
    
  } else if (hu <= 0 || hu >= 1) {
    
    return (message('The bandwidth for smoothing transformed densities should be chosen by a positive number less than 1.'))
    
  }
  
  if (is.null(hx)) {
    
    hx <- rep(0.25 * n^(-1/5), ncol(X)) * apply(apply(x, 2, range), 2, diff)
    
  } else if (length(hx) < 2) {
    
    return (message('The bandwidth must be multi-dimensional.'))
    
  } else if (min(hx) <= 0 || max(hx) >= 1) {
    
    return (message('Bandwidths for additive component functions should be chosen by a positive number less than 1.'))
    
  }
  
  if (is.null(dSup)) {
    
    message('The density response is assumed to have the common support. The observed min and max will be used for the lower/upper limits of the support.') 
    dSup <- range(unlist(Ly))
    
  }
  
  # common arguments
  M <- nrow(x)
  d <- ncol(X)
  
  # minimum bandwidth
  hxMin <- apply(apply(apply(X, 2, sort), 2, diff), 2, max)
  for (j in 1:d) {
    if (hx[j] < hxMin[j]) {
      hx[j] <- hxMin[j]
    }
  }
  
  # common evaluation grids for density and LQD responses
  densGridLen <- 101
  lqdGridLen <- 201
  densGrid <- seq(0, 1, length.out = densGridLen)
  lqdGrid <- seq(0, 1, length.out = lqdGridLen)
  
  # normalize random samples to be supported on [0,1]
  normalizeLy <- mclapply(1:n, 
                          function (i) {
                            normalizeY_i <- (Ly[[i]] - dSup[1]) / diff(dSup)
                          }
  )
  
  # density response reconstruction
  message('Estimating density resopnses...(1/4)')
  Ldens <- mclapply(1:n,
                    function (i) {
                      f_i <- frechet::CreateDensity(y = normalizeLy[[i]], 
                                                    optns = list(outputGrid = densGrid))$y  
                      return (f_i)
                    }
  )
  
  densMat <- matrix(unlist(Ldens), nrow = n, ncol = densGridLen, byrow = TRUE)
  
  # LQD transformation
  message('Transforming density resopnses...(2/4)')
  Llqd <- mclapply(1:n, 
                   function (i) {
                     
                     f_i <- Ldens[[i]]
                     if (min(f_i) < 1e-8) {
                       
                       f_i <- fdadensity::RegulariseByAlpha(densGrid, f_i)
                       
                     }
                     
                     lqd_i <- fdadensity::dens2lqd(dens = f_i, 
                                                   dSup = densGrid,
                                                   lqdSup = lqdGrid)
                     return (lqd_i)
                   }
  )
  
  lqdMat <- matrix(unlist(Llqd), nrow = n, ncol = lqdGridLen, byrow = TRUE)
  
  # smoothing LQD responses
  LlqdSmooth <- mclapply(1:n, 
                         function (i) {
                           lqdSmooth_i <- fdapace::Lwls1D(bw = hu,
                                                          kernel_type = 'epan',
                                                          xin = lqdGrid,
                                                          yin = Llqd[[i]],
                                                          xout = lqdGrid)
                           return (lqdSmooth_i)
                         }
  )
  
  lqdSmoothMat <- matrix(unlist(LlqdSmooth), nrow = n, ncol = lqdGridLen, byrow = TRUE)
  
  # smooth backfitting
  message('Smooth backfitting...(3/4)')
  g0Sbf <- c()
  gjSbf <- mclapply(1:d,
                    function (j) {
                      return (matrix(NA, nrow = lqdGridLen, ncol = M))
                    }
  ) 
  gjSbfAddMean <- mclapply(1:d,
                           function (j) {
                             return (matrix(NA, nrow = lqdGridLen, ncol = M))
                           }
  ) 
  
  
  for (l in 1:lqdGridLen) {
    
    sbfSurf <- fdapace::SBFitting(lqdSmoothMat[,l], x, X, h = hx)
    
    g0Sbf[l] <- sbfSurf$mY
    for (j in 1:d) {
      
      gjSbf[[j]][l,] <- sbfSurf$SBFit[,j]
      gjSbfAddMean[[j]][l,] <- sbfSurf$mY + sbfSurf$SBFit[,j]
      
    }
  }
  
  # LQD inversion to density
  message('Inverting to density resopnses...(4/4)')
  dens0Sbf <- fdadensity::lqd2dens(lqd = g0Sbf,
                                   lqdSup = lqdGrid,
                                   dSup = densGrid)
  densjSbf <- mclapply(1:d, 
                       function (j) {
                         return (matrix(NA, nrow = densGridLen, ncol = M))  
                       }
  )
  
  for (m in 1:M) {
    for (j in 1:d) {
      
      densSbf_j <- fdadensity::lqd2dens(lqd = gjSbfAddMean[[j]][, m],
                                        lqdSup = lqdGrid,
                                        dSup = densGrid)
      
      densjSbf[[j]][, m] <- fdapace::Lwls1D(bw = hu,
                                            kernel_type = 'epan',
                                            xin = densGrid,
                                            yin = densSbf_j,
                                            xout = densGrid)
    }
  }
  
  return (list(Ly = Ly,
               X = X,
               x = x,
               hu = hu,
               hx = hx,
               dSup = dSup,
               lqdGrid = lqdGrid,
               densGrid = (dSup[1] + densGrid*diff(dSup)),
               lqdSbfMean = g0Sbf,
               LlqdSbfComp = gjSbf,
               densSbfMean = dens0Sbf,
               LdensSbfComp = densjSbf
  )
  )
}

## Transport map with density responses with Wasserstein metric
Tf_den_list = function(Qout, Qin){
  
  ## Each column of Qout is the quantile function supported on qSup
  # d is the length of qSup and the number of rows of Qout
  
  if(!is.matrix(Qout)){
    
    stop("Qout should be matrix.")
    
  }
  
  n = dim(Qout)[2]
  
  if(is.vector(Qin)){
    
    Tf <- mclapply(c(1:n), function(i){splinefun(Qin, Qout[,i], method = "natural")})
    #result = sapply(1:nSample, function(i){return(Tf[[i]](qSup, 0))})
    
  } else if(is.matrix(Qin)){
    
    Tf <- mclapply(c(1:n), function(i){splinefun(Qin[,i], Qout[,i], method = "natural")})
    #result = sapply(1:nSample, function(i){return(Tf[[i]](qSup, 0))})
    
  } else{
    
    stop("Qin should be either vector or matrix.")
    
  }
  
  return(Tf)
  
}

## Inverse transport map with density responses with Wasserstein metric
Tf_inv_den_list = function(Qout, Qin){
  
  ## Each column of Qout is the quantile function supported on qSup
  # d is the length of qSup and the number of rows of Qout
  
  n = dim(Qout)[2]
  
  if(!is.matrix(Qout)){
    
    stop("Qout should be matrix.")
    
  }
  
  if(is.vector(Qin)){
    
    ITf <- mclapply(c(1:n), function(i){splinefun(Qout[,i], Qin, method = "natural")})
    #result = sapply(1:nSample, function(i){return(ITf[[i]](qSup, 0))})
    
  } else if(is.matrix(Qin)){
    
    ITf <- mclapply(c(1:n), function(i){splinefun(Qout[,i], Qin[,i], method = "natural")})
    #result = sapply(1:nSample, function(i){return(ITf[[i]](qSup, 0))})
    
  } else{
    
    stop("Qin should be either vector or matrix.")
    
  }
  
  return(ITf)
  
}

## ADOPT fitting and predict 
ADOPT_den = function(xin = NULL, qin = NULL, xout = NULL, max_iter = 10, tol = 1e-5, qSup = seq(0,1,length.out = 201)){
  
  n = nrow(xin)
  p = ncol(xin)
  l = length(qSup)
  
  lower_bdd = min(qin)
  upper_bdd = max(qin)
  
  if(is.null(xout)){
    xout = xin
  }
  xout = as.matrix(xout)
  
  
  Fmean = rowMeans(qin)
  tf_qin = Tf_den_list(qin,Fmean)
  
  ## Initialize geodesics (p = length(init_tf) = length(init_tf_inv))
  id = matrix(rep(qSup, n), l, n)
  tf = mclapply(1:p,function(j) {return(Tf_den_list(id,id))})
  tf_inv = mclapply(1:p,function(j) {return(Tf_inv_den_list(id,id))})
  
  ## Estimated
  tf_out = mclapply(1:p,function(j) {return(Tf_den_list(id,id))})
  
  ## Initialize fitted values
  prev = init_fitted = matrix(rep(Fmean, n), ncol = n)
  loss = Inf
  iter = 0
  
  while(iter < max_iter & loss > tol){
    
    for(j in 1:p){
      
      ## Obtain p-th partial residual
      curr = init_fitted
      
      if(j > 1){
        
        for(k in (j-1):1){
          
          curr = sapply(1:n, function(i) {tf_inv[[k]][[i]](curr[,i],0)})

        }
        
      }
      
      curr = sapply(1:n, function(i) {tf_qin[[i]](curr[,i],0)})

      if(j < p){
        
        for(k in p:(j+1)){
          
          curr = sapply(1:n, function(i) {tf_inv[[k]][[i]](curr[,i],0)})

        }
        
      }
      
      ## Sanity Check
      curr[curr < lower_bdd] = lower_bdd
      curr[curr > upper_bdd] = upper_bdd
      curr = apply(curr, 2, sort)
      
      ## Apply local-Freg to p-th partial residual
      g_j = t(LocDenReg(xin = xin[,j], qin = t(curr), xout = xin[,j],  optns = list(bw = 0.1,qSup = qSup))$qout)
      g_out = t(LocDenReg(xin = xin[,j], qin = t(curr), xout = xout[,j],  optns = list(bw = 0.1,qSup = qSup))$qout)
      
      ## F-mean of p-th additive term
      g_j_Fmean = rowMeans(g_j)
      
      ## Update tf and inv
      tf[[j]] = Tf_den_list(g_j,g_j_Fmean)
      tf_inv[[j]] = Tf_inv_den_list(g_j,g_j_Fmean)
      
      tf_out[[j]] = Tf_den_list(g_out,g_j_Fmean)
      
    }
    
    res_fitted = init_fitted
    for(j in 1:p){
      
      res_fitted = sapply(1:n, function(i) {tf[[j]][[i]](res_fitted[,i],0)})
      
    }
    
    ## Sanity Check
    res_fitted[res_fitted < lower_bdd] = lower_bdd
    res_fitted[res_fitted > upper_bdd] = upper_bdd
    res_fitted = apply(res_fitted, 2, sort)
    
    loss = sum((res_fitted-prev)^2)/((l-1)*n)
    iter = iter + 1
    prev = res_fitted
    
    print(iter)
    print(loss)
    
  }
  
  res_out = matrix(rep(Fmean, nrow(xout)), ncol = nrow(xout))
  for(j in 1:p){
    
    res_out = sapply(1:nrow(xout), function(i) {tf_out[[j]][[i]](res_out[,i],0)})
    
  }
  
  ## Sanity Check
  res_out[res_out < lower_bdd] = lower_bdd
  res_out[res_out > upper_bdd] = upper_bdd
  res_out = apply(res_out, 2, sort)
  
  return(list(xin = xin, qfitted = res_fitted, qout = res_out))
  
}

## ADOPT fitting only (faster than ADOPT_den)
ADOPT_den_fit = function(xin = NULL, qin = NULL, xout = NULL, max_iter = 10, tol = 1e-5, qSup = seq(0,1,length.out = 201)){
  
  n = nrow(xin)
  p = ncol(xin)
  l = length(qSup)
  
  lower_bdd = min(qin)
  upper_bdd = max(qin)
  
  if(is.null(xout)){
    xout = xin
  }
  xout = as.matrix(xout)
  
  
  Fmean = rowMeans(qin)
  tf_qin = Tf_den_list(qin,Fmean)
  
  ## Initialize geodesics (p = length(init_tf) = length(init_tf_inv))
  id = matrix(rep(qSup, n), l, n)
  tf = mclapply(1:p,function(j) {return(Tf_den_list(id,id))})
  tf_inv = mclapply(1:p,function(j) {return(Tf_inv_den_list(id,id))})
  
  ## Estimated
  #tf_out = mclapply(1:p,function(j) {return(Tf_den_list(id,id))})
  
  ## Initialize fitted values
  prev = init_fitted = matrix(rep(Fmean, n), ncol = n)
  loss = Inf
  iter = 0
  
  while(iter < max_iter & loss > tol){
    
    for(j in 1:p){
      
      ## Obtain p-th partial residual
      curr = init_fitted
      
      if(j > 1){
        
        for(k in (j-1):1){
          
          curr = sapply(1:n, function(i) {tf_inv[[k]][[i]](curr[,i],0)})
          
        }
        
      }
      
      curr = sapply(1:n, function(i) {tf_qin[[i]](curr[,i],0)})
      
      if(j < p){
        
        for(k in p:(j+1)){
          
          curr = sapply(1:n, function(i) {tf_inv[[k]][[i]](curr[,i],0)})
          
        }
        
      }
      
      ## Sanity Check
      curr[curr < lower_bdd] = lower_bdd
      curr[curr > upper_bdd] = upper_bdd
      curr = apply(curr, 2, sort)
      
      ## Apply local-Freg to p-th partial residual
      g_j = t(LocDenReg(xin = xin[,j], qin = t(curr), xout = xin[,j],  optns = list(bw = 0.1,qSup = qSup))$qout)
      #g_out = t(LocDenReg(xin = xin[,j], qin = t(curr), xout = xout[,j],  optns = list(bw = 0.1,qSup = qSup))$qout)
      
      ## F-mean of p-th additive term
      g_j_Fmean = rowMeans(g_j)
      
      ## Update tf and inv
      tf[[j]] = Tf_den_list(g_j,g_j_Fmean)
      tf_inv[[j]] = Tf_inv_den_list(g_j,g_j_Fmean)
      
      #tf_out[[j]] = Tf_den_list(g_out,g_j_Fmean)
      
    }
    
    res_fitted = init_fitted
    for(j in 1:p){
      
      res_fitted = sapply(1:n, function(i) {tf[[j]][[i]](res_fitted[,i],0)})
      
    }
    
    ## Sanity Check
    res_fitted[res_fitted < lower_bdd] = lower_bdd
    res_fitted[res_fitted > upper_bdd] = upper_bdd
    res_fitted = apply(res_fitted, 2, sort)
    
    loss = sum((res_fitted-prev)^2)/((l-1)*n)
    iter = iter + 1
    prev = res_fitted
    
    print(iter)
    print(loss)
    
  }
  
  #res_out = matrix(rep(Fmean, nrow(xout)), ncol = nrow(xout))
  #for(j in 1:p){
  
  #  res_out = sapply(1:nrow(xout), function(i) {tf_out[[j]][[i]](res_out[,i],0)})
  
  #}
  
  ## Sanity Check
  #res_out[res_out < lower_bdd] = lower_bdd
  #res_out[res_out > upper_bdd] = upper_bdd
  #res_out = apply(res_out, 2, sort)
  
  return(list(xin = xin, qfitted = res_fitted))
  #return(list(xin = xin, qfitted = res_fitted, qout = res_out))
  
}


###################################################################################
######### Simulation Data Generation analyzed in Section 5.1 (datgen) #############
###################################################################################

## Case I (ADOPT True Model, 3d Predictors, beta distribution)
ADOPT_DatGenDen1 = function(n = 100, qSup = seq(0,1,length.out=201)){
  
  ## Predictors X = (X1, X2, X3)
  # Define parameters
  mu <- c(0, 0, 0)  # Mean vector
  Sigma <- matrix(c(1, 0.5, 0.3, 0.5, 1, 0.5, 0.3, 0.5, 1), nrow = 3)  # Covariance matrix
  
  # Generate bivariate normal samples size of n
  V <- mvrnorm(n, mu, Sigma)
  
  # Apply the standard normal CDF (pnorm) to each element
  X <- pnorm(V)
  
  # Generate Response
  l = length(qSup)
  unif_qf = seq(0,1,length.out= l)
  
  g_1 = sapply(1:n, function(i) {qbeta(unif_qf,1+2*X[i,1],1)})
  g_1_Fmean = rowMeans(g_1)
  
  g_2 = sapply(1:n, function(i) {qbeta(unif_qf,1,2+3*X[i,2])})
  g_2_Fmean = rowMeans(g_2)
  
  g_3 = sapply(1:n, function(i) {qbeta(unif_qf,1/2+1/2*X[i,3],1/2+1/2*X[i,3])})
  g_3_Fmean = rowMeans(g_3)
  
  e = runif(n,min=-1,max=1)
  GenNoise <- function(u,e) u + 1/(2*pi)*e*sin(2*pi*u)
  
  tf = list()
  tf[[1]] = Tf_den_list(g_1,g_1_Fmean)
  tf[[2]] = Tf_den_list(g_2,g_2_Fmean)
  tf[[3]] = Tf_den_list(g_3,g_3_Fmean)
  
  Y = matrix(rep(unif_qf, n), length(unif_qf), n)
  for(j in 1:3){
    
    Y = sapply(1:n,function(i){tf[[j]][[i]](Y[,i],0)})
    
  }
  qin = sapply(1:n, function(i){GenNoise(Y[,i],e[i])})
  
  Ly <- list()
  for (i in 1:n) {
    
    Ly[[i]] <- qin[,i]
 
  }
  
  return(list(xin = X, qin = qin, Y_true = Y, Ly = Ly))
}

## Case II (ADOPT True Model, 3d Predictors, normal distribution)
ADOPT_DatGenDen2 = function(n = 100, qSup = seq(0,1,length.out=201)){
  
  ## Predictors X = (X1, X2, X3)
  # Define parameters
  mu <- c(0, 0, 0)  # Mean vector
  Sigma <- matrix(c(1, 0.5, 0.3, 0.5, 1, 0.5, 0.3, 0.5, 1), nrow = 3)  # Covariance matrix
  
  # Generate bivariate normal samples size of n
  V <- mvrnorm(n, mu, Sigma)
  
  # Apply the standard normal CDF (pnorm) to each element
  X <- pnorm(V)
  
  # Generate Response
  l = length(qSup)
  unif_qf = seq(0,1,length.out= l)
  
  g_1 = sapply(1:n, function(i) {qnorm(c(0.001, unif_qf[2:(l-1)],0.999), mean = X[i,1], 1)})
  g_1_Fmean = rowMeans(g_1)
  
  g_2 = sapply(1:n, function(i) {qnorm(c(0.001, unif_qf[2:(l-1)],0.999), mean = X[i,2]^2, 1)})
  g_2_Fmean = rowMeans(g_2)
  
  g_3 = sapply(1:n, function(i) {qnorm(c(0.001, unif_qf[2:(l-1)],0.999), mean = exp(-X[i,3]), 1)})
  g_3_Fmean = rowMeans(g_3)
  
  e = runif(n,min=-1,max=1)
  GenNoise <- function(u,e) u + 1/(2*pi)*e*sin(2*pi*u)
  
  tf = list()
  tf[[1]] = Tf_den_list(g_1,g_1_Fmean)
  tf[[2]] = Tf_den_list(g_2,g_2_Fmean)
  tf[[3]] = Tf_den_list(g_3,g_3_Fmean)
  
  Y = matrix(rep(unif_qf, n), length(unif_qf), n)
  for(j in 1:3){
    
    Y = sapply(1:n,function(i){tf[[j]][[i]](Y[,i],0)})
    
  }
  qin = sapply(1:n, function(i){GenNoise(Y[,i],e[i])})
  
  Ly <- list()
  for (i in 1:n) {
    
    Ly[[i]] <- qin[,i]

  }
  
  return(list(xin = X, qin = qin, Y_true = Y, Ly = Ly))
}


####################################################################
######### Simulations in Section 5.1 (simulations_den) #############
####################################################################

## Set number of cores
num_cores <- detectCores() - 1 


## quantile dSup
qSup = seq(0,1,length.out=201)
l = length(qSup)


## Conduct simulations
run_experiment_den = function(n, dat){
  
  ### ADOPT
  estADOPT = ADOPT_den_fit(xin = dat$xin, qin = dat$qin, max_iter = 5, tol = 1e-5, qSup = qSup)
  #estADOPT = ADOPT_den(xin = dat$xin, qin = dat$qin, max_iter = 1, tol = 1e-5, qSup = qSup)
  MISE_den_ADOPT = sum( (estADOPT$qfitted - dat$Y_true)^2 ) / (n * l)
  
  ### Global F Reg
  fitGloDensReg = GloDenReg(xin = dat$xin, qin = t(dat$qin), optns = list(qSup = qSup))
  estGloDensReg = t(fitGloDensReg$qout)
  MISE_den_GF =  sum( (estGloDensReg - dat$Y_true)^2 ) / (n * l)
  
  ### Additive Den Reg
  estAddDensReg <- AddDensReg(Ly = dat$Ly, X = dat$xin, dSup = c(0,1))
  estAddDensReg_lqd = matrix(0, nrow = l , ncol = n)
  for(i in 1:n){
    estAddDensReg_lqd[,i] = estAddDensReg$lqdSbfMean 
    for(j in 1:ncol(dat$xin)){
      estAddDensReg_lqd[,i] = estAddDensReg_lqd[,i] + estAddDensReg$LlqdSbfComp[[j]][,i]  
    }
  }
  
  QuanEstAdd = apply(exp(estAddDensReg_lqd),MARGIN = 2,FUN = cumsum) %*% diag(1/colSums(exp(estAddDensReg_lqd)))
  MISE_den_ADR = sum( ( QuanEstAdd - dat$Y_true)^2 ) / (n * l)
  
  # F mean
  MISE_Fmean =sum((dat$qin- matrix(rep(rowMeans(dat$qin), n), ncol = n))^2)  / (n * l)
  
  return(c(MISE_den_ADOPT, MISE_den_GF, MISE_den_ADR, MISE_Fmean))
  
}

generate_df_den <- function(ADOPT_func, B, num_cores, nvals = c(100, 300)) {
  out_list <- vector("list", length(nvals))
  
  for (i in seq_along(nvals)) {
    n_val   <- nvals[i]
    
    # Run in parallel for B replicates
    results <- mclapply(1:B, ADOPT_func, n = n_val, mc.cores = num_cores)
    
    # Combine into a matrix
    results_mat <- do.call(rbind, results)
    
    # Build a data frame with columns: pred_err, Method, N
    df_temp <- data.frame(
      pred_err = c(results_mat),
      Method   = rep(c("ADOPT","GFR","ADR","FM"), each = B),
      N        = factor(n_val, levels = nvals)
    )
    out_list[[i]] <- df_temp
  }
  
  # Merge all n-values into one data frame
  df_all <- do.call(rbind, out_list)
  
  # Ensure Method is a factor in desired order
  df_all$Method <- factor(df_all$Method, levels = c("ADOPT","GFR","ADR","FM"))
  
  return(df_all)
}

## Case I
ADOPT_den1 <- function(n, b) {
  
  set.seed(b)  # Set seed for reproducibility
  dat_den = ADOPT_DatGenDen1(n = n, qSup = qSup)
  
  res = run_experiment_den(n = n, dat = dat_den)
  
  return(res)
}


## Case II
ADOPT_den2 <- function(n, b) {
  
  set.seed(b)  # Set seed for reproducibility
  dat_den = ADOPT_DatGenDen2(n = n, qSup = qSup)
  
  res = run_experiment_den(n = n, dat = dat_den)
  
  return(res)
}


## Generate Table 1
summarize_blocks <- function(x, blocks, scale = 1000, effective_n = NULL) {
  do.call(rbind, lapply(names(blocks), function(lbl) {
    idx <- blocks[[lbl]]
    vals <- x[idx]
    n_se <- if (is.null(effective_n)) length(idx) else effective_n
    B = 200
    data.frame(
      Method      = lbl,
      n    = n_se,
      Mean_x1000  = mean(vals) * scale,
      SE_x1000    = (sd(vals) / sqrt(B)) * scale,
      check.names = FALSE,
      row.names   = NULL
    )
  }))
}


res_ADOPT_den1 = generate_df_den(ADOPT_den1, B = 200, num_cores = num_cores, nvals = c(100, 300))
res_ADOPT_den2 = generate_df_den(ADOPT_den2, B = 200, num_cores = num_cores, nvals = c(100, 300))

# Blocks for the two regimes
blocks_100 <- list(ADOPT = 1:200,   GF = 201:400, ADR = 401:600, FM = 601:800)
blocks_300 <- list(ADOPT = 801:1000, GF = 1001:1200, ADR = 1201:1400, FM = 1401:1600)

# Build the two tables (using your intended n for SE: 100 and 300)
summary_100_case1 <- summarize_blocks(res_ADOPT_den1[, 1], blocks_100, scale = 1000, effective_n = 100)
summary_300_case1 <- summarize_blocks(res_ADOPT_den1[, 1], blocks_300, scale = 1000, effective_n = 300)

# Build the two tables (using your intended n for SE: 100 and 300)
summary_100_case2 <- summarize_blocks(res_ADOPT_den2[, 1], blocks_100, scale = 1000, effective_n = 100)
summary_300_case2 <- summarize_blocks(res_ADOPT_den2[, 1], blocks_300, scale = 1000, effective_n = 300)

# Combine
summary_tbl_case1 <- rbind(summary_100_case1, summary_300_case1)
summary_tbl_case2 <- rbind(summary_100_case2, summary_300_case2)

# Optional rounding for presentation
summary_tbl_case1$Mean_x1000 <- round(summary_tbl_case1$Mean_x1000, 3)
summary_tbl_case1$SE_x1000   <- round(summary_tbl_case1$SE_x1000, 3)
summary_tbl_case1$Case = "I"

summary_tbl_case2$Mean_x1000 <- round(summary_tbl_case2$Mean_x1000, 3)
summary_tbl_case2$SE_x1000   <- round(summary_tbl_case2$SE_x1000, 3)
summary_tbl_case2$Case = "II"

summary_tbl_den = rbind(summary_tbl_case1,summary_tbl_case2)
print(summary_tbl_den)
