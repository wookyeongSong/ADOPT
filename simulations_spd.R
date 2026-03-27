######################################################################################################
##### Title: ADOPT: Additive Optimal Transport Regression ############################################
##### Description: Simulations on SPD matrices-valued response with log-Cholesky metric ##############
##### Section: Section 5.2.  #########################################################################
######################################################################################################
library(frechet)
library(fdadensity)
library(fdapace)
library(matrixcalc)
library(parallel)
library(MASS)
library(Matrix)


##################################################
######### Required Functions (utils) #############
##################################################

## Transport map with SPD responses with Log-Cholesky metric
Tf_spd_list = function(Mout, Min, corrOut){
  
  ## Min and Mout should be 3d array size of (m x m x n)
  
  n = dim(Mout)[3]
  
  if(!all(apply(Mout,3,is.positive.semi.definite))){
    
    stop("output matrix is not SPD.")
    
  }
  
  if(length(dim(Min))==2){
    
    if(!is.positive.semi.definite(Min)){
      
      stop("input matrix is not SPD.")
      
    }
    
    Tf = lapply(1:n, function(i){logcholfun(Min, Mout[,,i], corrOut)})
    
  } else if(length(dim(Min))==3){
    
    if(!all(apply(Min,3,is.positive.semi.definite))){
      
      stop("input matrix is not SPD.")
      
    }
    
    Tf = lapply(1:n, function(i){logcholfun(Min[,,i], Mout[,,i], corrOut)})
    
  } else{
    
    stop("input should be either matrix or 3d array.")
    
  }
  
  return(Tf)
  
}

## Inverse transport map with SPD responses with Log-Cholesky metric
Tf_inv_spd_list = function(Mout, Min, corrOut){
  
  ## Min and Mout should be 3d array size of (m x m x n)
  
  n = dim(Mout)[3]
  
  if(!all(apply(Mout,3,is.positive.semi.definite))){
    
    stop("output matrix is not SPD.")
    
  }
  
  if(length(dim(Min))==2){
    
    if(!is.positive.semi.definite(Min)){
      
      stop("input matrix is not SPD.")
      
    }
    
    ITf = lapply(1:n, function(i){logcholfun(Mout[,,i], Min, corrOut)})
    
  } else if(length(dim(Min))==3){
    
    if(!all(apply(Min,3,is.positive.semi.definite))){
      
      stop("input matrix is not SPD.")
      
    }
    
    ITf = lapply(1:n, function(i){logcholfun(Mout[,,i],Min[,,i], corrOut)})
    
  } else{
    
    stop("input should be either matrix or 3d array.")
    
  }  
  
  return(ITf)
  
}

## Geodesic optimal transport with Log-Cholesky metric 
logcholfun = function(Min, Mout, corrOut){
  
  lower.part = function(M){
    
    lower.triangle(M) - diag(diag(M))
    
  }
  
  Min = as.matrix(Matrix::nearPD(Min, corr = corrOut)$mat)
  Min = as.matrix(Matrix::forceSymmetric(Min))
  
  Mout = as.matrix(Matrix::nearPD(Mout, corr = corrOut)$mat)
  Mout = as.matrix(Matrix::forceSymmetric(Mout))
  
  res_fun = function(M){
    
    ## Geodesic from Min to Mout
    C1 = t(chol(Min))
    C2 = t(chol(Mout))
    
    L1 = lower.part(C1)
    L2 = lower.part(C2)
    
    D1 = diag(C1)
    D2 = diag(C2)
    
    geod = (L2 - L1) + diag(exp(log(D2) - log(D1)))
    
    M = as.matrix(Matrix::nearPD(M, corr = corrOut)$mat)
    M = as.matrix(Matrix::forceSymmetric(M))
    
    input = t(chol(M))
    res_lower_tri = lower.part(geod) + lower.part(input) + diag(exp(log(diag(geod)) + log(diag(input))))
    
    res = res_lower_tri %*% t(res_lower_tri)
    
    return(res)
    
  }
  
  return(res_fun)
  
}

## Loss functions with Log-Cholesky metric
loss_SPD = function(M1, M2){
  
  n = dim(M1)[3]
  res = 0
  
  for(i in 1:n){
    
    C1 = t(chol(M1[,,i]))
    C2 = t(chol(M2[,,i]))
    
    res = res + sqrt(sum((C1- C2)^2))
  }
  
  return(res)
  
}

## ADOPT fitting and predict 
ADOPT_spd = function(xin = NULL, Min = NULL, xout = NULL, max_iter = 10, tol = 1e-3, bw = 0.1, corrOut= FALSE){
  
  n = nrow(xin)
  p = ncol(xin)
  m = dim(Min)[1] ## dimension of SPD
  
  if(is.null(xout)){
    xout = xin
  }
  
  ## F mean of response
  Fmean = CovFMean(M=Min,optns=list(metric="log_cholesky",corrOut=corrOut))$Mout[[1]]
  tf_Min = Tf_spd_list(Min, Fmean, corrOut = corrOut)
  
  ## Initialization of geodesic transport
  id = array(diag(m),dim=c(m,m,n))
  tf = lapply(1:p, function(j){Tf_spd_list(id, id, corrOut = corrOut)})
  tf_inv = lapply(1:p, function(j){Tf_inv_spd_list(id, id, corrOut = corrOut)})
  
  ## Estimated
  id_out = array(diag(m),dim=c(m,m,nrow(xout)))
  tf_out = lapply(1:p, function(j){Tf_spd_list(id_out, id_out, corrOut = corrOut)})
  
  
  ## Initialize fitted values
  prev = init_fitted = array(Fmean, dim=c(m,m,n))
  loss = Inf
  iter = 0
  
  while(iter < max_iter & loss > tol){
    
    for(j in 1:p){
      
      curr = init_fitted
      
      if(j > 1){
        for(k in (j-1):1){
          
          curr = array(sapply(1:n, function(i) {tf_inv[[k]][[i]](curr[,,i])}),dim=c(m,m,n))
          
        }
      }
      
      curr = array(sapply(1:n, function(i) {tf_Min[[i]](curr[,,i])}), dim=c(m,m,n))
      
      if(j < p){
        for(k in p:(j+1)){
          
          curr = array(sapply(1:n, function(i) {tf_inv[[k]][[i]](curr[,,i])}),dim=c(m,m,n))
          
        }
      }
      
      ## Sanity Check
      for(i in 1:n){
        
        curr[,,i] = as.matrix(Matrix::nearPD(curr[,,i],corr = corrOut)$mat)
        curr[,,i] = as.matrix(Matrix::forceSymmetric(curr[,,i]))
        
      }
      
      ## Apply Backfitting (using local Frechet regression) to j-th partial residual
      g_j = LocCovReg(x = xin[,j] ,M=curr, xout = xin[,j],optns=list(corrOut=corrOut, bwCov = bw, metric="log_cholesky"))
      g_j = simplify2array(g_j$Mout)
      
      g_out = LocCovReg(x = xin[,j] ,M=curr, xout = xout[,j],optns=list(corrOut=corrOut, bwCov = bw, metric="log_cholesky"))
      g_out = simplify2array(g_out$Mout)
      
      ## F-mean of j-th additive term
      g_j_Fmean = CovFMean(M=g_j, optns=list(corrOut=corrOut,metric="log_cholesky"))$Mout[[1]]
      
      ## Update tf and tf_inv
      tf[[j]] = Tf_spd_list(g_j, g_j_Fmean, corrOut = corrOut)
      tf_inv[[j]] = Tf_inv_spd_list(g_j, g_j_Fmean, corrOut = corrOut)
      
      tf_out[[j]] = Tf_spd_list(g_out, g_j_Fmean, corrOut = corrOut)
      
    }
    
    res_fitted = init_fitted
    for(j in 1:p){
      
      res_fitted = array(sapply(1:n, function(i) {tf[[j]][[i]](res_fitted[,,i])}),dim=c(m,m,n))
      
    }
    
    ## Sanity Check
    for(i in 1:n){
      
      res_fitted[,,i] = as.matrix(Matrix::nearPD(res_fitted[,,i],corr = corrOut)$mat)
      res_fitted[,,i] = as.matrix(Matrix::forceSymmetric(res_fitted[,,i]))
      
    }
    
    loss = loss_SPD(res_fitted,prev)/n
    iter = iter + 1
    prev = res_fitted
    
    print(iter)
    print(loss)
    
  }
  
  res_out = array(Fmean, dim=c(m,m,nrow(xout)))
  for(j in 1:p){
    
    res_out = array(sapply(1:nrow(xout), function(i) {tf_out[[j]][[i]](res_out[,,i])}), dim=c(m,m,nrow(xout)))
    
  }
  
  ## Sanity Check
  for(i in 1:nrow(xout)){
    
    res_out[,,i] = as.matrix(Matrix::nearPD(res_out[,,i],corr = corrOut)$mat)
    res_out[,,i] = as.matrix(Matrix::forceSymmetric(res_out[,,i]))
    
  }
  
  #return(list(xin = xin, qfitted = res_fitted))
  return(list(xin = xin, Mfitted = res_fitted, Mout = res_out))
  
}


###################################################################################
######### Simulation Data Generation analyzed in Section 5.2 (datgen) #############
###################################################################################

## Case I (ADOPT True Model, 3d Predictors, linear model)
ADOPT_DatGenSPD1 = function(n = 100, m = 5, corrOut = FALSE){
  
  lower.part = function(M){
    
    lower.triangle(M) - diag(diag(M))
    
  }
  
  ## Predictors X = (X1, X2, X3)
  # Define parameters
  mu <- c(0, 0, 0)  # Mean vector
  Sigma <- matrix(c(1, 0.5, 0.3, 0.5, 1, 0.5, 0.3, 0.5, 1), nrow = 3)  # Covariance matrix
  
  # Generate bivariate normal samples size of n
  V <- mvrnorm(n, mu, Sigma)
  
  # Apply the standard normal CDF (pnorm) to each element
  X <- pnorm(V)
  
  # Generate Response
  g_1= lapply(1:n, function(i){
    
    L = matrix(0,nrow=m, ncol=m)
    for(row in 1:m){
      for(col in 1:row){
        L[row,col] = exp(-abs(row-col)/2)*(1/2 + 1/2*X[i,1])
      }
    }
    L %*% t(L)
    
  })
  g_1 = simplify2array(g_1)
  
  g_2= lapply(1:n, function(i){
    
    L = matrix(0,nrow=m, ncol=m)
    for(row in 1:m){
      for(col in 1:row){
        L[row,col] = exp(-abs(row-col)/3)*(1/4 + 1/4*X[i,2])
      }
    }
    L %*% t(L)
    
  })
  g_2 = simplify2array(g_2)
  
  
  g_3= lapply(1:n, function(i){
    
    L = matrix(0,nrow=m, ncol=m)
    for(row in 1:m){
      for(col in 1:row){
        L[row,col] = exp(-abs(row-col)/4)*(1/8 + 1/8*X[i,3])
      }
    }
    L %*% t(L)
    
  })
  g_3 = simplify2array(g_3)
  
  
  g_1_Fmean = CovFMean(M=g_1,optns=list(corrOut = corrOut, metric="log_cholesky"))$Mout[[1]]
  g_2_Fmean = CovFMean(M=g_2,optns=list(corrOut = corrOut,metric="log_cholesky"))$Mout[[1]]
  g_3_Fmean = CovFMean(M=g_3,optns=list(corrOut = corrOut,metric="log_cholesky"))$Mout[[1]]
  
  ### Generate error
  noise = mclapply(1:n, function(i){
    
    eps = matrix(rnorm(m*m,mean = 0, sd = 1/100), nrow=m)
    L = lower.part(eps) + diag(exp(diag(eps)))
    L %*% t(L)    
    
  })
  noise = simplify2array(noise)
  tf_noise = Tf_spd_list(noise, diag(m), corrOut = corrOut)
  
  tf = list()
  tf[[1]] = Tf_spd_list(g_1, g_1_Fmean, corrOut = corrOut)
  tf[[2]] = Tf_spd_list(g_2, g_2_Fmean, corrOut = corrOut)
  tf[[3]] = Tf_spd_list(g_3, g_3_Fmean, corrOut = corrOut)
  
  Y = array(diag(m),dim=c(m,m,n))
  for(j in 1:3){
    
    Y = array(sapply(1:n, function(i) {tf[[j]][[i]](Y[,,i])}),dim=c(m,m,n))
    
  }
  
  Min = array(sapply(1:n, function(i) {tf_noise[[i]](Y[,,i])}),dim=c(m,m,n))
  
  return(list(xin = X, Min = Min, Y_true = Y))
}

## Case II (ADOPT True Model, 3d Predictors, Non-linear model)
ADOPT_DatGenSPD2 = function(n = 100, m = 5, corrOut = FALSE){
  
  lower.part = function(M){
    
    lower.triangle(M) - diag(diag(M))
    
  }
  
  ## Predictors X = (X1, X2, X3)
  # Define parameters
  mu <- c(0, 0, 0)  # Mean vector
  Sigma <- matrix(c(1, 0.5, 0.3, 0.5, 1, 0.5, 0.3, 0.5, 1), nrow = 3)  # Covariance matrix
  
  # Generate bivariate normal samples size of n
  V <- mvrnorm(n, mu, Sigma)
  
  # Apply the standard normal CDF (pnorm) to each element
  X <- pnorm(V)
  
  # Generate Response
  g_1= lapply(1:n, function(i){
    
    L = matrix(0,nrow=m, ncol=m)
    for(row in 1:m){
      for(col in 1:row){
        L[row,col] = exp(-abs(row-col)/2)*sin(pi*(X[i,1])/4 + pi/4)
      }
    }
    L %*% t(L)
    
  })
  g_1 = simplify2array(g_1)
  
  g_2= lapply(1:n, function(i){
    
    L = matrix(0,nrow=m, ncol=m)
    for(row in 1:m){
      for(col in 1:row){
        L[row,col] = exp(-abs(row-col)/4)*(X[i,2]^2/2 + 1/2)
      }
    }
    L %*% t(L)
    
  })
  g_2 = simplify2array(g_2)
  
  
  g_3= lapply(1:n, function(i){
    
    L = matrix(0,nrow=m, ncol=m)
    for(row in 1:m){
      for(col in 1:row){
        L[row,col] = exp(-abs(row-col)/6)*exp(-X[i,3])
      }
    }
    L %*% t(L)
    
  })
  g_3 = simplify2array(g_3)
  
  
  g_1_Fmean = CovFMean(M=g_1,optns=list(corrOut = corrOut,metric="log_cholesky"))$Mout[[1]]
  g_2_Fmean = CovFMean(M=g_2,optns=list(corrOut = corrOut,metric="log_cholesky"))$Mout[[1]]
  g_3_Fmean = CovFMean(M=g_3,optns=list(corrOut = corrOut,metric="log_cholesky"))$Mout[[1]]
  
  ### Generate error
  noise = mclapply(1:n, function(i){
    
    eps = matrix(rnorm(m*m,mean = 0, sd = 1/100), nrow=m)
    L = lower.part(eps) + diag(exp(diag(eps)))
    L %*% t(L)    
    
  })
  noise = simplify2array(noise)
  tf_noise = Tf_spd_list(noise, diag(m),corrOut = corrOut)
  
  tf = list()
  tf[[1]] = Tf_spd_list(g_1, g_1_Fmean,corrOut = corrOut)
  tf[[2]] = Tf_spd_list(g_2, g_2_Fmean,corrOut = corrOut)
  tf[[3]] = Tf_spd_list(g_3, g_3_Fmean,corrOut = corrOut)
  
  Y = array(diag(m),dim=c(m,m,n))
  for(j in 1:3){
    
    Y = array(sapply(1:n, function(i) {tf[[j]][[i]](Y[,,i])}),dim=c(m,m,n))
    
  }
  
  Min = array(sapply(1:n, function(i) {tf_noise[[i]](Y[,,i])}),dim=c(m,m,n))
  
  return(list(xin = X, Min = Min, Y_true = Y))
}

####################################################################
######### Simulations in Section 5.2 (simulations_spd) #############
####################################################################

## Set number of cores
num_cores <- detectCores() - 1 

## Conduct simulations
run_experiment_spd = function(n, dat){
  
  ### ADOPT
  estADOPT = ADOPT_spd(xin = dat$xin, Min = dat$Min, xout = dat$xin, max_iter = 5, tol = 1e-3, corrOut = FALSE)
  MISE_cov_ADOPT = loss_SPD(estADOPT$Mfitted, dat$Y_true)/n
  
  ### Global F Reg
  fitGloCovReg =GloCovReg(x = dat$xin, M = dat$Min, xout = dat$xin, optns = list(corrOut= FALSE, metric="log_cholesky"))$Mout
  estGloCovReg = simplify2array(fitGloCovReg)
  MISE_cov_GF = loss_SPD(estGloCovReg, dat$Y_true)/n
  
  m = dim(dat$Min)[1]
  
  ### Fmean
  Fmean = CovFMean(M=dat$Min,optns=list( corrOut = FALSE, metric="log_cholesky"))$Mout[[1]]
  #Fmean = as.matrix(Matrix::nearPD(Fmean,corr = FALSE)$mat)
  #Fmean = as.matrix(Matrix::forceSymmetric(Fmean))
  
  estFmean = array(Fmean, dim=c(m,m,n))
  MISE_cov_Fmean =  loss_SPD(estFmean, dat$Y_true)/n
  
  return(c(MISE_cov_ADOPT, MISE_cov_GF, MISE_cov_Fmean))
  
}

generate_df_spd <- function(ADOPT_func, B, num_cores, nvals = c(100, 300)) {
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
      Method   = rep(c("ADOPT","GFR","FM"), each = B),
      N        = factor(n_val, levels = nvals)
    )
    out_list[[i]] <- df_temp
  }
  
  # Merge all n-values into one data frame
  df_all <- do.call(rbind, out_list)
  
  # Ensure Method is a factor in desired order
  df_all$Method <- factor(df_all$Method, levels = c("ADOPT","GFR","FM"))
  
  return(df_all)
}

## Case I with m by m SPD matrices with m = 10
ADOPT_spd1 <- function(n, b) {
  
  set.seed(b)  # Set seed for reproducibility
  dat_spd = ADOPT_DatGenSPD1(n = n, m = 10)
  
  res = run_experiment_spd(n = n, dat = dat_spd)
  
  return(res)
}

## Case I with m by m SPD matrices with m = 20
ADOPT_spd2 <- function(n, b) {
  
  set.seed(b)  # Set seed for reproducibility
  dat_spd = ADOPT_DatGenSPD1(n = n, m = 20)
  
  res = run_experiment_spd(n = n, dat = dat_spd)
  
  return(res)
}

## Case II with m by m SPD matrices with m = 10
ADOPT_spd3 <- function(n, b) {
  
  set.seed(b)  # Set seed for reproducibility
  dat_spd = ADOPT_DatGenSPD2(n = n, m = 10)
  
  res = run_experiment_spd(n = n, dat = dat_spd)
  
  return(res)
}

## Case II with m by m SPD matrices with m = 20
ADOPT_spd4 <- function(n, b) {
  
  set.seed(b)  # Set seed for reproducibility
  dat_spd = ADOPT_DatGenSPD2(n = n, m = 20)
  
  res = run_experiment_spd(n = n, dat = dat_spd)
  
  return(res)
}

## Generate Table 2
summarize_blocks <- function(x, blocks, m, scale = 100, effective_n = NULL) {
  do.call(rbind, lapply(names(blocks), function(lbl) {
    idx <- blocks[[lbl]]
    vals <- x[idx]
    n_se <- if (is.null(effective_n)) length(idx) else effective_n
    B = 200
    data.frame(
      Method      = lbl,
      n    = n_se,
      m = m,
      Mean_x100  = mean(vals) * scale,
      SE_x100    = (sd(vals) / sqrt(B)) * scale,
      check.names = FALSE,
      row.names   = NULL
    )
  }))
}



res_ADOPT_spd1 = generate_df_spd(ADOPT_spd1, B = 200, num_cores = num_cores, nvals = c(100,300))
#save(res_ADOPT_spd1, file = "res_ADOPT_spd_linear1.RData")

res_ADOPT_spd2 = generate_df_spd(ADOPT_spd2, B = 200, num_cores = num_cores, nvals = c(100,300))
#save(res_ADOPT_spd2, file = "res_ADOPT_spd_linear2.RData")

res_ADOPT_spd3 = generate_df_spd(ADOPT_spd3, B = 200, num_cores = num_cores, nvals = c(100,300))
#save(res_ADOPT_spd3, file = "res_ADOPT_spd_linear3.RData")

res_ADOPT_spd4 = generate_df_spd(ADOPT_spd4, B = 200, num_cores = num_cores, nvals = c(100, 300))
#save(res_ADOPT_spd4, file = "res_ADOPT_spd_nonlinear1.RData")


# Blocks for the two regimes
blocks_100 <- list(ADOPT = 1:200,   GF = 201:400, FM = 401:600)
blocks_300 <- list(ADOPT = 601:800, GF = 801:1000, FM = 1001:1200)

# Build the two tables (using your intended n for SE: 100 and 300)
summary_100_case1_m10 <- summarize_blocks(res_ADOPT_spd1[, 1], m = 10, blocks_100, scale = 100, effective_n = 100)
summary_300_case1_m10 <- summarize_blocks(res_ADOPT_spd1[, 1], m = 10, blocks_300, scale = 100, effective_n = 300)

summary_100_case1_m20 <- summarize_blocks(res_ADOPT_spd2[, 1], m = 20, blocks_100, scale = 100, effective_n = 100)
summary_300_case1_m20 <- summarize_blocks(res_ADOPT_spd2[, 1], m = 20, blocks_300, scale = 100, effective_n = 300)

summary_100_case2_m10 <- summarize_blocks(res_ADOPT_spd3[, 1], m = 10, blocks_100, scale = 100, effective_n = 100)
summary_300_case2_m10 <- summarize_blocks(res_ADOPT_spd3[, 1], m = 10, blocks_300, scale = 100, effective_n = 300)

summary_100_case2_m20 <- summarize_blocks(res_ADOPT_spd4[, 1], m = 20, blocks_100, scale = 100, effective_n = 100)
summary_300_case2_m20 <- summarize_blocks(res_ADOPT_spd4[, 1], m = 20, blocks_300, scale = 100, effective_n = 300)

# Combine
summary_tbl_case1_m10 <- rbind(summary_100_case1_m10, summary_300_case1_m10)
summary_tbl_case1_m20 <- rbind(summary_100_case1_m20, summary_300_case1_m20)

summary_tbl_case1 = rbind(summary_tbl_case1_m10,summary_tbl_case1_m20)
summary_tbl_case1$Case = "I"

summary_tbl_case2_m10 <- rbind(summary_100_case2_m10, summary_300_case2_m10)
summary_tbl_case2_m20 <- rbind(summary_100_case2_m20, summary_300_case2_m20)

summary_tbl_case2 = rbind(summary_tbl_case2_m10,summary_tbl_case2_m20)
summary_tbl_case2$Case = "II"

summary_tbl_spd = rbind(summary_tbl_case1,summary_tbl_case2)
summary_tbl_spd$Mean_x100 = round(summary_tbl_spd$Mean_x100, 3)
summary_tbl_spd$SE_x100 = round(summary_tbl_spd$SE_x100, 3)

print(summary_tbl_spd)


