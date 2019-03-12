library(phangorn)
library(ape)
library(plyr)
library(phytools)

####### Pau's Hansen model functions prior to Java implementation ##########

####### UNIVARIATE CASE

# Description: The function computes the covariance matrix for an OU unitrait model in the case
# in which the root is supposed to be fixed and in the case it is supposed to be random.

# Remark: The output of varOU should be multiplied by sigma^2

varOU <- function(tree, alpha, fixedRoot = T){
  
  if(alpha <= 0) stop("alpha should be positive for an OU matrix") #Because the model requires it to be positive
  bmvar <- vcv.phylo(tree)
  n <- length(tree$tip.label)
  ouvar <- matrix(data = NA, nrow = n, ncol = n)
  
  if (fixedRoot){
    for (row in 1:n){
      for (col in 1:n){
        ouvar[row,col] = exp(-alpha*(bmvar[row,row] + bmvar[col,col] - 2*bmvar[row,col]))*(1 - exp(-2*alpha*bmvar[row,col]))
      }
    }
  } else{
      for (row in 1:n){
        for (col in 1:n){
          ouvar[row,col] = exp(-alpha*(bmvar[row,row] + bmvar[col,col] - 2*bmvar[row,col]))
        }
      }
  }
  ouvar = ouvar/(2*alpha)

  return(ouvar)
}


# Description: The function outputs the weight matrix for an OU unitrait model in the case the root
#   is merged with the eldest selective regime and in the case one want to estimate it isolated.
weightMat <- function(tree, alpha, n.regimes, branch.regimes, mergeRoot = F){
  
  n <- length(tree$tip.label)
  WMat <- matrix(data = 0, nrow = n , ncol = n.regimes)
  regimes <- mapvalues(branch.regimes, from = levels(branch.regimes), to = 1:length(levels(branch.regimes)))
 
  for(sp in 1:n){
    ancestors <- c(sp, Ancestors(x = tree, node = sp)) # require library (phangorn)
   
     for(node in 1:(length(ancestors) - 1)){
      WMat[sp,regimes[ancestors[node]]] = WMat[sp,regimes[ancestors[node]]] + exp(alpha*nodeheight(tree, ancestors[node])) - exp(alpha*nodeheight(tree, ancestors[node + 1]))
    }
  }
  
  if (mergeRoot){
    rootlabel <- length(tree$tip.label) + 1 # To identify the label of the oldest selective regime
    WMat[, regimes[rootlabel]] = WMat[, regimes[rootlabel]] + 1
    
  } else {
    WMat = cbind(1, WMat)
  }
  return(exp(-alpha*diag(vcv(tree)))*WMat)
}

# Description: Likelihood function of the multivariate normal distribution
computeLk <- function(Data, param, Sigma, sigmasq){
  n <- dim(Sigma)[1]
  likelihood = exp(-0.5*t(Data - param)%*%solve(Sigma)%*%(Data - param)/sigmasq)/((2*pi*sigmasq)^(n/2)*det(Sigma)^(1/2))
  likelihood
}

####### MULTIVARIATE CASE

# Description: Computes exp(-A*t) for a given time t
expA <- function(eigen.values, eigen.Matvector, t){
  
  eigen.Matvalues.exp <- diag(exp(- eigen.values * t))
  
  res <- eigen.Matvector %*% eigen.Matvalues.exp %*% solve(eigen.Matvector)
  
  return(res)
}

# Description: Auxiliar function for computing the blocks of the multi trait covariance matrix
computEigenChunk <- function(eigenValues, cii, cjj, cij, fixedRoot = T){
  m <- length(eigenValues)
  res <- matrix(data = 0, nrow = m, ncol = m)
  
  if(fixedRoot){
    for(k in 1:m){
      for(l in 1:m){
        res[k,l] = (1/(eigenValues[k] + eigenValues[l]))*(1 - exp(-(eigenValues[k] + eigenValues[l])*cij))
      }
    }
  }
  
  else{
    for(k in 1:m){
      for(l in 1:m){
        res[k,l] = (1/(eigenValues[k] + eigenValues[l]))*exp(-eigenValues[k]*(cii - cij))*exp(-eigenValues[l]*(cjj - cij))
      }
    }
  }
  
  res
}

# Description: Multivariate covariance matrix for OU model (from mvMORPH)
# It computes the eigen descomposition of Alpha and the Cholesky decomposition of R
MultivarOU <- function(tree, Alpha, Drift, fixedRoot = T){
  
  # Veryfing dimension matching between Alpha and Drift
  if(dim(Alpha)[1] != dim(Alpha)[2]) stop("Alpha must be a square matrix")
  if(dim(Drift)[1] != dim(Drift)[2]) stop("Drift must be a square matrix")
  if(dim(Drift)[1] != dim(Alpha)[1]) stop("Alpha and Drift must have same dimensions")
  m <- dim(Drift)[1] # Number of traits
  
  # A Eigen value decomposition
    eigenA <- eigen(Alpha)
    P <- eigenA$vectors # Using the same notation as the mvMORPH article
    eigenVal <- eigenA$values
    if(any(is.complex(eigenVal))){ stop(paste("A matrix has complex number eigen values. List of eigenvalues: \n", eigenVal)) }
    Pinv <- solve(P)
      # Because the model requires it to have real eigen values
  
  # R matrix (drift) cholesky decomp
    #drift.chol <- t(chol(Drift)) # Take the lower triangular as default to match it with the formulas
    
  # Phylogenetic matrix
    bmvar <- vcv.phylo(tree)
    n <- length(tree$tip.label)
  
  ouvar <- matrix(data = NA, nrow = n*m, ncol = n*m)
  
  commonChunk <- Pinv %*% Drift %*% t(Pinv)
  
  if (fixedRoot){
    for (row in 1:n){
      for (col in 1:n){
        row.elements <- seq(from = ((row - 1)*m + 1), to = row*m)
        col.elements <- seq(from = ((col - 1)*m + 1), to = col*m)
        bigChunk <- computEigenChunk(eigenVal, bmvar[row,row], bmvar[col,col], bmvar[row,col], T)*commonChunk
        bigChunk <- P %*% bigChunk %*% t(P)
        ouvar[row.elements, col.elements] = expA(eigenVal, P, bmvar[row,row] - bmvar[row, col]) %*% bigChunk %*% expA(eigenVal, t(Pinv), bmvar[col,col] - bmvar[row, col])
        
        print(ouvar[row.elements, col.elements])
      }
    }
  } else{
    for (row in 1:n){
      for (col in 1:n){
        row.elements <- seq(from = ((row - 1)*m + 1), to = row*m)
        col.elements <- seq(from = ((col - 1)*m + 1), to = col*m)
        bigChunk <- computEigenChunk(eigenVal, bmvar[row,row], bmvar[col,col], bmvar[row,col], F)*commonChunk
        ouvar[row.elements, col.elements] = P %*% bigChunk %*% t(P)
      }
    }
  }
  
  return(ouvar)
}

# Description: Weight matrix of the multi Trait case allowing merging the root with
# the primary optimum value (the optimum value associated with the eldest selective regime)
# or not. It uses functions expA to make the calculations
multiweightMat <- function(tree, Alpha, branch.regimes, mergeRoot = F){
  
  if(dim(Alpha)[1] != dim(Alpha)[2]) stop("Alpha must be a square matrix")
  m <- dim(Alpha)[1] # Number of traits
  
  # A Eigen value decomposition
  eigenA <- eigen(Alpha)
  P <- eigenA$vectors # Using the same notation as the mvMORPH article
  eigenVal <- eigenA$values
  if(any(is.complex(eigenVal))){ stop(paste("A matrix has complex number eigen values. List of eigenvalues: \n", eigenVal)) }
  Pinv <- solve(P)
    # Because the model requires it to have real eigen values
  
  # Phylogenetic matrix
  bmvar <- vcv.phylo(tree)
  n <- length(tree$tip.label)
    
  regimes <- mapvalues(branch.regimes, from = levels(branch.regimes), to = 1:length(levels(branch.regimes)))
  r <- length(levels(regimes))
  WMat <- matrix(data =0, nrow = n * m, ncol = r * m)
  
  for (sp in 1:n) {
    print(paste("sp = ", sp))
    ancestors <- c(sp, Ancestors(x = tree, node = sp))
    print(paste("Ancestors of species ", sp))
    print(ancestors)
    for (k in 1:r) {
      print(paste("k = ", k))
      row.elements <- seq(from = m * (sp - 1) + 1, to = m * sp, by = 1) # The select the rows of the current species 
      col.elements <- seq(from = m * (k - 1) + 1, to = m * k, by = 1) # To select the columns of the selective regime
      
      W_spk <- matrix(data =0, nrow = m, ncol = m)
      
      for (node in 1:(length(ancestors) - 1)) {
        print(paste("ancestor = ", ancestors[node]))
        print(paste("regime of node =", regimes[ancestors[node]]))
        if (regimes[ancestors[node]] == k){
          W_spk = W_spk + expA(eigenVal, P, nodeheight(tree, sp) - nodeheight(tree, ancestors[node])) - expA(eigenVal, P, nodeheight(tree, sp) - nodeheight(tree, ancestors[node + 1]))
        }
        print("W_spk = ")
        print(W_spk)
      }
      WMat[row.elements, col.elements] = W_spk
    }
  }
  
  if (mergeRoot) {
    rootlabel <- length(tree$tip.label) + 1 # To identify the label of the eldest selective regime
    eldest.regim.labels <- seq(from = (as.numeric(regimes[rootlabel]) - 1)*m + 1, to = as.numeric(regimes[rootlabel])*m, by = 1)
    
    for (sp in 1:n){
      row.elements <- seq(from = m*(sp - 1) + 1, to = m*sp, by = 1) # The select the rows of the current species 
      WMat[row.elements, eldest.regim.labels] = WMat[row.elements, eldest.regim.labels] + expA(eigenVal, P, nodeheight(tree, sp))
    }
  } else{
      W0 = matrix(data =0, nrow = n*m, ncol = m)
      
      for (sp in 1:n){
        selrow <- seq(from = (sp - 1)*m + 1, to = sp*m, by = 1)
        
        W0[selrow,] = expA(eigenVal, P, nodeheight(tree, sp))
      }
      print("W0 = ")
      print(W0)
      
      WMat = cbind(W0, WMat)
    }
  
  return(WMat)
}


computeMultiLk <- function(Data, param, weightMat, covMat){
  
  n <- dim(Data)[1]
  m <- dim(Data)[2]
  
  param <- as.vector(t(param))
  Data <- as.vector(t(Data))
  print(Data)
  print(param)
  denomchunk <- sqrt(det(covMat)*(2*pi)^(n*m))
  print(denomchunk)
  likelihood = exp(-0.5*t(Data - weightMat %*% param) %*% solve(covMat) %*% (Data - weightMat %*% param))/denomchunk
  
  likelihood
}