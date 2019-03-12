## Early-Burst functions

# covariance matrix
covmatEB <- function(tree, rateMat, g){
  
  phylomat <- vcv(tree)
  transphylomat <- (exp(g*phylomat) - 1)/g
  res <- kronecker(X = transphylomat, Y = rateMat, make.dimnames = T)
  
  return(res)
}

# likelihood for EB (withouth including weight matrix because we are not using natural selection here)
computeMultiLk <- function(Data, meanvec, covMat){
  
  n <- dim(Data)[1]
  m <- dim(Data)[2]
  
  meanvec <- as.vector(t(meanvec))
  Data <- as.vector(t(Data))
  denomchunk <- sqrt(det(covMat)*(2*pi)^(n*m))
  likelihood = exp(-0.5*t(Data - meanvec) %*% solve(covMat) %*% (Data - meanvec))/denomchunk
  
  likelihood
}

