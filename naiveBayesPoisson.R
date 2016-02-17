library(matrixStats)
naiveBayesPoisson<-function(x,z){
  delta=1e-3
  z = model.matrix(~0+z)
  #browser()
  alpha = t(x)%*%z+delta
  beta = colSums(z)+delta
  return(list(alpha,beta))
}

naiveBayesPoissonPredict<-function(x,z,xstar){
  delta=1e-3
  z = model.matrix(~0+z)
  #browser()
  alpha = t(z)%*%as.matrix(x)+delta
  beta = colSums(z)+delta
  
  K = length(beta)
  lnp = matrix(nrow = dim(xstar)[1], ncol = K)
  for (k in 1:K){
    xnd_alpha = sweep(xstar,2,alpha[k,],'+')
    t1 = alpha[k,]*log(beta[k])-lgamma(alpha[k,])
    t2 = lgamma(xnd_alpha)-lgamma(xstar+1)-log(beta[k]+1)*xnd_alpha
    lnp[,k] = rowSums(sweep(t2,2,t1,'+'))
    #browser()
  }
  
  lnZ = apply(lnp, 1,logSumExp)
  p = exp(sweep(lnp,1,lnZ,'-'))
}