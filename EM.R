library(mvtnorm)
library(foreach)
library(doParallel)
library(parallel)
library(tidyverse)

.cores = ifelse(any(str_detect(Sys.info(),"indows")),1,5)

# Object:
#   Create a EM function fit Gaussian-Mixture algorithm that take data and K cluster
# as input and produce K's mean, K's Sigma and K'p as output, as well as the total
# objective funtion(log-likelihood)
#
#    Create a predict.function that input data and give log-likelihood correspond to
# Each K cluster


objective = 
  function(data,cluster,p,mu,sigma,method = "AIC"){
    if (length(cluster) != nrow(data)) {
      return(NA)}
    
    if (method == "deviance"|method == "AIC"){
      k = length(unique(cluster))
      Pr =
        mclapply(1:nrow(data),
                 FUN = function(x){
                   #k = cluster[[x]]
                   L = lapply(1:k,
                              function(y) p[[y]]*dmvnorm(data[x,],mu[y,],sigma[[y]]))
                   L = do.call(sum,L)
                   return(L)
                 }) %>% 
        unlist()
      
      Pr[Pr<1e-10] = 1e-10
      
      L = sum(log(Pr))
      
      if (method == "AIC") {
        Df = sum(ncol(mu)*(k+1),k)-1
        return(-2*(L-Df))}
      
      return(L)}
    
    if(method == "CH"){
      k = sum(unique(cluster))
      Ce = colSums(data) %>% as.vector()
      
      Wk = mclapply(1:nrow(data), function(x){
        k = cluster[[x]]
        d = data[x,]-mu[k,]
        W = (d)%*%t(d)
        tr = sum(diag(W))
        return(tr)
      },
      mc.cores = .cores) %>% 
        do.call(sum,.)
      
      Bk = mclapply(1:nrow(data), function(x){
        k = cluster[[x]]
        d = mu[k,]-Ce
        B = d%*%t(d)
        tr = sum(diag(B))
        return(tr)
      },
      mc.cores = .cores) %>% 
        do.call(sum,.)
      
      return(Bk/Wk*((nrow(data)-k)/(k-1)))
    }
  }


# EM function

gaussian_mixture =
  function(data,
           k = 1,
           max_iter = 100,
           method = "deviance") {
    
    data = as.matrix(data) %>% scale()
    
    #test total cluster number
    if (nrow(data) < k)
      stop("Total amount of cluster cannot exceed row of data")
    
    #test if pca-ed
    #.tol = 1e-10
    #if (any((cor(data) - diag(ncol(data))) > .tol))
    #  stop("input 'data' must be a pca result")
    
    # set-up
    
    ## row & col/parameters
    
    N = nrow(data)
    
    C = ncol(data)
    
    ## latent variable/ cluster
    
    p = rep(1 / k, k)
    
    ## list of mu (with PCA, assuming all mu are important)
    
    #mu = data[sample(1:N, k),] %>% as.matrix()
    mu = kmeans(data,k)$center
    
    ## list of matrix Sigma, which should be diagonal matrix b
    
    Sigma = rerun(k, diag(C))
    
    cluster = rep(1,N)
    
    obj = objective(data,cluster,p,mu,Sigma,method)
    
    obj0 = -Inf
    
    # evaluation
    
    iter = 1
    
    while (abs(obj-obj0)>1e-2 & iter <= max_iter) {
      iter = iter + 1
      obj0 = obj
      ## E step
      Q =
        mclapply(
          X = 1:k,
          FUN =
            function(x)
              apply(data, 1, dmvnorm, mean = mu[x,], sigma = Sigma[[x]]),
          mc.cores = .cores
        ) %>% 
        do.call(cbind,.)
      
      tempmat <- matrix(rep(p,N),nrow=N,byrow = T)
      Q = (Q * tempmat) / rowSums(Q * tempmat)

      
      # M step
      ## mu_update
      mu0 = mu
      mu = t(Q) %*% data / colSums(Q)
      
      mu[is.na(mu)] = 0
      
      Sigma = 
        mclapply(1:k,
                 FUN = 
                   function(x){
                     data_ = data - matrix(rep(mu0[x,],N),nrow = N, byrow = T)
                     Sigma = matrix(0,C,C)
                     for (i in 1:N){
                       Sigma = Sigma + Q[i,x] * data_[i,] %*% t(data_[i,])
                     }
                     #Sigma = t(data_) %*% sqrt((Q[,x])%*%t(Q[,x])) %*% data_
                     Sigma = Sigma/sum(Q[,x])
                     #Sigma[Sigma<1e-6] = 0
                     Sigma[is.na(Sigma)] = 0
                     return(Sigma)
                   })
    

      p = colSums(Q)/N
      
      cluster = which(Q == apply(Q, 1, max), arr.ind = T)
      cluster = cluster[order(cluster[,1]),2]
      
      obj = objective(data,cluster,p,mu,Sigma,method)
      if (is.na(obj)) return(gaussian_mixture(data = data,k =k,
                                              max_iter = max(max_iter -iter,2)))
        
    }
    
    return(list(
      obj = obj,
      k = k,
      mu = mu,
      sigma = Sigma,
      p = p,
      cluster = cluster
    ))
  }
