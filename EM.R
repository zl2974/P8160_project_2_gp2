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


# EM function

gaussian_mixture =
  function(data,
           k = 1,
           max_iter = 1e+5) {
    
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
    
    mu = data[sample(1:N, k),] %>% as.matrix()
    
    
    ## list of matrix Sigma, which should be diagonal matrix b
    
    Sigma = rerun(k, diag(C))
    
    # evaluation
    
    iter = 1
    
    while (iter <= max_iter) {
      iter = iter + 1
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
                     return(Sigma)
                   })

      p = colSums(Q)/N
      
      
    }
    
    return(list(
      Q = Q,
      mu = mu,
      sigma = Sigma,
      p = p
    ))
  }
