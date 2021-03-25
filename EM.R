library(mvtnorm)
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
           max_iter = 10) {
    data = as.matrix(data)
    
    #test total cluster number
    if (nrow(data) < k)
      stop("Total amount of cluster cannot exceed row of data")
    
    #test if pca-ed
    .tol = 1e-10
    if (any((cor(data) - diag(ncol(data))) > .tol))
      stop("input 'data' must be a pca result")
    
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
        )
      
      Q = lapply(X=1:k,FUN = function(x) p[[x]]*Q[[x]]) %>% 
        unlist() %>% 
        matrix(.,nrow = N,byrow = F)
      
      Q[Q<1e-10] = 0
      
      Q = Q/sum(Q) #N*k matrix
      
      
      # M step
      ## mu_update
      mu =
        mclapply(
          X = 1:k,
          FUN =
            function(x) {
              q = Q[,x]
              nk = sum(q)
              mu = (t(q)%*%data) / nk #1*C matrix
              return(mu)
            },
          mc.cores = .cores
        ) %>%
        unlist() %>%
        matrix(., nrow = k,byrow = T) # k*C matrix
      
      mu[mu<1e-10]=0
      
      Sigma =
        mclapply(
          X = 1:k,
          FUN =
            function(x) {
              q = Q[, x]
              nk = sum(q)
              g = data - matrix(rep(mu[x,],N),nrow = N,byrow = T)
              Sigma = t(g)%*% sqrt(q%*%t(q)) %*% (g) # C*N %*% N*1 %*% 1*N %*% N*C 
              Sigma = Sigma / nk
              Sigma[Sigma<1e-10] = 0
              return(Sigma)
            },
          mc.cores = .cores
        )
      
      p = 
        mclapply(
          X=1:k,
          FUN = function(x) sum(Q[,x]),
          mc.cores = .cores
        )
      
      
    }
    
    return(list(
      Q = Q,
      mu = mu,
      sigma = Sigma,
      p = p
    ))
  }

gaussian_mixture(sngcll_pca, 2)
