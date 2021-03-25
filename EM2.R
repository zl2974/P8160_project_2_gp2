library(foreach)
library(doParallel)
library(parallel)
library(tidyverse)
library(mvtnorm)

EM_MG_algrm <- function(data, ncluster){
  
  #setting
  data <- as.matrix(data) %>% scale()
  N <- nrow(data)
  q <- ncol(data)
  p_j <- rep(1/ncluster, ncluster)
  mu <-  data[sample(N, ncluster),  ] %>% as.matrix()   
  covmat <- diag(ncol(data))
  
  covList <- list()
  for(i in 1:ncluster){
    covList[[i]] <- covmat
  }
  
  count=1
  while(count <100){
    mu0 <- mu
    
    # E-step: Evaluate posterior probability, gamma
    gamma <- c()
    for(j in 1:ncluster){
      gamma2 <- apply(data,1, dmvnorm, mean = mu[j,], sigma = covList[[j]])
      gamma <- cbind(gamma, gamma2)
    }
    
    #gamma = foreach(j=1:ncluster,.combine = cbind,.inorder=T)%do%{
    #  apply(data,1, dmvnorm, mean = mu[j,], sigma = covList[[j]])
    #}

    # M- step: Calculate mu
    tempmat <- matrix(rep(p_j,N),nrow=N,byrow = T)
    r <- (gamma * tempmat) / rowSums(gamma * tempmat)  
    mu <- t(r) %*% data / colSums(r) 
    
    # M- step: Calculate Sigma and p
    for(j in 1:ncluster){
      sigma <- matrix(rep(0,q^2),ncol=q)
      for(i in 1:N){
        sigma = sigma + r[i,j] * (data[i,]-mu0[j,]) %*% t(data[i,]-mu0[j,])   
      }
      covList[[j]] <- sigma/sum(r[,j]) 
    }

    p_j <- colSums(r)/N
    count = count + 1
  }
  
  cluster <- which(r == apply(r, 1, max), arr.ind = T)
  cluster <- cluster[order(cluster[,1]),]
  return(list(mu = mu,covList = covList, p_j = p_j,cluster = cluster))
}