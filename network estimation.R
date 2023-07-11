

###################################################################################
## Simulate networks for specified pairwise probability values p_ij
## to draw an edge from hotel i to j (i-> j)
## Select the network which minimize the Root Mean Squared Error;
## this is done for each advance booking in k=1,...,K
## Input: - N, number of hotels - TT, temporal size, - K total number of advance
## booking - S, total number of simulations - inv_dist, matrix of probabilities
## p_ij to draw an edge from hotel i to j (i-> j)
## - ho, list of matrices of advance booking formatted using function 'data_format'
## - z covariates (if no covariates put z=FALSE) 
## - dum dummy variables for weekend and holidays
###################################################################################


NetEstgam <- function(N, TT, K, S, ho, z, p_dist){
  
  net <- list()
  par <- list()
  lam <- list()
  res <- list()
  std <- list()
  nete <- list()
  pare <- list()
  lame <- list()
  stde <- list()
  net.se <- list()
  netA.se <- list()
  par.se <- list()
  lam.se <- list()
  std.se <- list()
  ss<- vector()
  map <- vector()
  sse <- vector()
  mse <- vector()
  p.se <- list()
  a <- vector()
  tlm <- vector()
  lm <- list()
  netej <- list()
  parej <- list()
  lamej <- list()
  stdej <- list()
  beta1 <- vector()
  se1 <- vector()
  
  m <- 3+ncol(z)
    
    for(k in 1:K){ 
      
      for (i in 1:S){ 
        
        net[[i]] <- adja_dist(p=p_dist, N=N)
        res[[i]] <- ols.nar1(N=N, TT=TT, y=ho[[k]], x=ho[[k+1]], z=z, W=net[[i]], m=m)
        par[[i]] <- res[[i]]$beta
        lam[[i]] <- genNAR(N=N, TT=TT, ho[[k+1]], z, par[[i]], net[[i]], m=m)
        ss[i] <- mean((ho[[k]]-lam[[i]])^2)                    ## variance
        std[[i]] <- diag((sqrt(ss[i]))^2*solve(res[[i]]$XX))   ## variance betas
        map[i] <- sqrt(ss[i])                                  ## RMSE
      }
      
      map <- as.matrix(map)
      best <- which(map==min(map), arr.ind=T)[1]
      sse[k] <- min(map)
      
      nete[[k]] <- net[[best]] 
      pare[[k]] <- par[[best]]
      lame[[k]] <- lam[[best]]
      stde[[k]] <- std[[best]]
      netA.se[[k]] <- ifelse(nete[[k]]==0,0,1)
    }
  
      results <- matrix(0, K, 22)
      colnames(results) <- c("Intercept",  "Netowrk", "Lag", "Restaurant", "Meeting", "Size", "Geo",
                             "std_Int",  "std_Net", "std_Lag", "std_Rest", "std_Meet", "std_Size", "std_Geo",
                             "t_Int", "t_Net", "t_Lag", "t_Rest", "t_Meet", "t_Size", "t_Geo",
                             "Min. RMSE")
  
    
  for (k in 1:K){
    results[k,] <- c(pare[[k]], sqrt(stde[[k]]), pare[[k]]/sqrt(stde[[k]]), sse[k])
  }
  
  rownames(results) <- seq(from=0, to=K-1, by=1)
  return(list(results=results, hat_W=nete, hat_A=netA.se, hat_lambda=lame))
  
}
