
###############################################
##  Generate NAR(1) model with covariates z
###############################################

genNAR <- function(N, TT, x, z, par, W, m){
  
  lambda <- matrix(0, N, TT)
  
   G <- par[2]*W+par[3]*diag(N)
   
   d <-par[1]*rep(1, length=N)+z%*%par[4:m]
  
  for (t in 1:TT){
    lambda[,t]=d+G%*%x[,t]
  }
  
  return(lambda)
}
