
###################################################################################
##  Function for the OLS estimation of linear NAR(1) model with covariates
###################################################################################

ols.nar1 <- function(N, TT, y, x, z, W, m){

  z <- as.matrix(z)
  
  XX <- matrix(0, m, m)
  Xy <- matrix(0, m, 1)
  
  x <- as.matrix(x)
  y <- as.matrix(y)
    
    for (t in 1:TT){
      XXt <- cbind(rep(1,N), W%*%x[,t], x[,t], z)
      XX <- XX + t(XXt)%*%XXt
      Xy <- Xy + t(XXt)%*%y[,t]
    }
  
  theta <- solve(XX)%*%Xy
  return(list(beta=theta, XX=XX))
  
}
