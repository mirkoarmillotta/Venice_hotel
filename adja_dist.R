
##########################################################################################
##  Generate an adjacency matrix form independent Bernoulli trials
##  with probability p(i,j) for each pair of elements (i,j)
##  input: p matrix of probabilities, N number of elements
##  the output is the scaled adjacency matrix
##########################################################################################

adja_dist <- function(p, N){
  
  p_v <- as.matrix(vec(p))
  A_vec <- matrix(0, N*N)
  
  for( r in 1:(N*N)) {
    A_vec[r] <- rbinom(1, 1, p_v[r])
  }
  
  A <- invvec(A_vec, N, N)
  D <- diag(rowSums(A)^(-1))
  D[which(D==Inf)] <- 0
  W <- D%*%A                
  W <- as.matrix(W)
  
  return(W)
  
}
