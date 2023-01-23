beta_cox <-
function(y, x, L, delta, sigma, beta.mle = NULL){
  
  n <- length(y)
  
  mu <- x %*% beta.mle
  w = 1 - exp(-L ^ (1 / sigma) * exp(-mu / sigma))
  w1 <- -(1 / sigma) * L ^ (1 / sigma) * 
    exp(-L ^ (1 / sigma) * exp(-mu / sigma) - mu / sigma)
  w2 <- -(1 / sigma ^ 2) * L ^ (1 / sigma) * 
    exp(-L ^ (1 / sigma) * exp(-mu / sigma) - mu / sigma) *
    (L ^ (1 / sigma) * exp(-mu / sigma) - 1)
  W <- diag(c(w))
  W1 <- diag(c(w1))
  W2 <- diag(c(w2))
  K <- (1 / sigma ^ 2) * t(x) %*% W %*% x
  Z <- x %*% solve(K) %*% t(x)
  Zd <- diag(diag(Z))
  P <- solve(K) %*% t(x)
  one <- matrix(1, nrow = n)
  
  #estimated bias for mle of beta
  bias <- c(-(1 / (2 * sigma ^ 3)) * P %*% Zd %*% (W + 2 * sigma * W1) %*% one)
  
  out <- beta.mle - bias
  
  return(out)
}
