var_firth <-
function(y, x, L, delta, sigma, beta.firth = NULL){
  
  if (is.null(beta.firth))
    beta.firth <- beta_firth(y, x, L, delta, sigma)
  
  # MLE
  mu <- x %*% beta.firth
  w = 1 - exp(-L ^ (1 / sigma) * exp(-mu / sigma))
  W <- diag(c(w))
  K <- (1 / sigma ^ 2) * t(x) %*% W %*% x
  
  firth <- solve(K)
  
  out <- list(cov = firth)
  
  return(out)
}
