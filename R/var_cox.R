var_cox <-
function(y, x, L, delta, sigma, beta.cox = NULL, corrected = TRUE){
  
  one<-matrix(1,nrow=length(y))
  
  diago<-function(x)
  {diag(c(diag(x)))}
  
  if (is.null(beta.cox))
    beta.cox <- beta_cox(y, x, L, delta, sigma)
  
  # BCE first-order covariance matrix
  mu <- x %*% beta.cox
  w = 1 - exp(-L ^ (1 / sigma) * exp(-mu / sigma))
  W <- diag(c(w))
  K <- (1 / sigma ^ 2) * t(x) %*% W %*% x
  
  cox <- solve(K)
  
  
  if (corrected){
    # BCE second-order covariance matrix
    w1<--(1/sigma)*L^(1/sigma)*exp(-L^(1/sigma)*exp(-mu/sigma)-mu/sigma)
    w2<--(1/sigma^2)*L^(1/sigma)*exp(-L^(1/sigma)*exp(-mu/sigma)-mu/sigma)*(L^(1/sigma)*exp(-mu/sigma)-1)
    W1<-diag(c(w1))
    W2<-diag(c(w2))
    
    tau <- c(0,-1)
    tau1 <- tau[1]
    tau2 <- tau[2]  
    W.est <- diag(c(w * (w - 2) - 2 * sigma * w1 + sigma * tau1 * (w1 + 2 * sigma * w2)))
    
    Z<-x%*%solve(K)%*%t(x)
    Zd<-diago(Z)
    Z2 <- Z * Z
    W.est.est <- diag(c(Z %*% (W + 2 * sigma * W1) %*% Zd %*% one))
    Delta.1 <- -(1 / sigma ^ 4) * t(x) %*% W.est %*% Zd %*% x
    Delta.2 <-
      -(1 / sigma ^ 6) * t(x) %*% (W %*% Z2 %*% W - 2 * sigma * W %*% Z2 %*% W1 -
                                     6 * sigma ^ 2 * W1 %*% Z2 %*% W1) %*% x
    Delta.3 <- (1 / sigma ^ 5) * t(x) %*% W1 %*% W.est.est %*% x
    Delta <- -0.5 * Delta.1 + 0.25 * Delta.2 + 0.5 * tau2 * Delta.3
    
    cox <- solve(K)+solve(K)%*%(Delta+t(Delta))%*%solve(K)
  }
  
  out <- list(cov = cox)
  
  return(out)
}
