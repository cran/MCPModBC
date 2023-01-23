beta_firth <-
function(y, x, L, delta, sigma, beta.mle = NULL) {
  
  n <- length(y)
  
  mod.score <- function(x, y, X, L, sigma, delta)
  {
    diago <- function(x)
    {
      diag(c(diag(x)))
    }
    
    beta <- x
    x <- X
    mu <- x %*% beta
    w = 1 - exp(-L ^ (1 / sigma) * exp(-mu / sigma))
    w1 <-
      -(1 / sigma) * L ^ (1 / sigma) * exp(-L ^ (1 / sigma) * exp(-mu / sigma) -
                                             mu / sigma)
    W <- diag(c(w))
    W1 <- diag(c(w1))
    K <- (1 / sigma ^ 2) * t(x) %*% W %*% x
    Z <- x %*% solve(K) %*% t(x)
    Zd <- diago(Z)
    P <- solve(K) %*% t(x)
    one <- matrix(1, nrow = n)
    bias <- c(-(1 / (2 * sigma ^ 3)) * P %*% Zd %*% (W + 2 * sigma * W1) %*%
                one)
    v <- (-delta + exp((y - mu) / sigma)) * (1 / sqrt(w))
    U <- c((1 / sigma) * t(x) %*% sqrt(W) %*% v)
    c(U - K %*% bias)
  }
  
  out <-
    try(nleqslv(
      x = beta.mle,
      fn = mod.score,
      y = y,
      X = x,
      L = L,
      sigma = sigma,
      delta = delta
    )$x,
    silent = TRUE)
  return(out)
}
