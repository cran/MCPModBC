sigma.jackknife <-
function(y, x, delta){
  
  llike <- function(beta, y, x, delta)
  {
    sigma <- exp(beta[1])
    beta <- beta[-1]
    mu <- x %*% beta
    out <- -sum(delta * (-log(sigma) + (y - mu) / sigma) - 
                  exp((y - mu) / sigma))
    
    return(out)
  }
  
  
  sigma.jackk <- c()
  n <- length(y)
  p <- ncol(x)
  for (b in 1:n)
  {
    aux.fit = nlminb(rep(0, p + 1),
                     llike,
                     y = y[-b],
                     x = x[-b, ],
                     delta = delta[-b])
    sigma.jackk[b] = exp(aux.fit$par)[1]
  }
  
  sigma.p <- mean(sigma.jackk)
  aux.fit = nlminb(rep(0, p + 1),
                   llike,
                   y = y,
                   x = x,
                   delta = delta)
  out <- n * exp(aux.fit$par[1]) - (n - 1) * sigma.p
  
  return(out)
}
