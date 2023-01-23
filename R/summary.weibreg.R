#' Print a summary for a object of the \code{weibreg} class.
#'
#'@description Summarizes the results for a object of the \code{weibreg} class.
#'
#'@param object an object of the \code{weibreg} class.
#'@param ... 	additional arguments affecting the summary produced.
#'
#'@return A complete summary for the coefficients extracted from a \code{weibreg} object. 
#'
#'@aliases summary.weibreg
#'
#'@examples 
#'\donttest{
#'require(survival)
#'set.seed(2100)
#'
#'##Generating covariates
#'n=20 
#'x<-runif(n, max=10)
#'lambda<-exp(1.2-0.5*x)
#'sigma<-1.5
#'
#'##Drawing T from Weibull model and fixing censoring at 1.5
#'T<-rweibull(n, shape=1/sigma, scale=lambda)
#'L<-rep(1.5, n)
#'
#'##Defining the observed times and indicators of failure
#'y<-pmin(T,L)
#'delta<-ifelse(T<=L, 1, 0)
#'data=data.frame(y=y, delta=delta, x=x)
#'
#'##Fitting for Weibull regression model
#'
#'##Traditional MLE with corrected variance
#'ex1=weibfit(Surv(y,delta)~x, data=data, L=L, estimator="MLE", corrected.var=TRUE)
#'summary(ex1)
#'
#'##BCE without corrected variance
#'ex2=weibfit(Surv(y,delta)~x, data=data, L=L, estimator="BCE", corrected.var=FALSE)
#'summary(ex2)
#'
#'##BCE with corrected variance
#'ex3=weibfit(Surv(y,delta)~x, data=data, L=L, estimator="BCE", corrected.var=TRUE)
#'summary(ex3)
#'
#'##Firth's correction without corrected variance
#'ex4=weibfit(Surv(y,delta)~x, data=data, L=L, estimator="BCE", corrected.var=FALSE)
#'summary(ex4)
#'}
#'
#'@export
summary.weibreg <-
  function (object, ...) 
  {
    asterisk <- function(x) {
      if (x > 0.1) {
        ast = " "
      }
      else {
        if (x > 0.05) {
          ast = "."
        }
        else {
          if (x > 0.01) {
            ast = "*"
          }
          else {
            if (x > 0.001) {
              ast = "**"
            }
            else {
              {
                ast = "***"
              }
            }
          }
        }
      }
      return(ast)
    }
           cat("Weibull regression model\n", 
                "Method to obtain estimators: maximum likelihood estimators\n","Method to reduce bias in regression coefficients: ",ifelse(object$estimator=="MLE","none\n",
ifelse(object$estimator=="BCE","Cox and Snell's method\n","Firth's method\n")),
"Method to correct variance in finite samples: ",ifelse(object$corrected.var,"yes\n","no\n"), sep = "")
        tt <- cbind(object$coefficients, sqrt(diag(object$var)), 
        object$coefficients/sqrt(diag(object$var)), pnorm(abs(object$coefficients/sqrt(diag(object$var))), 
            lower.tail = FALSE))
    ast = sapply(tt[, 4], FUN = asterisk)
    tt = data.frame(round(tt, 5), ast)
    colnames(tt) <- c("coef", "s.e.", 
        "z value", "Pr(>|z|)", "")
    
    cat("-------------------------------------------------------------------------\n")
    cat("Regression Coefficients\n")
    print(tt[, , drop = FALSE])
    cat("---\n")
    cat("Signif. codes:  0 \"***\" 0.001 \"**\" 0.01 \"*\" 0.05 \".\" 0.1 \" \" 1\n")
    cat("-------------------------------------------------------------------------\n")
    cat("\n---\n")
  }

