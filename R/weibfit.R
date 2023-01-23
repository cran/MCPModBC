#' Computes different estimators for the censored Weibull regression model
#' 
#' Computes the maximum likelihood estimators (MLE) for the censored Weibull
#' regression model. The bias-corrected estimators based on the Cox and Snell's
#' and Firth's methods also are available. In addition, for the covariance
#' matrix the corrected estimators discussed in Magalhaes et al. 2021 also are
#' available.
#' 
#' The Weibull distribution considered here has probability density function
#' \deqn{ f(t;\lambda, \sigma)=\frac{1}{\sigma
#' \lambda^{1/\sigma}}t^{1/\sigma-1}\exp\left\{-\left(\frac{t}{\lambda}\right)^{1/\sigma}\right\},
#' \quad t, \sigma, \lambda>0. } The regression structure is incorporated as
#' \deqn{ \log(\lambda_i)={\bm x}_i^\top {\bm \beta}, \quad i=1,\ldots,n. } For
#' the computation of the bias-corrected estimators, \eqn{\sigma} is assumed as
#' fixed in the jackknife estimator based on the traditional MLE.
#' 
#' The Fisher information matrix for \eqn{\bm \beta} is given by \eqn{{\bm
#' K}=\sigma^{-2} {\bm X}^\top {\bm W} {\bm X}}, where \eqn{{\bm X}=({\bm
#' x}_1,\ldots,{\bm x}_n)^\top}, \eqn{{\bm W}=\mbox{diag}(w_1,\ldots,w_n)}, and
#' \deqn{ w_i=E\left[\exp\left(\frac{y_i-\log
#' \lambda_i}{\sigma}\right)\right]=q \times \left\{ 1 - \exp\left[
#' -L_i^{1/\sigma} \exp(-\mu_i/\sigma) \right] \right\} + (1-q)\times
#' \left(r/n\right), } with \eqn{q = P\left(W_{(r)}\leq \log L_i\right)} and
#' \eqn{W_{(r)}} denoting the \eqn{r}th order statistic from \eqn{W_1, \ldots,
#' W_n}, with \eqn{q=1} and \eqn{q=0} for types I and II censoring,
#' respectively. (See Magalhaes et al. 2019 for details).
#' 
#' The bias-corrected maximum likelihood estimator based on the Cox and Snell's
#' method (say \eqn{\widetilde{\bm \beta}}) is based on a corrective approach
#' given by \eqn{\widetilde{\beta}=\widehat{\beta}-B(\widehat{\beta})}, where
#' \deqn{ B({\bm \beta})= - \frac{1}{2 \sigma^3} {\bm P} {\bm Z}_d \left({\bm
#' W} + 2 \sigma {\bm W}^{\prime}\right) {\bm 1}, } with \eqn{{\bm P} = {\bm
#' K}^{-1} {\bm X}^{\top}}, \eqn{{\bm Z} = {\bm X} {\bm K}^{-1} {\bm
#' X}^{\top}}, \eqn{{\bm Z}_d} is a diagonal matrix with diagonal given by the
#' diagonal of \eqn{{\bm Z}}, \eqn{{\bm W}^{\prime} =} diag\eqn{(w_1^{\prime},
#' \ldots, w_n^{\prime})}, \eqn{w_i^{\prime} = - \sigma^{-1} L_i^{1/\sigma}
#' \exp\{ -L_i^{1/\sigma} \exp(-\mu_i/\sigma) - \mu_i/\sigma \}} and \eqn{{\bm
#' 1}} is a \eqn{n}-dimensional vector of ones.
#' 
#' The bias-corrected maximum likelihood estimator based on the Firth's method
#' (say \eqn{\check{\bm \beta}}) is based on a preventive approach, which is
#' the solution for the equation \eqn{{\bm U}_{{\bm \beta}}^{\star} = {\bm 0}},
#' where \deqn{ {\bm U}_{{\bm \beta}}^{\star} = {\bm U}_{{\bm \beta}} - {\bm
#' K}_{{\bm \beta} {\bm \beta}} B({\bm \beta}). }
#' 
#' The covariance correction is based on the general result of Magalhaes et al.
#' 2021 given by \deqn{ \mbox{{\bf Cov}}_{\bm 2}^{\bm \tau}({\bm
#' \beta}^{\star}) = {\bm K}^{-1} + {\bm K}^{-1} \left\{ {\bm \Delta} + {\bm
#' \Delta}^{\top} \right\} {\bm K}^{-1} + \mathcal{O}(n^{-3}) } where
#' \eqn{{\Delta} = -0.5 {\Delta}^{(1)} + 0.25 {\Delta}^{(2)} + 0.5 \tau_2
#' {\Delta}^{(3)}}, with \deqn{ \Delta^{(1)} = \frac{1}{\sigma^4} {\bm
#' X}^{\top} {\bm W}^{\star} {\bm Z}_{d} {\bm X}, } \deqn{ \Delta^{(2)} = -
#' \frac{1}{\sigma^6} {\bm X}^{\top} \left[ {\bm W} {\bm Z}^{(2)} {\bm W} - 2
#' \sigma {\bm W} {\bm Z}^{(2)} {\bm W}^{\prime} - 6 \sigma^2 {\bm W}^{\prime}
#' {\bm Z}^{(2)} {\bm W}^{\prime} \right] {\bm X}, } and \deqn{ \Delta^{(3)} =
#' \frac{1}{\sigma^5} {\bm X}^{\top} {\bm W}^{\prime} {\bm W}^{\star\star} {\bm
#' X}, } where \eqn{{\bm W}^{\star} =} diag\eqn{(w_1^{\star}, \ldots,
#' w_n^{\star})}, \eqn{w_i^{\star} = w_i (w_i -2) - 2 \sigma w_i^{\prime} +
#' \sigma \tau_1 (w_i^{\prime} + 2 \sigma w_i^{\prime\prime})}, \eqn{{\bm
#' Z}^{(2)} = {\bm Z} \odot {\bm Z}}, with \eqn{\odot} representing a direct
#' product of matrices (Hadamard product), \eqn{{\bm W}^{\star\star}} is a
#' diagonal matrix, with \eqn{{\bm Z} ( {\bm W} + 2 \sigma {\bm W}^{\prime})
#' {\bm Z}_{d} {\bm 1}} as its diagonal, \eqn{{\bm W}^{\prime\prime} =}
#' diag\eqn{(w_1^{\prime\prime}, \ldots, w_n^{\prime\prime})},
#' \eqn{w_i^{\prime\prime} = - \sigma^{-1} w_i^{\prime} \left[ L_i^{1/\sigma}
#' \exp(-\mu_i/\sigma) - 1 \right]}, \eqn{{\bm \tau} = (\tau_1, \tau_2) = (1,
#' 1)} indicating the second-order covariance matrix of the MLE \eqn{{\bm
#' \beta}^{\star} = \widehat{{\bm \beta}}} denoted by \eqn{\mbox{{\bf
#' Cov}}_{\bm 2}(\widehat{{\bm \beta}})} and \eqn{{\bm \tau} = (0, -1)}
#' indicating the second-order covariance matrix of the BCE \eqn{{\bm
#' \beta}^{\star} = \widetilde{{\bm \beta}}} denoted by \eqn{\mbox{{\bf
#' Cov}}_{\bm 2}(\widetilde{{\bm \beta}})}.
#' 
#' @param formula A formula that contains on the left hand side an object of
#' the type Surv and on the right hand side the covariates definition
#' @param data A data.frame in which the formula argument can be evaluated
#' @param L the prefixed censoring times. \eqn{L=\infty} by default.
#' @param estimator the class of estimator used: MLE (maximum likelihood
#' estimator, by default), BCE (bias-corrected estimator based on the Cox and
#' Snell's method) and Firth (bias-corrected estimator based on the Firth's
#' method).
#' @param corrected.var should the covariance-corrected estimator be used?
#' (FALSE by default). See details.
#' @return \item{coefficients}{a vector with the estimated coefficients for
#' \eqn{\bm \beta}.} \item{var}{a matrix with the estimated covariance matrix
#' for the estimates of the regression coefficients \eqn{\bm \beta}}
#' \item{scale}{the estimated scale parameter \eqn{\sigma}} \item{loglik}{the
#' value for the logarithm of the likelihood function evaluated in the
#' estimates of \eqn{\bm \beta} and \eqn{\sigma}} \item{linear.predictors}{a
#' vector with the estimated linear predictor \eqn{{\bm x}_i^\top {\bm \beta}}}
#' \item{y}{a vector with the observed times (possibly censored)}
#' \item{estimator}{the estimator used for \eqn{\bm \beta}: MLE, BCE or Firth.}
#' \item{corrected.var}{logical. TRUE if a correction for the covariance was
#' used, FALSE otherwise.}
#' 
#' @author Gallardo D.I., Diniz, M.A., Magalhaes, T.M.
#' 
#' @references Cox, D.R., Snell E.J. A general definition of residuals Journal
#' of the Royal Statistical Society. Series B (Methodological).
#' 1968;30:248-275.
#' 
#' Firth, D. Bias reduction of maximum likelihood estimates Biometrika.
#' 1993;80:27-38.
#' 
#' Magalhaes Tiago M., Botter Denise A., Sandoval Monica C. A general
#' expression for second- order covariance matrices - an application to
#' dispersion models Brazilian Journal of Probability and Statistics.
#' 2021;35:37-49.
#' 
#' @examples
#' 
#' \donttest{
#' require(survival)
#' set.seed(2100)
#' 
#' ##Generating covariates
#' n=20; 
#' x<-runif(n, max=10)
#' lambda<-exp(1.2-0.5*x); sigma<-1.5
#' 
#' ##Drawing T from Weibull model and fixing censoring at 1.5
#' T<-rweibull(n, shape=1/sigma, scale=lambda); L<-rep(1.5, n)
#' 
#' ##Defining the observed times and indicators of failure
#' y<-pmin(T,L); 
#' delta<-ifelse(T<=L, 1, 0)
#' data=data.frame(y=y, delta=delta, x=x)
#' 
#' ##Fitting for Weibull regression model
#' 
#' ##Traditional MLE with corrected variance
#' ex1=weibfit(Surv(y,delta)~x, data=data, L=L, estimator="MLE", corrected.var=TRUE)
#' summary(ex1)
#' 
#' ##BCE without corrected variance
#' ex2=weibfit(Surv(y,delta)~x, data=data, L=L, estimator="BCE", corrected.var=FALSE)
#' summary(ex2)
#' 
#' ##BCE with corrected variance
#' ex3=weibfit(Surv(y,delta)~x, data=data, L=L, estimator="BCE", corrected.var=TRUE)
#' summary(ex3)
#' 
#' ##Firth's correction without corrected variance
#' ex4=weibfit(Surv(y,delta)~x, data=data, L=L, estimator="BCE", corrected.var=FALSE)
#' summary(ex4)
#' }
#' 
#' @importFrom nleqslv nleqslv
#' @importFrom stats model.frame nlminb pnorm
#' 
#' @export
weibfit <-
  function(formula, data, L=Inf, estimator="MLE", corrected.var=FALSE)
  {
    if (!inherits(formula, "formula")) {
      if (inherits(formula, "data.frame")) 
        warning("You gave a data.frame instead of a formula.")
      stop("formula is not an object of type formula")
    }
    if (!inherits(data, "data.frame")) {
      if (inherits(data, "formula")) 
        warning("You gave a formula instead of a data.frame.")
      stop("data is not an object of type data.frame")
    }
    if (missing(formula)) 
      stop("Missing formula")
    if (missing(data)) 
      stop("Missing data")
    mf <- model.frame(formula, data)
    Y <- mf[[1]]
    if (!inherits(Y, "Surv")) 
      stop("left hand side not a survival object")
    x <- model.matrix(formula, data)
    t <- Y[, 1]
    delta <- Y[, 2]
    #if (is.null(prec)) 
    #    stop("prec must be specified")
    #if (is.null(max.iter)) 
    #    stop("max.iter must be specified")
    y <- c(log(t))
    delta <- c(delta)
    cluster <- c(cluster)
    #max.iter <- round(max.iter)
    if (!any(estimator == c("MLE", "BCE", "Firth"))) 
      stop("estimator argument is not recognized")
    if (length(t) != length(delta)) 
      stop("t and delta don't have the same length")
    #if (prec > 1) 
    #    stop("prec is too high")
    #if (max.iter <= 0) 
    #    stop("max.iter at least 1")
    

    llike <- function(beta, y, x, delta)
    {
      sigma <- exp(beta[1])
      beta <- beta[-1]
      mu <- x %*% beta
      out <- -sum(delta * (-log(sigma) + (y - mu) / sigma) - 
                    exp((y - mu) / sigma))
      return(out)
    }

    surv.mod=survreg(Surv(t, delta) ~ x[,-1], dist='weibull')
    
    ### MLE
    beta <- coef(surv.mod)
    names(beta)<-colnames(x)
    sigma<-surv.mod$scale
    if(estimator=="BCE")
    {
      sigma <- sigma.jackknife(y, x, delta)
      beta <- beta_cox (y, x, L, delta, sigma=sigma, beta.mle = beta)
    }
    if(estimator=="Firth")
    {
      sigma <- sigma.jackknife(y, x, delta)
      beta <- beta_firth (y, x, L, delta, sigma=sigma, beta.mle = beta)
    }
    vari=surv.mod$var[1:ncol(x),1:ncol(x)]
    if(estimator=="MLE" & corrected.var==TRUE)
    {vari=var_mle2 (y, x, L, delta, sigma, beta.mle = beta)$cov}
    if(estimator=="BCE")
    {vari=var_cox (y, x, L, delta, sigma, beta.cox = beta,corrected=corrected.var)$cov}
    if(estimator=="Firth")
    {vari=var_firth (y, x, L, delta, sigma, beta.firth = beta)$cov}
    colnames(vari)=rownames(vari)=colnames(x)
    
    
    out<-list(coefficients=beta)
    out$var<-vari
    out$scale<-sigma
    out$loglik<-llike(c(log(sigma),beta), y, x, delta)
    out$linear.predictors<-c(x%*%beta)
    out$y<-c(y)
    out$estimator<-estimator
    out$corrected.var<-corrected.var
    
    class(out)<-"weibreg"
    out
  }
