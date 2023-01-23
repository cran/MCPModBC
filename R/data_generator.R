#' Data Generator
#' 
#' @description It generates data for a dose-finding trial.
#'
#' @param doses a numeric vector indicating the doses that will be considered in the clinical trial.
#' @param sample.size a numeric value indicating the sample size per dose in the clinical trial.
#' @param distr a character value indicating the distribution of the response variable. Currently, the only option available is 'weibull'.
#' @param parm a named list of true values for the simulation. See mode in Details.
#' @param censoring.rate a numeric value between 0 and 1 indicating the censoring rate when generated data. It is required when \code{distr = "weibull"}.
#'
#'@details If \code{distr = "weibull"}, the list \code{parm} should contain two components - \code{lambda} and \code{sigma} - that are 
#'the scale and shape parameters in the following parametrization of the Weibull distribution: 
#'
#' \deqn{
#'f(t;\lambda,\sigma) = \frac{1}{\sigma\lambda^{1/\sigma}} t^{1/\sigma-1} \exp\left\{ -\left( t/\lambda \right)^{1/\sigma} \right\}, \ t > 0,
#'}
#'
#'with hazard rate given by
#'
#'\deqn{
#'h(t) = \frac{1}{\lambda^{1/\sigma}\sigma}t^{1/\sigma - 1}
#'}
#'
#'and regression structure
#'\deqn{
#'log(\lambda_i) = d_i\beta_i.
#'}
#'where \eqn{log(\lambda_i)} represents the model effect for dose i, \code{doses[i]}.
#'
#' @return a data frame of dimension length(doses)*sample.size \eqn{\times} 3 when \code{distr = "weibull"} 
#' containing time-to-event, censoring indicator and dose;
#'
#' @examples
#' \donttest{ 
#' library(DoseFinding)
#' library(MCPModBC)
#' 
#' ## doses scenarios 
#' doses <- c(0, 5, 25, 50, 100)
#' nd <- length(doses)
#' 
#' # median survival time for placebo dose
#' mst.control <- 4 
#' 
#' # shape parameter
#' sigma.true <- 0.5
#'
#' # maximum hazard ratio between active dose and placebo dose 
#' hr.ratio <- 4  
#' # minimum hazard ratio between active dose and placebo dose
#' hr.Delta <- 2 
#' 
#' # hazard rate for placebo dose
#' placEff <- log(mst.control/(log(2)^sigma.true)) 
#' 
#' # maximum hazard rate for active dose
#' maxEff <- log((mst.control*(hr.ratio^sigma.true))/(log(2)^sigma.true))
#' 
#' # minimum hazard rate for active dose
#' minEff.Delta <- log((mst.control*(hr.Delta^sigma.true))/(log(2)^sigma.true))
#' Delta <- (minEff.Delta - placEff)
#' 
#' ## MCP Parameters  
#' emax <- guesst(d = doses[4], p = 0.5, model="emax")
#' exp <- guesst(d = doses[4], p = 0.1, model="exponential", Maxd = doses[nd])
#' logit <- guesst(d = c(doses[3], doses[4]), p = c(0.1,0.8), "logistic",  Maxd= doses[nd])
#' betam <- guesst(d = doses[2], p = 0.3, "betaMod", scal=120, dMax=50, Maxd= doses[nd])
#' 
#' models.candidate <- Mods(emax = emax, linear = NULL,
#'                          exponential = exp, logistic = logit,
#'                          betaMod = betam, doses = doses,
#'                          placEff = placEff, maxEff = (maxEff- placEff))
#' plot(models.candidate)
#' 
#' ## True Model
#' model.true <- "emax"
#' response <- model_response(doses = doses,
#'                            distr = "weibull", 
#'                            model.true = model.true, 
#'                            models.candidate = models.candidate) 
#' lambda.true <- response$lambda
#' parm <- list(lambda = lambda.true, sigma = sigma.true)
#' 
#' ## Scenario: Censoring 10%
#' censoring.rate <- 0.1
#' 
#' dt <- data_generator(doses = doses,
#'                      sample.size = 20,
#'                      distr = "weibull",
#'                      parm = parm,                    
#'                      censoring.rate = censoring.rate)
#' }
#'
#'@importFrom stats quantile runif
#' 
#'@export
data_generator <-
  function(doses,
           sample.size,
           distr,
           parm,
           censoring.rate = NULL){
    
    if (distr == "weibull"){
      
      if (is.list(parm) & all(names(parm) == c("lambda", "sigma"))){
        lambda.true <- parm$lambda
        sigma.true <- parm$sigma
      } else {
        stop("If distr = 'weibull', the argument `parm` should be a named list with true values for lambda and sigma.")
      }
      
      if (is.null(censoring.rate)){
        stop("If distr = 'weibull', the argument `censoring.rate` should be informed.")
      }
      
      lambda <- rep(lambda.true, each = sample.size)
      n <- length(lambda)
      t <- as.numeric(lambda * (-log(runif(n))) ^ sigma.true)
      censor <- as.numeric(rep(quantile(t, 1 - censoring.rate), n))
      
      response <- pmin(t, censor)
      status <- ifelse(t > censor, 0, 1)
      
      out <- data.frame(time = response, status, censor = censor, 
                        doses = as.factor(rep(doses, each = sample.size)))
      
    }
    return(out)
  }
