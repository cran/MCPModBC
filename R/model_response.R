#' Model Responser
#' 
#' It calculates the model response and parameters of interest for a given distribution.
#' 
#'@param doses a numeric vector indicating the doses that will be considered in the clinical trial.
#'@param distr a character value indicating the distribution of the response variable. Currently, the only option available is 'weibull'.
#'@param model.true a character value indicating the functional form of the true dose-response curve. Options are "constant", "linear", "linlog",  
#'"quadratic", "exponential", "emax", "sigmaEmax", "betaMod", "logistic",  "linInt".
#'@param models.candidate an object of class Mods. See more details in \code{\link[DoseFinding]{Mods}}.
#' 
#' @return a data frame with dimension length(doses) \eqn{\times} 3 with the following columns: (1) model response (2) model parameter and (3) doses
#' 
#' @author Diniz, M.A., Gallardo D.I., Magalhaes, T.M.
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
#' response
#' 
#' lambda.true <- response$lambda
#' parm <- list(lambda = lambda.true, sigma = sigma.true)
#' }
#' 
#' @importFrom DoseFinding quadratic exponential emax sigEmax betaMod logistic linear linlog linInt getResp
#' 
#' @export
model_response <-
  function(doses, distr, model.true, models.candidate){
    
 aux <- getResp(models.candidate, doses)
 mu <- aux[, model.true]

 if (distr == "weibull") {
   lambda <- exp(mu)
   out <- data.frame(mu, lambda, doses)
 }
    
  return(out)
}