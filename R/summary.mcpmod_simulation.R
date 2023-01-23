#' Summary of simulation results
#' 
#' @description It summarizes results of simulations of dose-finding trials 
#' following the MCP-Mod approach with bias-corrected and second-order 
#' covariance matrices.
#' 
#' @param object an object of the "mcpmod_simulation" class.
#' @param ... additional arguments affecting the summary produced.
#' 
#' @examples 
#' \donttest{
#' library(DoseFinding)
#' library(MCPModBC)
#' 
#' ## doses scenarios 
#' doses <- c(0, 5, 25, 50, 100)
#' nd <- length(doses)
#' sample.size <- 25
#' 
#' # shape parameter
#' sigma.true <- 0.5
#' 
#' # median survival time for placebo dose
#' mst.control <- 4 
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
#' significance.level <- 0.05
#' selModel <- "AIC"
#' 
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
#' ## Simulation Parameters
#' n.sim <- 10
#' seed <- 1234 
#' n.cores <- 1
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
#' test <- mcpmod_simulation(doses = doses,
#'            parm = parm,
#'            sample.size = sample.size,
#'            model.true = model.true,
#'            models.candidate = models.candidate,
#'            selModel = selModel,
#'            significance.level = significance.level,
#'            Delta = Delta,
#'            distr = "weibull",
#'            censoring.rate = censoring.rate,
#'            sigma.estimator = "jackknife",
#'            n.cores = n.cores, seed = seed, n.sim = n.sim)
#' summary(test)
#' }
#' 
#' @importFrom dplyr mutate group_by summarise
#' @importFrom rlang .data
#' 
#' @export 
summary.mcpmod_simulation <- function(object, ...){
  
  df <- rbind(object$mle, 
              object$mle2, 
              object$bce, 
              object$bce2, 
              object$firth)
  
  out <- df |> 
    mutate(bias = .data$estimated.td - object$conditions$td.true,
           mse = (.data$estimated.td - object$conditions$td.true)^2) |> 
    group_by(.data$method, .data$sample.size) |>
    summarise(convergence_rate = 100*mean(!is.na(.data$poc)),
              poc = mean(.data$poc, na.rm = TRUE), 
              ps = mean(.data$ps, na.rm = TRUE), 
              estimated_td = mean(.data$estimated.td, na.rm = TRUE),
              bias = mean(.data$bias, na.rm = TRUE),
              rmsq = sqrt(mean(.data$mse, na.rm = TRUE)))

  return(data.frame(out))
}