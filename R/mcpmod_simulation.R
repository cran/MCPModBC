#' Simulation to obtain operating characteristics for MCP-Mod design
#' 
#' @description It simulates dose-finding trials using MCP-Mod design with Maximum Likelihood Estimator and Fisher Information (MLE), 
#' Maximum Likelihood Estimator and second-order Fisher Information (MLE2), Cox and Snell's Bias-Corrected Estimator and 
#' Fisher Information (BCE), Cox and Snell's Bias-Corrected Estimator and second-order Fisher Information (BCE2), and 
#' Firth Bias-Corrected estimators (Firth) as discussed in Diniz, Magalhães and Gallardo. 
#' 
#' @param doses a numeric vector indicating the doses that will be considered in the clincial trial.
#' @param parm a named list of true values for the simulation. See more details in \code{\link{data_generator}}.
#' @param sample.size a numeric vector indicating the sample sizes per dose in the clinical trial to be evaluated in the simulation study.
#' @param model.true a character value indicating the true functional form of dose-response curve. See more details in \code{\link{model_response}}.
#' @param models.candidate an object of class 'Mods'. See more details in \code{\link[DoseFinding]{Mods}}.
#' @param selModel a character value indicating the model selection criterion for dose estimation. See more details in \code{\link[DoseFinding]{MCPMod}}.
#' @param significance.level a numeric value indicating the significance level when evaluating proof-of-concept based on an one-sided Wald test.
#' @param Delta a numerical value indicating the target effect size used for the target dose. See \code{\link[DoseFinding]{TD}}.
#' @param censoring.rate a numeric value between 0 and 1 indicating the censoring rate when generated data. It is required when \code{distr = "weibull"}.
#' @param distr a character value indicating the distribution of the response variable. Currently, the only option available is 'weibull'.
#' @param sigma.estimator a character value indicating whether the estimator for sigma should be a maximum likelihood or 
#' jackknife estimator. It is required when \code{distr = "weibull"}. Options are "mle" and "jackknife".
#' @param n.cores a numeric value indicating the number of cores to be used in the 'simulation performed in parallel. 
#' Use parallel::detectCores() to check the number of cores available.
#' @param seed an integer value, containing the random number generator (RNG) state for random number generation.
#' @param n.sim a numerical value indicating the number of simulated trials.
#' 
#' @return \code{mle}: a matrix of dimension n.sim x 4 with results when using the MCP-Mod approach with MLE;
#' @return \code{mle2}: a matrix of dimension n.sim x 4 with results when using the MCP-Mod approach with MLE2;
#' @return \code{bce}: a matrix of dimension n.sim x 4 with results when using the MCP-Mod approach with BCE;
#' @return \code{bce2}: a matrix of dimension n.sim x 4 with results when using the MCP-Mod approach with BCE2;
#' @return \code{firth}: a matrix of dimension n.sim x 4 with results when using the MCP-Mod approach with Firth's estimator;
#' 
#' All matrices contain the following columns:
#' (1) the first column indicates whether proof-of-concept (1 = "yes", 0 = "no"), in other words, the p-value of Wald test was statistically significant;
#' (2) the second column indicates whether the true model was selected to estimate the dose-response curve (1 = "yes", 0 = "no") when proof-of-concept is demonstrated;
#' (3) the third column contains the estimated target dose;
#' (4) the fourth column contains the sample size considered in the trial.
#' 
#' @return \code{conditions}: a list containing the conditions of the simulation.
#' 
#' @author Diniz, M.A., Gallardo D.I., Magalhaes, T.M.
#' 
#' @references 
#' Bretz F, Pinheiro JC, Branson M. Combining multiple comparisons and modeling techniques in dose‐response studies. Biometrics. 2005 Sep;61(3):738-48.
#' 
#' Bornkamp B, Pinheiro J, Bretz F. MCPMod: An R package for the design and analysis of dose-finding studies. Journal of Statistical Software. 2009 Feb 20;29:1-23.
#' 
#' Pinheiro J, Bornkamp B, Glimm E, Bretz F. Model‐based dose finding under model uncertainty using general parametric models. Statistics in medicine. 2014 May 10;33(10):1646-61.
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
#' }
#' 
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom doRNG %dorng%
#' @importFrom survival survreg Surv
#' @importFrom DoseFinding TD MCPMod
#' @importFrom stats coef model.matrix vcov
#' 
#' @export
mcpmod_simulation <-
  function(doses, parm, sample.size,
           model.true, models.candidate, selModel = "AIC",
           significance.level = 0.025,
           Delta, 
           distr = "weibull",
           censoring.rate = NULL,
           sigma.estimator = NULL,
           n.cores, seed, n.sim){
  
  # Summary of parallel simulation
  comb <- function(...) {
    mapply('rbind', ..., SIMPLIFY=FALSE)
  }
  
  # Summary for simulations
  aux <-
    function(mcp.out, significance.level, model.true, 
             models.candidate){
      
      theta <- models.candidate[[model.true]]
      if (model.true == "constant")
        theta <- rep(attr(models.candidate, "placEff"), 3)

      if (all(!is.na(mcp.out)) & all(class(mcp.out) != "try-error")){
        
        p.values <- attr(mcp.out$MCTtest$tStat,"pVal")
        
        if (any(p.values < significance.level)){
          poc <- 1
          if (mcp.out$selMod == model.true){
            ps <- 1
          } else {
            ps <- 0
          }
          estimated.td <- as.numeric(mcp.out$doseEst[mcp.out$selMod])
          
          if (!is.na(estimated.td) & estimated.td > doses[length(doses)])
            estimated.td <- doses[length(doses)]
          
          if (!is.na(estimated.td) & estimated.td < doses[1])
            estimated.td <- doses[1]
          
        } else {
          poc <- 0
          ps <- NA
          estimated.td <- NA
        }
        
      } else {
        poc <- NA
        ps <- NA
        estimated.td <- NA
      }
      
      out <- c(poc, ps, estimated.td)
      names(out) <- c("poc", "ps", "estimated.td")
      
      return(out)
    }

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

    if (is.null(sigma.estimator)){
      warning("If distr = 'weibull', the argument `sigma.estimator` should be informed.")
    } else {
      if (sigma.estimator != "jackknife" & sigma.estimator != "mle") 
        stop("Please define the argument sigma.estimator as either `jackknife` or `model`.")
    }
    
    mu.index <- 1:length(lambda.true)
    cov.index <- 1:length(lambda.true)
  }
  
  # Serial Loop
  mle <- list()
  mle2 <- list()
  bce <- list()
  bce2 <- list()
  firth <- list()
  
  for (i in 1:length(sample.size)){
  
    # Parallel Loop
    registerDoParallel(n.cores)
    set.seed(seed)
    results <- foreach (i = 1:n.sim,
                             .combine='comb',
                             .multicombine=TRUE,
                             .errorhandling = "remove",
                             .packages = c("survival",
                                           "DoseFinding",
                                           "nleqslv")) %dorng% {
                                             # for(i in 1:n.sim){
  
                               dt <- data_generator(doses = doses,
                                                    sample.size = sample.size,
                                                    parm = parm,
                                                    distr = "weibull",
                                                    censoring.rate = censoring.rate)
                               
                               mcp <- list()
                               
                               surv.mod <- survreg(Surv(time, status) ~ doses - 1, 
                                                     data = dt, dist = "weibull")
                             
                               ### MLE
                               beta.mle <- coef(surv.mod)
                               
                               mu.hat <- beta.mle
                               S.hat <- vcov(surv.mod)
                               design.matrix <- model.matrix(surv.mod)
                               
                               # sigma as fixed, estimated via jackknife
                               if (sigma.estimator == "jackknife"){
                                 sigma.estimated <- 
                                   sigma.jackknife(y = log(dt$time),
                                                   x = design.matrix,
                                                   delta = dt$status)
                               } else {
                                 sigma.estimated <- surv.mod$scale
                               } 
                               
                               mcp$mle <- try(MCPMod(dose = doses[mu.index], 
                                                      resp = mu.hat[mu.index], 
                                                      models = models.candidate, 
                                                      S = S.hat[cov.index, 
                                                                cov.index], 
                                                      type = "general", 
                                                      placAdj = FALSE, 
                                                      selModel = selModel, 
                                                      alpha = significance.level, 
                                                      df = NULL, critV = NULL, 
                                                      doseType = "TD", 
                                                      Delta = Delta,  
                                                      pVal = TRUE,
                                                      alternative = "one.sided"))
  
                               ### Corrections
  
                               # MLE second-order covariance matrix
                               vcov.mle2 <- try(
                                 var_mle2(y = log(dt$time),
                                          x = design.matrix,
                                          L = dt$censor,
                                          delta = dt$status,
                                          sigma = sigma.estimated,
                                          beta.mle = beta.mle)
                               )
                               
                               if (!inherits(vcov.mle2, "try-error")){
  
                                 mu.hat <- beta.mle
                                 S.hat <- vcov.mle2$cov
                                 mcp$mle2 <- try(
                                   MCPMod(dose = doses[mu.index],
                                          resp = mu.hat[mu.index],
                                          models = models.candidate,
                                          S = S.hat[cov.index,
                                                    cov.index],
                                          type = "general",
                                          placAdj = FALSE,
                                          selModel = selModel,
                                          alpha = significance.level,
                                          df = NULL, critV = NULL,
                                          doseType = "TD",
                                          Delta = Delta,
                                          pVal = TRUE,
                                          alternative = "one.sided")
                                 )
                               } else {
                                 mcp$mle2 <- NA
                               }
  
                               # Cox estimator for beta
                               beta.cox <- try(
                                 beta_cox(y = log(dt$time),
                                          x = design.matrix,
                                          L = dt$censor,
                                          delta = dt$status,
                                          sigma = sigma.estimated,
                                          beta.mle = beta.mle)
                               )
  
                               if(!inherits(beta.cox, "try-error")){
                                   
                                   vcov.cox <- try(
                                     var_cox(y = log(dt$time),
                                             x = design.matrix,
                                             L = dt$censor,
                                             delta = dt$status,
                                             sigma = sigma.estimated,
                                             beta.cox = beta.cox,
                                             corrected = FALSE)
                                   )
                                   
                                if (!inherits(vcov.cox, "try-error")){
                                   mu.hat <- beta.cox
                                   S.hat <- vcov.cox$cov
  
                                   mcp$bce <- try(
                                     MCPMod(dose = doses[mu.index],
                                            resp = mu.hat[mu.index],
                                            models = models.candidate,
                                            S = S.hat[cov.index,
                                                      cov.index],
                                            type = "general",
                                            placAdj = FALSE,
                                            selModel = selModel,
                                            alpha = significance.level,
                                            df = NULL, critV = NULL,
                                            doseType = "TD",
                                            Delta = Delta,
                                            pVal = TRUE,
                                            alternative = "one.sided")
                                   )
                                } else {
                                  mcp$bce <- NA
                                }
                               } else {
                                 mcp$bce <- NA
                               }
                                   
                               if(!inherits(beta.cox, "try-error")){
                                   vcov.cox <- try(
                                     var_cox(y = log(dt$time),
                                             x = design.matrix,
                                             L = dt$censor,
                                             delta = dt$status,
                                             sigma = sigma.estimated,
                                             beta.cox = beta.cox,
                                             corrected = TRUE)
                                   )
  
                                if (!inherits(vcov.cox, "try-error")){
                                   mu.hat <- beta.cox
                                   S.hat <- vcov.cox$cov
  
                                   mcp$bce2 <- try(
                                     MCPMod(dose = doses[mu.index],
                                            resp = mu.hat[mu.index],
                                            models = models.candidate,
                                            S = S.hat[cov.index,
                                                      cov.index],
                                            type = "general",
                                            placAdj = FALSE,
                                            selModel = selModel,
                                            alpha = significance.level,
                                            df = NULL, critV = NULL,
                                            doseType = "TD",
                                            Delta = Delta,
                                            pVal = TRUE,
                                            alternative = "one.sided")
                                   )
  
                                 
                               } else { 
                                 mcp$bce2 <- NA
                               }
                              } else {
                                 mcp$bce2 <- NA
                              }
  
                               #Firth estimator
                               beta.firth <- try(beta_firth(y = log(dt$time),
                                                            x = design.matrix,
                                                            L = dt$censor,
                                                            delta = dt$status,
                                                            sigma = sigma.estimated,
                                                            beta.mle = beta.mle))
  
                               if(!inherits(beta.firth, "try-error")){
                                 vcov.firth <- try(var_firth(y = log(dt$time),
                                                             x = design.matrix,
                                                             L = dt$censor,
                                                             delta = dt$status,
                                                             sigma = sigma.estimated,
                                                             beta.firth = beta.firth))
  
                                 if (!inherits(vcov.firth, "try-error")){
  
                                   mu.hat <- beta.firth
                                   S.hat <- vcov.firth$cov
  
                                   mcp$firth <- try(MCPMod(dose = doses[mu.index],
                                                          resp = mu.hat[mu.index],
                                                          models = models.candidate,
                                                          S = S.hat[cov.index,
                                                                    cov.index],
                                                          type = "general",
                                                          placAdj = FALSE,
                                                          selModel = selModel,
                                                          alpha = significance.level,
                                                          df = NULL, critV = NULL,
                                                          doseType = "TD",
                                                          Delta = Delta,
                                                          pVal = TRUE,
                                                          alternative = "one.sided"))
  
                                 } else {
                                   mcp$firth <- NA
                               } 
                              } else {
                                 mcp$firth <- NA
                              }
                               
                               out <- lapply(mcp, aux,
                                             significance.level = significance.level,
                                             model.true = model.true,
                                             models.candidate = models.candidate)
                               
                             }
    stopImplicitCluster()
  
    
      
    mle[[i]] <- data.frame(results$mle, method = "mle", sample.size = sample.size[i])
    mle2[[i]] <- data.frame(results$mle2, method = "mle2", sample.size = sample.size[i])
    bce[[i]] <- data.frame(results$bce, method = "bce", sample.size = sample.size[i])
    bce2[[i]] <- data.frame(results$bce2, method = "bce2", sample.size = sample.size[i])
    firth[[i]] <- data.frame( results$firth, method = "firth", sample.size = sample.size[i])
  
  }
  
  mle <- Reduce(rbind, mle)
  mle2 <- Reduce(rbind, mle2)
  bce <- Reduce(rbind, bce)
  bce2 <- Reduce(rbind, bce2)
  firth <- Reduce(rbind, firth)
  
  td.true <- TD(models.candidate, Delta = Delta)[model.true]
  conditions <- list(doses = doses,
                     parm = parm, 
                     sample.size = sample.size,
                     model.true = model.true, 
                     models.candidate = models.candidate,
                     selModel = selModel,
                     significance.level = significance.level,
                     Delta = Delta, 
                     distr = distr,
                     censoring.rate = censoring.rate,
                     sigma.estimator = sigma.estimator,
                     n.cores = n.cores, 
                     seed = seed, 
                     n.sim = n.sim,
                     td.true = td.true)

  
  out <- list(mle = mle, mle2 = mle2, bce = bce, bce2 = bce2, firth = firth, 
              conditions = conditions)
  class(out) <- "mcpmod_simulation"
  
  return(out)
  
}
