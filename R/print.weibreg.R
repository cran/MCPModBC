print.weibreg <-
function (x, digits = max(3L, getOption("digits") - 3L), 
    ...) 
{
            cat("Weibull regression model\n", 
                "Method to obtain estimators: maximum likelihood estimators\n","Method to reduce bias in regression coefficients: ",ifelse(x$estimator=="MLE","none\n",
ifelse(x$estimator=="BCE","Cox and Snell's method\n","Firth's method\n")),
"Method to correct variance in finite samples: ",ifelse(x$corrected.var,"yes\n","no\n"), sep = "")
    cat("Sample size: ",length(x$y))
    cat("\n")
    cat("\n")
    cat("Coefficients:\n")
    print.default(format(x$coefficients, digits = digits), print.gap = 2L, 
        quote = FALSE)
    cat("\n")
    cat("Scale: ",x$scale)
    cat("\n")
    invisible(x)
}
