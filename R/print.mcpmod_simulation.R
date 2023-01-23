print.mcpmod_simulation <-
function (x, digits = max(3L, getOption("digits") - 3L), 
    ...) 
{
x=x$conditions
hh=function(aa){
bb=paste(aa[1],sep="")
for(i in 2:(length(aa)-1))
{
	bb=paste(bb,", ",aa[i],sep="")
}
bb=paste(bb," and ",aa[length(aa)],sep="")
bb
}
 cat("Obtaining operating characteristics for MCP-Mod design") 
cat("\n------")
   cat("\ndoses: ",hh(x$doses),"\nsample size: ",x$sample.size,
"\ntrue model :",x$model.true,"\nMethod to select model: ",x$selModel,"\ndistribution of the response variable: ",
x$distr,"\ncensoring rate: ",x$censoring.rate,"\nsignificance level :",x$significance.level, sep = "")
    cat("\n")
    invisible(x)
}


