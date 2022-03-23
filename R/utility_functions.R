################################################################################
#
#Plot and Print Functions for rank.two.samples and rankFD Functions
#Author: Frank Konietschke and Sarah Friedrich
#
################################################################################
#Plot Functions
################################################################################

#' @export 
plot.rankFD <- function (x, ...) {
  
  object <- x
  dots <- list(...)
  a <- object$plotting
  # default values
  args <- list(nf = a$nf, fac_names = a$fac_names, n.hypotheses = a$n.hypotheses,
               Descriptive.Factors = a$Descriptive.Factors, CI.method = a$CI.method, alpha=a$alpha)
  args[names(dots)] <- dots
  do.call(plotting, args = args)}
  
  #' @export
plot.sci <- function (x, ...) {
  object <- x
  dots <- list(...)
  a <- object$plotsci
  # default values
  args <- list(resmctp=a$resmctp,hypothesis=a$hypothesis,sci.method=a$sci.method,alpha=a$alpha)
  args[names(dots)] <- dots
  do.call(plotting.sci, args = args)}
  
  #' @export
plot.ranktwosamples <- function (x, ...) {
  object <- x
  dots <- list(...)
  a <- object$plotting
  args <- list(method=a$method,Effects=a$Effects,cmpid=a$cmpid, alpha=a$alpha)
  args[names(dots)] <- dots
  do.call(plot2samples, args = args)}

################################################################################
#Print Functions
################################################################################

printmultisamples <-function(output.hypothesis,output.ranks,output.confidence,
                  output.info,output.contrast,output.sci.method, output.alpha)   {
   if(output.info){
   cat("\n \n               Nonparametric Methods for General Factorial Designs      \n \n ")
   cat("---------------------------------------------------------------------------\n")
   cat(" #Hypotheses: Tested in",output.hypothesis,"\n")
   cat(" #Ranking Method:", output.ranks,"\n")
   cat(" #Confidence Intervals:",(1-output.alpha)*100,"% with",output.confidence, "\n \n \n")
   if(!is.null(output.contrast) && output.sci.method=="multi.t"){
   cat(" #MCTP: Multivariate T Approximation \n")}
   if(!is.null(output.contrast) && output.sci.method=="fisher"){
   cat(" #MCTP: Fisher Transformation and multivariate T-Approximation\n")}
   cat("---------------------------------------------------------------------------\n") }}

#' @export
print.rankFD <- function (x, ...) {
  object <- x
  dots <- list(...)
  a <- object$output 
  args <- list(output.hypothesis=a$output.hypothesis,output.ranks=a$output.ranks,output.confidence=a$output.confidence,
  output.info=a$output.info,output.contrast=a$output.contrast,output.sci.method=a$output.sci.method,output.alpha=a$output.alpha)
  do.call(printmultisamples, args = args)
  cat("\n", "Call:", "\n",sep="")
  print(x$Call)
  cat("\n", "Descriptive:", "\n", sep = "")
  print(x$Descriptive)
  cat("\n", "Wald.Type.Statistic:", "\n", sep = "")
  print(x$Wald.Type.Statistic)
  cat("\n", "ANOVA.Type.Statistic:", "\n", sep = "")
  print(x$ANOVA.Type.Statistic)
  cat("\n", "Kruskal-Wallis Test:", "\n", sep = "")
  print(x$Kruskal.Wallis.Test)
  cat("\n", "MCTP:", "\n", sep = "")
  print(x$MCTP)
  cat("\n", "Covariance Matrix:", "\n", sep = "")
  print(x$Covariance.Matrix)
  cat("\n", "Factor.Information:", "\n", sep = "")
  print(x$Factor.Information)
  }


print2samples <-function(output.alternative,output.interpretation,output.method,permu,nperm,info,wilcoxon,alpha,delta.interpretation){

  if(info){
  cat("\n \n               Nonparametric Methods for 2 Independent Samples      \n \n ")
  cat("#Alternative:",output.alternative,"\n")
  cat(" #Method:", output.method,"\n")
  cat(" #Interpretation:", output.interpretation, "\n")
  cat(" #Confidence Level:",(1-alpha)*100,"%", "\n ")
  if(permu){cat("#Number of permutations:", nperm, "\n \n")}
  cat("#Wilcoxon-Mann-Whitney Test:", wilcoxon,"\n")
  cat(" #Shift-Effect:", delta.interpretation,"\n")
  cat("---------------------------------------------------------------------------\n")}}

#' @export
print.ranktwosamples <- function (x, ...) {
  object <- x
  dots <- list(...)
  a <- object$output 
  args <- list(output.alternative=a$output.alternative,output.interpretation=a$output.interpretation,
  output.method=a$output.method,permu=a$permu,nperm=a$nperm,info=a$info,wilcoxon=a$wilcoxon,alpha=a$alpha,delta.interpretation=a$delta.interpretation)
  args[names(dots)] <- dots
  do.call(print2samples, args = args)
  cat("\n", "Call:", "\n",sep="")
  print(x$Call)
  cat("\n", "Descriptive:", "\n", sep = "")
  print(x$Descriptive)
  cat("\n----------------------Analysis of Relative Effects-------------------------")
  cat("\n", "Test Results:", "\n", sep = "")
  print(x$Analysis)
  cat("\n", "Studentized Permutation Test:", "\n", sep = "")
  print(x$Studentized.Permutation)
  cat("\n\n-------------------Analysis of Distribution Functions-----------------------", "\n \n")
  cat("Wilcoxon-Mann-Whitney Test:","\n", sep = "")
  print(x$Wilcoxon)
  cat("\n", "Shift Effects:","\n", sep = "")
  print(x$Shift.Effects)
  cat("")}
