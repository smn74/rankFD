VV <- NULL
#' Rank-based tests for general factorial designs
#' 
#' The function implements purely nonparametric rank-based methods for the analysis 
#' of general factorial designs. You can chose to use either classical ranks (mid-ranks) 
#' (\code{effect="weighted"}) or pseudo-ranks (\code{effect="unweighted"})
#' for making inference. Pseudo-ranks are used by default. 
#' The package implements point estimators of relative effects (weighted and unweighted) as
#' well as test procedures (Wald-Type and ANOVA-Type statistics) for testing global null hypotheses
#' formulated in either (i) distribution functions \code{hypothesis="H0F"} 
#' or (ii) relative effects \code{hypothesis="H0p"}. In case of one-way factorial 
#' designs, the function additionally computes 
#' the Kruskal-Wallis test either with ranks or pseudo-ranks. In addition, multiple 
#' contrast tests (and simultaneous confidence intervals)
#' for the main or interaction effects are implemented within the \code{contrast} 
#' statement. You can either choose from pre-defined contrasts (options see below) or
#'  you can provide your own user-defined contrast matrix. Both the
#' Fisher-transformation (\code{sci.method="fisher"}) as well as a multivariate
#' t-approximation (\code{sci.method="multi.t"}) are implemented.
#' The Fisher approximation is used by default. To visualize the results, you can plot
#' the simultaneous confidence intervals using the \code{plot.sci} function.
#' Furthermore, confidence interval plots for the main or interaction relative effects
#' (not simultaneous) are available within the \code{plot} function. 
#' 
#' @param formula       A model \code{\link{formula}} object. The left hand side
#'                      contains the response variable and the right hand side contains 
#'                      the factor variables of interest. An interaction term must be specified.
#' @param  data         A data.frame, list or environment containing the variables in 
#'                      \code{formula}. The default option is \code{NULL}.
#' @param alpha         A number specifying the significance level; the default is 0.05.
#' @param CI.method     Either "logit" or "normal", specifying the method used for
#'                      calculation of the confidence intervals.
#' @param effect        In case of weighted, then weighted (by sample sizes) relative 
#'                      effects are estimated using classical ranks (mid-ranks) of the data. 
#'                      Otherwise, in case of effect="unweighted", unweighted relative 
#'                      effects are estimated with pseudo-ranks. The default option is "unweighted" resulting in pseudo-rank statistics. 
#' @param hypothesis    The null hypothesis to be tested, either "H0F" or "H0p". The option "H0F" computes tests for testing hypotheses
#'                      formulated in terms of distribution functions. Otherwise, 
#'                      hypotheses in relative effects are tested. The latter allows for variance heteroscedasticity even
#'                      under the null hypothesis of no treatment effect and thus covers the Nonparametric Behrens-Fisher problem. 
#' @param Factor.Information Logical. If TRUE, descriptive statistics with point estimators, standard error as well as confidence intervals for 
#'                      each main and interaction effect in the model are printed. The results can furthermore be plotted with the plot function.
#' @param contrast      a list containing the name of the main or interaction effect (written as group1:group2), 
#'                      a pre-defined contrast ("Dunnett", "Tukey", "Sequen", "AVE", "Changepoint", 
#'                      "Williams", "Marcus", "McDermott", "UmbrellaWilliams", "GrandMean") or a user-defined contrast matrix. If the contrast coefficients do not
#'                      sum up to 0, or if their sum of absolutes differs from 2, the coefficients are normalized. 
#' @param sci.method    Either "fisher" or "multi.t" as approximation method for the multiple contrast tests and simultaneous confidence intervals. 
#'                      The default option is "fisher".
#' @param info          Logical. If TRUE, additional output information and explanation is printed to the console.
#' @param rounds        Number of  decimals of the output values. The default option is rounds=4 (4 decimals).
#' @param covariance    Logical. If TRUE, the estimated covariance matrix of the vector of relative effects is computed.
#'
#' @details 
#' The rankFD() function calculates the Wald-type statistic (WTS), ANOVA-type 
#' statistic (ATS) as well as multiple contrast tests for general factorial designs 
#' for testing the null hypotheses \eqn{H_0^F: CF = 0} or \eqn{H_0^p: Cp = 0}. 
#' Almost every method explained in the comprehensive textbook from Brunner et al. (2019) 
#' is implemented in rankFD. The test procedures for testing null hypotheses in distribution 
#' functions have initially been proposed by Akritas et al. (1997), whereas methods for testing null 
#' hypotheses formulated in relative effects have been proposed by Brunner et al. (2017). 
#' We note that the multiple contrast test procedure using Fisher approximation computes critical and 
#' p-values from a multivariate t-distribution with respective degrees of freedom. Simulation studies by
#' Konietschke et al. (2012) demonstrated an accurate control of the type-1 error rate and
#' the procedure is therefore recommended.
#
#'   
#' @return A \code{rankFD} object containing the following components:
#' \item{Call}{Given response and factor names (formula)}           
#' \item{Descriptive}{Descriptive statistics of the data for all factor
#'                    level combinations. Displayed are the number of individuals per factor
#'                    level combination (size), the relative effect (Rel.Effect), Standard Error and  100*(1-alpha)\% confidence
#'                    intervals.}
#' \item{WTS}{The value of the WTS along with degrees of freedom of the central chi-square distribution and p-value.}
#' \item{ATS}{The value of the ATS, degrees of freedom of the central F distribution and the corresponding p-value.}
#' \item{Kruskal-Wallis-Test}{The value of the Kruskal-Wallis test along with degrees of freedom and p-value. If effect="unweighted", the Kruskal-Wallis test
#'                    using pseudo-ranks is computed. Otherwise, if effect="unweighted", the "established" Kruskal-Wallis test based on ranks is returned.} 
#' \item{MCTP}{Contrast matrix, local Results in terms of point estimates, standard error, value of the test statistic, (1-alpha)100% simultaneous confidence
#'                    intervals as well as adjusted p-values. As a summary, the function also returns the global test decision by printing the maximum test
#'                    statistic (in absolute value) as well as the (1-alpha) critical value from the multivariate T-distribution.}
#' \item{Covariance.Matrix}{The estimated covariance matrix of the vector of the estimated relative effects. Note that the vector is multiplied by root N.}
#' \item{Factor.Information}{Descriptive tables containing the point estimators, standard errors as well as (1-alpha)100% confidence intervals for all possible main
#'                     and interaction effects in the model. The confidence intervals are not simultaneous and for data descriptive purpose only.}
#'
#'
#' @examples
#' data(Coal)
#' model <- rankFD(Acidity ~ NaOH * Type, data = Coal, CI.method = "normal",
#' effect = "unweighted", hypothesis = "H0F")
#' 
#' data(Muco)
#' model.oneway <- rankFD(HalfTime ~ Disease, data = Muco, CI.method = "logit",
#' effect = "weighted", hypothesis = "H0p")
#' plot(model.oneway)
#' 
#' 
#' @references 
#' Brunner, E., Bathke, A.C., Konietschke, F. Rank and Pseudo-Rank Procedures 
#' for Independent Observations in Factorial Designs. Springer International Publishing, 2018.
#' 
#' Brunner, E., Konietschke, F., Pauly, M., Puri, M. L. (2017). Rank-based procedures in factorial designs: 
#' Hypotheses about non-parametric treatment effects. Journal of the Royal Statistical Society: Series B 
#' (Statistical Methodology), 79(5), 1463-1485.
#' 
#' Akritas, M. G., Arnold, S. F., and Brunner, E. (1997). Nonparametric hypotheses and rank statistics for unbalanced factorial designs.
#' Journal of the American Statistical Association 92, 258-265.
#' 
#' Brunner, E., Dette, H., and Munk, A. (1997). Box-Type Approximations in Nonparametric Factorial Designs. Journal
#' of the American Statistical Association 92, 1494-1502.
#'
#'Konietschke, F., Hothorn, L. A., Brunner, E. (2012). Rank-based multiple test procedures and simultaneous confidence intervals. Electronic Journal of Statistics,
#' 6, 738-759.

#' 
#'
#' @importFrom lattice xyplot panel.superpose panel.arrows panel.points panel.xyplot panel.abline
#' @importFrom stats formula model.frame pchisq pf qnorm terms var aggregate as.formula confint cov median pnorm pt qt quantile sd cov2cor
#' @importFrom utils read.table
#' @importFrom MASS ginv
#' @importFrom coin wilcox_test pvalue statistic
#' @importFrom graphics abline axis box plot points polygon title
#' @importFrom multcomp contrMat 
#' @importFrom mvtnorm qmvnorm pmvt qmvt pmvnorm
#' 
#' @export



rankFD <- function(formula, data,alpha=0.05, CI.method=c("logit","normal"),
                     effect=c("unweighted","weighted"),hypothesis=c("H0F","H0p"),
                    Factor.Information=FALSE, contrast = NULL, 
                   sci.method=c("fisher", "multi.t"),info = TRUE,covariance=FALSE, rounds=4){
  
  effect=match.arg(effect)
  hypothesis = match.arg(hypothesis)
  CI.method = match.arg(CI.method)
  sci.method = match.arg(sci.method)
  #-----------------------Determine the model----------------------------------#
  
  dat.Model0 <- model.frame(formula, data)
  
  #------------------------Numbers of factors----------------------------------#
  
  for(i in length(dat.Model0):2){ dat.Model0[,i] <-as.factor(dat.Model0[,i]) }
  nf = ncol(dat.Model0)-1
  n.levels = c()
  names.levels=list()
  for (i in 2:(nf+1)){n.levels[i-1]= nlevels(dat.Model0[,i])
                      names.levels[[i-1]] = levels(dat.Model0[,i])}
  names(names.levels) <- names(dat.Model0[,2:(nf+1)])
  #-----------------------Hypotheses matrices----------------------------------#

  perm_names <- t(attr(terms(formula), "factors")[-1, ])
  nr_hypo <- attr(terms(formula), "factors")
  fac_names <- colnames(nr_hypo)
  Hypotheses0 = HC(n.levels,"Hyp",perm_names,fac_names)
  Hypotheses=Hypotheses0[[1]]
  CI.Matrices = HC(n.levels,"CI",perm_names,fac_names)[[1]]
  n.hypotheses = length(Hypotheses)
  n.levels.prod=prod(n.levels)
  Output.names <- Hypotheses0[[2]]
  #--------------------Sort Data according to Factors--------------------------#
  for (i in length(dat.Model0):2) {
    dat.Model0 <- dat.Model0[order(dat.Model0[,i]),]}
  #-----------------Introduce Pseudo Factor -----------------------------------#
  dat.Model0$Interaction = interaction(dat.Model0[,2:length(dat.Model0)],sep=":")
  dat.response <- dat.Model0[,1]
  #---------------------Sizes and Factor Constellations------------------------#
  n <- aggregate(formula,data=dat.Model0,length)
  for(i in (length(n)-1):1) {n <-n[order(n[,i]),]}
  colnames(n)[ncol(n)]<-"Size"
  #---------------------------Compute Inference Methods------------------------#
  dat.Model0$INum <- as.factor(rep(1:n.levels.prod, n$Size))
  #---------------------Compute the Estimators---------------------------------#
  H0pW <-Effects(dat.response, dat.Model0$INum,effect)
  WTS = WTSp=matrix(0,n.hypotheses,3)
  ATS = matrix(0,n.hypotheses, 4)
  ATSp = matrix(0,n.hypotheses,4)
  KW = matrix(0,1,3)
  Descriptive.Factors = list()
  Levels.Factors = list()
  if(n.hypotheses==1){KW[1,] = KWTEST(c(H0pW$pd), H0pW$varKW, n$Size )}
   for(i in 1:n.hypotheses){
    WTS[i,] = Wald(c(H0pW$pd),Hypotheses[[i]],H0pW$VH0F)
    WTSp[i,] = Wald(c(H0pW$pd),Hypotheses[[i]],H0pW$VBF/sum(n$Size))
    ATS[i,] =ANOVATYP(c(H0pW$pd),Hypotheses[[i]],H0pW$VH0F,n$Size)
    ATSp[i,]= ANOVATYPH0P(c(H0pW$pd),Hypotheses[[i]],H0pW$VBF,n$Size,H0pW$dfATS)
    CILimits <-Limits(c(CI.Matrices[[i]]%*%c(H0pW$pd)),CI.Matrices[[i]]%*%H0pW$VBF%*%t(CI.Matrices[[i]]),alpha,H0pW$N,CI.method)
    Descriptives <-data.frame(Rel.Effect=round(CI.Matrices[[i]]%*%H0pW$pd,rounds),
                              Std.Error= round(sqrt(c(diag(CI.Matrices[[i]]%*%H0pW$VBF%*%t(CI.Matrices[[i]])))/H0pW$N),rounds),
                              Lower=round(CILimits[,1],rounds),
                              Upper = round(CILimits[,2],rounds))
    Output.namesi <-Output.names[i] 
    formula.act <- as.formula(paste(names(dat.Model0)[1], Output.namesi, sep=" ~ "))
    aha <- data.frame(aggregate(formula.act,data=dat.Model0,mean))
    for(ii in (length(aha)-1):1) {aha <-aha[order(aha[,ii]),]}
    Descriptive.Factors[[i]] <-data.frame(aha,Descriptives)
    Descriptive.Factors[[i]] <-Descriptive.Factors[[i]][,-length(aha)]
    if(length(grep(":", Output.names[i]))<1){
      pos <- which(names(dat.Model0)==Output.names[i])
      Levels.Factors[[i]] = data.frame(X=levels(dat.Model0[,pos]))}
    if(length(grep(":", Output.names[i]))>=1){
      facs.singles <- c(strsplit(Output.names[i], ":")[[1]])
      Levels.Factors[[i]] = data.frame(n[,facs.singles])}
  }
  
  names(Descriptive.Factors) <- Output.names
  rownames(WTS) <- Output.names
  rownames(WTSp) <- Output.names
  rownames(ATS) <- Output.names
  rownames(ATSp)<- Output.names
  colnames(WTS) <- c("Statistic", "df", "p-Value")
  colnames(WTSp) <- c("Statistic", "df", "p-Value")
  colnames(ATS) <- c("Statistic", "df1", "df2", "p-Value")
  colnames(ATSp) <- c("Statistic", "df1", "df2", "p-Value")
  #---------------------Compute MCTP for Main Effects--------------------------#
 if(is.null(contrast)){res.mctp <- NULL}
 if(!is.null(contrast)){
 if(class(contrast)!="list"){stop("The contrast statement must be given as a list!")}
 if(!any(Output.names==contrast[[1]])){stop("Spelling error in the contrast statement. The factor is not contained in the data set. Please check.")}
 pos.contr=which(Output.names==contrast[[1]])
 
 if(length(contrast)==1){
 CC <- diag(nrow(Hypotheses[[pos.contr]]))
 switch(hypothesis,H0p={
 res.mctp <- MCTP(H0pW$pd,Hypotheses[[pos.contr]],H0pW$VBF,CC,n$Size,H0pW$placements,alpha,sci.method)},
   H0F = {res.mctp <- MCTP.H0F(H0pW$pd,Hypotheses[[pos.contr]],H0pW$VH0F,CC,n$Size,alpha,sci.method)})}
 if(length(contrast)==2){

if(any(class(contrast[[2]])=="character")){

 CC <- contrMat(n=rep(10,nrow(CI.Matrices[[pos.contr]])),contrast[[2]]) 
 switch(hypothesis,H0p={
 res.mctp <- MCTP(H0pW$pd,CI.Matrices[[pos.contr]],H0pW$VBF,CC,n$Size,H0pW$placements,alpha,sci.method)},
 H0F={res.mctp <- MCTP.H0F(H0pW$pd,CI.Matrices[[pos.contr]],H0pW$VH0F,CC,n$Size,alpha,sci.method)})}  
 if(any(class(contrast[2])!="character")){

 if(any(class(contrast[[2]])=="numeric")){
 if(length(contrast[[2]])!=nrow(CI.Matrices[[pos.contr]]) ){
 stop(paste("The number of coefficients in the contrast statement differs from the number of levels from the
 factor", contrast[[1]],". Please check."))}
  CC <- matrix(contrast[[2]],ncol=nrow(CI.Matrices[[pos.contr]])) 
}

 if(any(class(contrast[[2]])%in%"matrix")){
 CC<-contrast[[2]]   }
if(any(rowSums(CC)!=0)){
cpos <- which(rowSums(CC)!=0)
CC[cpos,]<-CC[cpos,]-rowMeans(CC)[cpos]
}

if(any(rowSums(abs(CC))!=2)){

rows_crit <- which(rowSums(abs(CC))!=2)

for(arg in 1:length(rows_crit)){
posi_neg <-which(c(CC[rows_crit[arg],])<0)
CC[rows_crit[arg], posi_neg]<-CC[rows_crit[arg], posi_neg]/sum(abs(CC[rows_crit[arg], posi_neg]))
posi_pos <-which(c(CC[rows_crit[arg],])>0)
CC[rows_crit[arg], posi_pos]<-CC[rows_crit[arg], posi_pos]/sum(abs(CC[rows_crit[arg], posi_pos]))

}

}

 }
 
if(covariance){
switch(hypothesis,H0F={VV <- H0pW$VH0F*H0pW$N},H0p={VV <- H0pW$VBF})}
if(!covariance){VV<-NULL}

 switch(hypothesis,H0p={
 res.mctp <- MCTP(H0pW$pd,CI.Matrices[[pos.contr]],H0pW$VBF,CC,n$Size,H0pW$placements,alpha,sci.method) },
 H0F={res.mctp <- MCTP.H0F(H0pW$pd,CI.Matrices[[pos.contr]],H0pW$VH0F,CC,n$Size,alpha,sci.method) })}
    colnames(res.mctp[[1]]) <- 1:ncol(res.mctp[[1]])
   rownames(res.mctp[[1]]) <- paste("C",1:nrow(CC),sep="")
   rownames(res.mctp[[2]]) <- paste("C",1:nrow(CC),sep="")
   res.mctp[[2]]<-round(res.mctp[[2]],rounds)
   res.mctp[[3]]<-round(res.mctp[[3]],rounds)
   res.mctp[[5]]<-round(res.mctp[[5]],rounds)} 
  #-----------------------------Descriptive Output------------------------------#
  n$Rel.Effect <- round(c(H0pW$pd),rounds)
  n$Std.Error <- round(sqrt(c(diag(H0pW$VBF))/sum(n$Size)),rounds)
  CI <- Limits(c(H0pW$pd),H0pW$VBF,alpha,c(H0pW$N),CI.method)
  n$Lower<- round(CI[,1],rounds);n$Upper <- round(CI[,2],rounds)
  rownames(n) <- 1:nrow(n)
  if (Factor.Information){res.factor.information <-Descriptive.Factors}
  if (!Factor.Information){res.factor.information <-NULL}

  if(hypothesis=="H0F"){
  if (nf ==1){
  colnames(KW) <- c("Statistic", "df", "p-Value")
  rownames(KW) <- Output.names
 result <- list(Call=formula,Descriptive=n, Wald.Type.Statistic = round(WTS,rounds), ANOVA.Type.Statistic=round(ATS,rounds), 
 Kruskal.Wallis.Test = round(KW,rounds),  MCTP=res.mctp, Factor.Information=res.factor.information, Covariance.Matrix=VV)}
 if(nf!=1){result <- list(Call=formula,Descriptive=n, Wald.Type.Statistic =round(WTS,rounds), ANOVA.Type.Statistic=round(ATS,rounds), MCTP=res.mctp,Factor.Information=res.factor.information, Covariance.Matrix=VV)}}
  if(hypothesis=="H0p"){
    result <- list(Call=formula,Descriptive=n, Wald.Type.Statistic=round(WTSp,rounds), ANOVA.Type.Statistic=round(ATSp,rounds), MCTP=res.mctp, Factor.Information=res.factor.information,Covariance.Matrix=VV)} 
    result$plotting <- list(nf = nf, fac_names = fac_names, n.hypotheses = n.hypotheses,
                          Descriptive.Factors = Descriptive.Factors, CI.method = CI.method, alpha=alpha)
    result$plotsci <- list(resmctp=res.mctp,hypothesis=hypothesis,alpha=alpha, sci.method=sci.method)
   
#---------------------------------RESULT and Output----------------------------#
 output.hypothesis=switch(hypothesis, H0F={"Distribution Functions"}, H0p={"Relative Effects"})
 output.ranks=switch(effect, weighted={"Global Ranks"}, unweighted={"Pseudo-Ranks"})
 output.confidence=switch(CI.method, logit={"Logit-Transformation"}, normal={"Normal Approximation"})

 
 result$output <- list(output.hypothesis=output.hypothesis,output.ranks=output.ranks,output.confidence=output.confidence,
 output.info=info,output.contrast=contrast,output.sci.method=sci.method,output.alpha=alpha)


class(result) <- "rankFD"
  return(result)
}


