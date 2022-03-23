#' A function for analyzing two-sample problems
#' 
#' The \code{rank.two.samples()} function calculates purely nonparametric rank-based 
#' methods for the analysis of two independent samples. Specifically, it
#' implements the Brunner-Munzel test and its generalizations for the Nonparametric Behrens-Fisher Problem,
#' that is, testing whether the relative effect p=P(X<Y)+1/2*P(X=Y) of the two independent samples X and Y 
#' is equal to 1/2.  Range preserving confidence intervals (and corresponding test statistics)
#' are available using Logit or Probit transformations. The function also implements studentized permutation 
#' tests and permutation based confidence intervals for p using any of the method above (see the details below). 
#' Furthermore, the Wilcoxon-Mann-Whitney test (exact and asymptotic) can be used to test the equality of
#' the two distribution functions of the two samples. The user can specify whether confidence intervals for shift
#'  effects shall be computed. The \code{rank.two.samples()} function implements one-sided and two-sided tests 
#'  and confidence intervals. You can plot the confidence intervals (for the relative
#' effects) with the \code{plot()} function. 
#' 
#' @param formula A model \code{\link{formula}} object. The left hand side
#'    contains the response variable and the right hand side contains the factor
#'    variables of interest. 
#' @param data A data.frame, list or environment containing the variables in 
#'    \code{formula}. The default option is \code{NULL}.
#' @param conf.level A number specifying the confidence level; the default is 0.95.
#' @param alternative A character string specifying the alternative hypothesis. One of "two.sided", "less", "greater". You can specify just the initial letter.
#' @param rounds Value specifying the number of digits the results are rounded to. Default is 4 decimals.
#' @param method specifying the method used for calculation of the confidence intervals.
#'    One of "t.app", "logit", "probit" or "normal".

#' @param permu A logical variable indicating whether you want to compute a studentized permutation test.
#' @param info A logical variable. Here, info = FALSE suppresses the output of additional information
#'    concerning e.g. the interpretation of the test results.
#' @param wilcoxon asymptotic or exact calculation of Wilcoxon test.
#' @param shift.int Logical, indicating whether or not shift effects should be considered.
#' @param nperm Number of permutations used, default is 10000.
#' 
#' @author Frank Konietschke
#' Brunner, E., Bathke, A. C., Konietschke, F. (2018). Rank and Pseudo-Rank Procedures for Independent Observations in Factorial Designs. Springer International
#' Publishing.
#' 
#' @references Brunner, E. and Munzel, U. (2000). The nonparametric Behrens-Fisher problem: Asymptotic
#' theory and a small-sample approximation. Biometrical Journal 1, 17 - 21.
#' 
#' @references Kaufmann, J., Werner, C., and Brunner, E. (2005). Nonparametric methods for analysing the
#' accuracy of diagnostic tests with multiple readers. Statistical Methods in Medical Research 14, 129 - 146
#' 
#' @references Pauly, M., Asendorf, T.,  Konietschke, F. (2016). Permutation-based inference for the AUC: a unified approach for continuous and discontinuous data.##' Biometrical Journal, 58(6), 1319 -- 1337.
#'
#' @details  The \code{rank.two.samples()} function calculates both transformed (logit or probit) and untransformed statistics 
#' (normal or t.app) for testing the null hypothesis p=1/2. If a studentized permutation test is performed, then the
#' permutation distribution of the respective statistics are computed, see Pauly et
#' al.(2016) for details. In any case, the function reports the point estimator and its estimated standard error, 
#' value of the test statistic, confidence interval and p-value. In case of separated samples, point estimator and standard error
#' would be 0 and thus, test statistics would not be defined. In such a case, point
#' estimator and its standard error are replaced by the numbers one would obtain if samples overlapped in a single point. 
#' A plot of the confidence interval can be obtained with the plot function.

#' 
#' @examples
#' data(Muco)  
#' Muco2 <- subset(Muco, Disease != "OAD")
#' Muco2$Disease <- droplevels(Muco2$Disease)
#'
#' twosample <- rank.two.samples(HalfTime ~ Disease, data = Muco2, 
#' wilcoxon = "exact", permu = TRUE, shift.int = TRUE, nperm = 1000)
#' twosample <- rank.two.samples(HalfTime ~ Disease, data = Muco2, 
#'    alternative = "greater", method = "probit", wilcoxon = "exact", permu = TRUE,
#'    shift.int = FALSE, nperm = 1000)
#' plot(twosample)
#' 
#' @seealso \code{\link{rankFD}}
#' 

#' @importFrom coin wilcox_test pvalue statistic
#' @export


rank.two.samples<-function (formula, data, conf.level = 0.95, 
alternative = c("two.sided", "less", "greater"), rounds = 4, method = c("t.app","logit", "probit", 
"normal"), permu=TRUE,  info = TRUE, 
wilcoxon=c("asymptotic","exact"),shift.int=TRUE, nperm = 10000) 
{  
alpha <- 1 - conf.level
if (alpha >= 1 || alpha <= 0){stop("The confidence level must be between 0 and 1!")
if (is.null(alternative)){stop("Please declare the alternative! (two.sided, less, greater)")}}
alternative <- match.arg(alternative)
method <- match.arg(method)
wilcoxon <- match.arg(wilcoxon)
if (length(formula) != 3) {stop("You can only analyse one-way layouts!")}
dat <- model.frame(formula, droplevels(data))
if (ncol(dat) != 2) {
   stop("Specify one response and only one class variable in the formula")}
if (is.numeric(dat[, 1]) == FALSE) {stop("Response variable must be numeric")}
response <- dat[, 1]
factorx <- as.factor(dat[, 2])
fl <- levels(factorx)
a <- nlevels(factorx)
if (a > 2){stop("You want to perform a contrast test (the factor variable has more than two levels)!")}
samples <- split(response, factorx)
n <- sapply(samples, length)
n1 <- n[1]
n2 <- n[2]
if (any(n == 1)){stop(paste("The factor level", fl[n == 1], "has only one observation!"))}
N <- sum(n)
data.info <- data.frame(Sample = fl, Size = n)
cmpid <- paste("p(", fl[1], ",", fl[2], ")", sep = "")
   
#------------------Compute the values of the Test Statistics-------------------#
  
BMstats <- BMstat(samples[[1]],samples[[2]],n1,n2,method)

#-----------------Compute Test for Relative Effects----------------------------#
crit1 <- switch(alternative,two.sided={qnorm(1-alpha/2)},less={qnorm(1-alpha)},greater={qnorm(1-alpha)})
crit2 <- switch(alternative,two.sided={qnorm(alpha/2)},less={qnorm(1-alpha)},greater={qnorm(1-alpha)})
ci.limits <- limits2(BMstats[,6],crit1,crit2,BMstats[,7],alternative,method)
pvalues=switch(alternative,two.sided={2*min(pnorm(BMstats[,5]),1-pnorm(BMstats[,5]))},
less={pnorm(BMstats[,5])},greater={1-pnorm(BMstats[,5])})
if(method=="t.app"){
pvalues=switch(alternative,two.sided={2*min(pt(BMstats[,5],BMstats[,9]),1-pt(BMstats[,5],BMstats[,9]))},
less={pt(BMstats[,5],BMstats[,9])},greater={1-pt(BMstats[,5],BMstats[,9])})
crit1 <- switch(alternative,two.sided={qt(1-alpha/2,BMstats[,9])},less={qt(1-alpha,BMstats[,9])},greater={qt(1-alpha,BMstats[,9])})
crit2 <- switch(alternative,two.sided={qt(alpha/2,BMstats[,9])},less={qt(1-alpha,BMstats[,9])},greater={qt(1-alpha,BMstats[,9])})
ci.limits <- limits2(BMstats[,6],crit1,crit2,BMstats[,7],alternative,method)}  
Rel.Effects <- data.frame(Effect=cmpid, Estimator=BMstats[,1], Std.Error=BMstats[,10], 
T=BMstats[,5],Lower=ci.limits[,1],Upper=ci.limits[,2],p.Value=pvalues)
Rel.Effects[,2:7] <- round(Rel.Effects[,2:7],rounds)
rownames(Rel.Effects)<-""

#-------------------Compute the Studentized Permutation Tests------------------#
if(permu){
permustat <-c()
Tpermu<-sapply(1:nperm,function(arg){
   resperm <- sample(response)
   BMstat(resperm[1:n1],resperm[(n1+1):N],n1,n2,method)[,5]})
p.permu <- switch(alternative, two.sided={2*min(mean(Tpermu<=BMstats[,5]),mean(Tpermu>=BMstats[,5]))},
  less={mean(Tpermu<=BMstats[,5])}, greater={mean(Tpermu>=BMstats[,5])})
crit1permu <-switch(alternative,two.sided={quantile(Tpermu,1-alpha/2)},less={quantile(Tpermu,(1-alpha))},greater={quantile(Tpermu,(1-alpha))})
crit2permu <- switch(alternative,two.sided={quantile(Tpermu,alpha/2)},less={quantile(Tpermu,(1-alpha))},greater={quantile(Tpermu,(1-alpha))})
ci.limits.permu <- limits2(BMstats[,6],crit1permu,crit2permu,BMstats[,7],alternative,method)
Stud.Perm <- data.frame(Effect=cmpid, Estimator=BMstats[,1], Std.Error=BMstats[,10], T=BMstats[,5],
Lower=ci.limits.permu[,1],Upper=ci.limits.permu[,2],p.Value=p.permu)
Stud.Perm[,2:7] <- round(Stud.Perm[,2:7],rounds)
rownames(Stud.Perm)<-""}
if(!permu){Stud.Perm<-NULL}


#-----------------------------Wilcoxon-Mann-Whitney ---------------------------#
alternativeakt <- switch(alternative,two.sided={"two.sided"}, less={"greater"}, greater="less")
Wilcox <- wilcox_test(response~factorx,distribution=wilcoxon,alternative=alternativeakt,conf.int=TRUE,conf.level=(1 - alpha))
if(wilcoxon=="asymptotic"){
WilcoxonTest <- data.frame( Effect=cmpid, Estimator=BMstats[,1], Std.Error=sqrt(BMstats[,8]/N),
Statistic=-1*statistic(Wilcox), p.Value=round(pvalue(Wilcox),rounds))
WilcoxonTest[,2:5]<-round(WilcoxonTest[,2:5],rounds)}
if(wilcoxon=="exact"){
WilcoxonTest <- data.frame( Effect=cmpid, Estimator=BMstats[,1],
Statistic=BMstats[,11], p.Value=round(pvalue(Wilcox),rounds))
WilcoxonTest[,2:4]<-round(WilcoxonTest[,2:4],rounds)}
rownames(WilcoxonTest)<-"" 

#------------------------------------Shift Effects-----------------------------#
if(shift.int){
HL.help=expand.grid(samples[[1]],samples[[2]])
effect.wilcoxon=median(HL.help[,2]-HL.help[,1])
shiftint=sort(-1*Wilcox@confint(1-alpha)$conf.int)
Lower.Shift=shiftint[1]
Upper.Shift=shiftint[2]
cmpidWilcoxon <- paste("delta","(",fl[2], "-", fl[1], ")", sep = "")
delta.interpretation="delta(.): Hodges-Lehmann Estimator"
shift.result <- data.frame(Effect=cmpidWilcoxon, Estimator=effect.wilcoxon,  Lower=Lower.Shift,Upper=Upper.Shift)
shift.result[,2:4]<-round(shift.result[,2:4],rounds)
rownames(shift.result)<-""}
if(!shift.int){
Lower.Shift = NA
Upper.Shift=NA
effect.wilcoxon = BMstats[,1]
cmpidWilcoxon <- paste("p(", fl[1], ",", fl[2], ")", sep = "")  
delta.interpretation=NA
shift.result<-NA}

#---------------------------------RESULT and Output----------------------------#
result <-list(Call=formula,Descriptive=data.info, Analysis=Rel.Effects, Studentized.Permutation=Stud.Perm, Wilcoxon=WilcoxonTest, Shift.Effects=shift.result)
output.alternative=switch(alternative, two.sided={"Relative Effect is unequal to 1/2"},less={"Relative Effect is less than 1/2"},
greater={"Relative Effect is greater than 1/2"})
output.interpretation=paste("If", cmpid, ">1/2, then data in group",fl[2],"tend to be larger than those in group",fl[1])
output.method=switch(method,normal={"Normal Approximation"},logit={"Logit Transformation"},probit={"Probit Transformation"},
t.app={paste("T-Approximation with", round(BMstats[,9],rounds), "DF")})
#-------------------------Prepare Plot and Print Tables------------------------#
result$plotting <- list(method=method, Effects=Rel.Effects, alpha=alpha)
result$output <-list(output.alternative=output.alternative,output.interpretation=output.interpretation,
output.method=output.method,permu=permu,nperm=nperm,info=info,wilcoxon=wilcoxon,alpha=alpha, delta.interpretation=delta.interpretation)
  
class(result) <- "ranktwosamples"
return(result)
}


