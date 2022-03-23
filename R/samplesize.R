

#' Sample size computation methods for the Wilcoxon-Mann-Whitney test based on pilot data
#'
#' The function implements the sample size formula proposed by Happ et al. (see reference below). It estimates the sample size needed to detect the effect with
#' pre-defined power at significance level alpha based on pilot data.  
#'  
#' @param x1 advance information for the first group
#' @param x2 advance information for the second group
#' @param alpha two sided type I error rate
#' @param power power with the sample sizes of each group
#' @param t proportion of subjects in the first group
#' @return Returns a data frame
#' @examples
#' x1 = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2)
#' x2 = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
#' 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3)
#' WMWSSP(x1,x2,0.05,0.8,0.5)
#' @references Brunner, E., Bathke A. C. and Konietschke, F. Rank- and Pseudo-Rank Procedures in Factorial Designs - Using R and SAS. Springer Verlag. 
#' Happ, M., Bathke, A. C., & Brunner, E. (2019). Optimal sample size planning for the Wilcoxon-Mann-Whitney test. Statistics in medicine, 38(3), 363-375.
#' @export
WMWSSP=function(x1,x2,alpha=0.05,power=0.8,t=1/2){
stopifnot(all(is.finite(x1)), all(is.finite(x2)),alpha>0, alpha<1, power>0, t>0, t<1)
m1 <- length(x1)
m2 <- length(x2)
R <- rank(c(x1,x2), ties.method="average")
R1 <- R[1:m1]
R2 <- R[m1+(1:m2)]
R11 <- rank(x1, ties.method="average")
R22 <- rank(x2, ties.method="average")
P1 <- R1 - R11
P2 <- R2 - R22
pStar <- (mean(R2) - mean(R1))/(m1 + m2) + 0.5
sigmaStar <- sqrt(sum((R - (m1 + m2 + 1)/2)^2)/(m1 + m2)^3)
sigma1Star <- sqrt(sum((P1 - mean(P1))^2)/(m1 * m2^2))
sigma2Star <- sqrt(sum((P2 - mean(P2))^2)/(m1^2 * m2))
N <- (sigmaStar * qnorm(1 - alpha/2) + qnorm(power) * sqrt(t * 
        sigma2Star^2 + (1 - t) * sigma1Star^2))^2/(t * (1 - t)*(pStar - 0.5)^2)
n1 <- N*t
n2 <- N*(1-t)
output<-matrix(c(alpha,power,pStar,N,t,n1,n2,ceiling(n1)+ceiling(n2),ceiling(n1),ceiling(n2)),ncol=1)
rownames(output)<-c("alpha (2-sided)","Power", "Estimated relative effect p", "N (total sample size needed)", "t=n1/N", "n1 in Group 1", "n2 in Group 2", "N rounded", "n1 rounded", "n2 rounded")
colnames(output)<-""
return (output)}

#' Sample size calculation for the Wilcoxon-Mann-Whitney test using the Noether formula. The function estimates the sample size needed to detect the effect with
#' pre-defined power at significance level alpha using Noether's formula'.    
#'
#' @param x1 advance information is only needed in case of ties
#' @param alpha two sided type I error rate
#' @param power power: detect a relative effect p at least with probability power
#' @param p relative effect
#' @param t proportion of subjects in the first group (between 0 and 1)
#' @param ties TRUE if ties are possible (non continuous distribution), otherwise FALSE
#' @return Returns a data frame with the sample sizes for each group
#' @examples
#' noether(0.05,0.8,1/2, 0.75)
#' @references Noether, G. E. (1987). Sample Size Determination for Some Common Nonparametric Tests. Journal of the American Statistical Association 85, 645.647.
#' @export
noether=function(alpha,power,t, p, x1=c(0), ties=FALSE){
stopifnot(all(is.finite(x1)), alpha>0, alpha<1, power>0, t>0, t<1, p<1, p>0)
if(ties==TRUE & length(x1)<=1){stop("Advance information is needed in case of ties!")}
if(ties==FALSE & length(x1)>1){print("Advance information is not used in case of a continuous distribution!")}
# case 1: advance information x1, non continuous
m1 <- length(x1)
if(m1>1){
R <- rank(c(x1), ties.method="average")
R1 <- R[1:m1]
sigma2M <- 1/m1^3*sum((R1-(m1+1)/2)^2)
Nv <- sigma2M/(p-1/2)^2*(qnorm(1-alpha/2)+qnorm(power))^2*1/(t*(1-t))
n1v <- Nv*t
n2v <- Nv*(1-t)
output<-matrix(c(alpha,power,p,Nv,t,n1v,n2v,ceiling(n1v)+ceiling(n2v),ceiling(n1v),ceiling(n2v)),ncol=1)}
else{
Nu <- (qnorm(1-alpha/2)+qnorm(power))^2*1/(12*t*(1-t)*(p-1/2)^2)
n1u <- Nu*t
n2u <- Nu*(1-t)
output<-matrix(c(alpha,power,p,Nu,t,n1u,n2u,ceiling(n1u)+ceiling(n2u),ceiling(n1u),ceiling(n2u)),ncol=1)}
rownames(output)<-c("alpha (2-sided)","Power", "relevant relative effect p", "N (total sample size needed)", "t=n1/N", "n1 in Group 1", "n2 in Group 2", "N rounded", "n1 rounded", "n2 rounded")
colnames(output)<-""
return(output)}
