% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/samplesize.R
\name{WMWSSP}
\alias{WMWSSP}
\title{Sample size computation methods for the Wilcoxon-Mann-Whitney test based on pilot data}
\usage{
WMWSSP(x1, x2, alpha = 0.05, power = 0.8, t = 1/2)
}
\arguments{
\item{x1}{advance information for the first group}

\item{x2}{advance information for the second group}

\item{alpha}{two sided type I error rate}

\item{power}{power with the sample sizes of each group}

\item{t}{proportion of subjects in the first group}
}
\value{
Returns a data frame
}
\description{
The function implements the sample size formula proposed by Happ et al. (see reference below). It estimates the sample size needed to detect the effect with
pre-defined power at significance level alpha based on pilot data.
}
\examples{
x1 = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2)
x2 = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3)
WMWSSP(x1,x2,0.05,0.8,0.5)
}
\references{
Brunner, E., Bathke A. C. and Konietschke, F. Rank- and Pseudo-Rank Procedures in Factorial Designs - Using R and SAS. Springer Verlag. 
Happ, M., Bathke, A. C., & Brunner, E. (2019). Optimal sample size planning for the Wilcoxon-Mann-Whitney test. Statistics in medicine, 38(3), 363-375.
}
