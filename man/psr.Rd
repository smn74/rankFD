% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/psr.R
\name{psr}
\alias{psr}
\title{A function for computing pseudo-ranks of data}
\usage{
psr(formula, data, psranks = "pseudorank")
}
\arguments{
\item{formula}{A model \code{\link{formula}} object. The left hand side
contains the response variable and the right hand side contains the factor
 variables of interest. Please use one-way layouts  for the computation of the pseudo-ranks only. In case of higher-way layouts, please use a 'help factor' that#' shrinks the layout to a one-way design.}

\item{data}{A data.frame, list or environment containing the variables in 
   
\code{formula}. The default option is \code{NULL}.}

\item{psranks}{A header specifying the name of the pseudo ranks in the output data set.}
}
\description{
The \code{psr()} function calculates pseuo-ranks of data in general factorial designs. It returns the input data set complemented by an additional variable containing the pseudo-ranks. We note that more efficient algorithms for the computation of pseudo-ranks are implemented within the package pseudorank.
}
\details{
The pseudo-ranks are exported within a new column attached to the given data set.
}
\examples{
data(Muco)
Muco2 <- psr(HalfTime~Disease,data=Muco, psranks="Mypseudos")

}
\references{
Konietschke, F., Hothorn, L. A., & Brunner, E. (2012). 
Rank-based multiple test procedures and simultaneous confidence intervals. Electronic Journal of Statistics, 6, 738-759.

Brunner, E., Bathke, A. C., Konietschke, F. (2018). Rank and pseudo-rank procedures for independent observations in factorial designs. Springer International
Publishing.

Happ, M., Zimmermann, G., Brunner, E., Bathke, A. C. (2020). Pseudo-ranks: How to calculate them efficiently in R. Journal of Statistical Software, 95(1), 1-22.
}
\seealso{
\code{\link{rankFD}}
}
