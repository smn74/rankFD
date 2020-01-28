
#' A function for computing pseudo-ranks of data
#' 
#' The \code{psr()} function calculates the pseuo-ranks of data in general factorial designs.

#' 
#' @param formula A model \code{\link{formula}} object. The left hand side
#'    contains the response variable and the right hand side contains the factor
#'    variables of interest. Please use one-way layouts only for the computation of the pseudo-ranks.
#' @param data A data.frame, list or environment containing the variables in 
#'    \code{formula}. The default option is \code{NULL}.
#' @param psranks A header specifying the name of the pseudo ranks in the output data set. 

#' @references Konietschke, F., Hothorn, L. A., & Brunner, E. (2012). 
#'Rank-based multiple test procedures and simultaneous confidence intervals. Electronic Journal of Statistics, 6, 738-759.
#' 
#' Kaufmann, J., Werner, C., and Brunner, E. (2005). Nonparametric methods for analysing the
#' accuracy of diagnostic tests with multiple readers. Statistical Methods in Medical Research 14, 129 - 146
#' 
#' 
#' @details The pseudo-ranks are exported within a new column attached to the given data set. 
#' @examples
#' data(Muco)
#' Muco2 <- psr(HalfTime~Disease,data=Muco, psranks="Mypseudos")
#' 
#' @seealso \code{\link{rankFD}}
#' 
#' @export
psr=function(formula, data, psranks){
 dat.Model0 <- model.frame(formula, data)
dat.Model0[,2] = as.factor(dat.Model0[,2])
dat.Model0 <- dat.Model0[order(dat.Model0[, 2]), ]
n <- aggregate(formula, data = dat.Model0, length)[,2]
gn=sum(n)
a=length(n)
n1     = n[1]                                                            
vsum   = matrix(1/n1,nrow=n1,ncol=1)                                                 
for(i in 2:a){                                                             
   ni = n[i]                                                             
 vsum = rbind(vsum, matrix(1/ni, nrow=ni,ncol=1) )                                    
}

X=dat.Model0[,1]
r=X
gn=length(X)
pr=c()
for ( k in 1:gn){                                                            
   pr[k] = ((sign(r[k]*matrix(1,nrow=1,ncol=gn) - t(r)) + matrix(1,nrow=1,ncol=gn))*0.5)%*%vsum*gn/a + 0.5 
} 

dat.Model0$psr=pr
names(dat.Model0)[3] = psranks
dat.Model0
}

