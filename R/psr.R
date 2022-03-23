
#' A function for computing pseudo-ranks of data
#' 
#' The \code{psr()} function calculates pseuo-ranks of data in general factorial designs. It returns the input data set complemented by an additional variable containing the pseudo-ranks. We note that more efficient algorithms for the computation of pseudo-ranks are implemented within the package pseudorank. 

#' 
#' @param formula A model \code{\link{formula}} object. The left hand side
#'   contains the response variable and the right hand side contains the factor
#'    variables of interest. Please use one-way layouts  for the computation of the pseudo-ranks only. In case of higher-way layouts, please use a 'help factor' that#' shrinks the layout to a one-way design. 
#' @param data A data.frame, list or environment containing the variables in 
#'    
#'\code{formula}. The default option is \code{NULL}.
#' @param psranks A header specifying the name of the pseudo ranks in the output data set. 

#' @references 
#'Konietschke, F., Hothorn, L. A., & Brunner, E. (2012). 
#'Rank-based multiple test procedures and simultaneous confidence intervals. Electronic Journal of Statistics, 6, 738-759.
#' 
#'Brunner, E., Bathke, A. C., Konietschke, F. (2018). Rank and pseudo-rank procedures for independent observations in factorial designs. Springer International
#' Publishing.
#'
#'Happ, M., Zimmermann, G., Brunner, E., Bathke, A. C. (2020). Pseudo-ranks: How to calculate them efficiently in R. Journal of Statistical Software, 95(1), 1-22.
#' 
#' @details The pseudo-ranks are exported within a new column attached to the given data set. 
#' @examples
#' data(Muco)
#' Muco2 <- psr(HalfTime~Disease,data=Muco, psranks="Mypseudos")
#' 
#' @seealso \code{\link{rankFD}}
#' 
#' @export
psr=function(formula, data,psranks="pseudorank"){
pr=c()
dat.Model0 <- model.frame(formula, data)
X=dat.Model0[,1]
dat.Model0[,2] = as.factor(dat.Model0[,2])
dat.Model0 <- dat.Model0[order(dat.Model0[, 2]), ]
n <- aggregate(formula, data = dat.Model0, length)[,2]
gn=sum(n)
a=length(n)
vsum =matrix(rep(1/n, n), ncol = 1) 
for (k in 1:gn){                                                            
pr[k] = ((sign(X[k]*matrix(1,nrow=1,ncol=gn) - t(X)) + matrix(1,nrow=1,ncol=gn))*0.5)%*%vsum*gn/a + 0.5 } 
dat.Model0$pseudorank=pr
if(!is.null(psranks)){
names(dat.Model0)[3] = psranks}
return(dat.Model0)
}



