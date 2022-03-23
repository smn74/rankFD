#' Steel-type multiple contrast tests
#' 
#' The function implements purely nonparametric Steel-type multiple contrast tests for either making
#' many-to-one (Dunnett-type) or all pairwise (Tukey-type) comparisons.
#' Null hypotheses are formulated in terms of the distribution functions. 
#' 
#' @param formula       A model \code{\link{formula}} object. The left hand side
#'                      contains the response variable and the right hand side contains 
#'                      the factor variable of interest. 
#' @param  data         A data.frame, list or environment containing the variables in 
#'                      \code{formula}. The default option is \code{NULL}.
#' @param control       Specification of the control group for making many-to-one-comparisons. If NULL, all-pairwise comparisons are performed.
#' @param alternative   Specification of the direction of the alternative. Default is two-sided.
#' @param info          Logical. If TRUE, additional output information and explanation is printed to the console.
#' @param correlation   Logical. If TRUE, the correlation matrix is printed.

#'
#' @details 
#' The steel() function calculates the Steel-type tests as explained by Munzel, U., Hothorn, L. A. (2001). A unified approach to simultaneous rank test procedures
#' in the unbalanced one-way layout. Biometrical Journal: Journal of Mathematical Methods in Biosciences, 43(5), 553-569.
#
#'   
#' @return A list containing the following components:
#' \item{Data.Info}{Groups and sample sizes of the data}           
#' \item{Analysis}{Data frame containing the test results (comparison, relative effect estimator, standard error, test statistic and p-value.)}
#' \item{Correlation}{Estimated correlation matrix}

#'
#'
#' @examples
#' 
#' data(Muco)
#' model.oneway <- steel(HalfTime ~ Disease, data = Muco,info=TRUE,correlation=TRUE)
#' 
#' 
#' @references 
#' Brunner, E., Bathke, A.C., Konietschke, F. Rank and Pseudo-Rank Procedures 
#' for Independent Observations in Factorial Designs. Springer International Publishing, 2018.
#' 
#'Munzel, U., Hothorn, L. A. (2001). A unified approach to simultaneous rank test procedures in the unbalanced one-way layout. Biometrical Journal: Journal of
#' Mathematical Methods in Biosciences, 43(5), 553-569.
#' 
#'Konietschke, F., Hothorn, L. A., Brunner, E. (2012). Rank-based multiple test procedures and simultaneous confidence intervals. Electronic Journal of Statistics,
#' 6, 738-759.
#' 
#'
#' @importFrom lattice xyplot panel.superpose panel.arrows panel.points panel.xyplot panel.abline
#' @importFrom stats formula model.frame pchisq pf qnorm terms var aggregate as.formula confint cov median pnorm pt qt quantile sd cov2cor
#' @importFrom utils read.table
#' @importFrom MASS ginv 
#' @importFrom mvtnorm qmvnorm pmvt qmvt
#' 
#' @export
steel<-function(formula,data, control = NULL,alternative=c("two.sided","less","greater"),info=TRUE,correlation=TRUE)
{
  alternative <- match.arg(alternative)
  ssq <- function(x) sum(x * x)
 #-----------------------Determine the model-----------------------------------#
  if (length(formula) != 3) {
  stop("You can only analyse one-way layouts!")}
  dat <- model.frame(formula, droplevels(data))
  if (ncol(dat) != 2) {
  stop("Specify one response and only one class variable in the formula")}
    if (is.numeric(dat[, 1]) == FALSE) {stop("Response variable must be numeric")}
    res <- dat[, 1]
    grp <- factor(dat[, 2])
    fl <- levels(grp)
    a <- nlevels(grp)
    samples <- split(res, grp)
    if (a <= 2) {stop("You want to perform a two-sample test. Please use the function rank.two.samples")}
    n <- sapply(samples, length)
    if (is.null(control)) {
    tmp <- expand.grid(1:a, 1:a)
    ind <- tmp[[1]] > tmp[[2]]
    vi <- tmp[[2]][ind]
    vj <- tmp[[1]][ind]}
    else {if (!any(fl == control)) {
            msg <- paste("Wrong control-group specification\n",
                "The data does not contain a group with factor-level ",
                control, sep = "")
            stop(msg, FALSE)}
        cg <- which(fl == control)
        vi <- which((1:a) != cg)
        vj <- rep(cg, a - 1)}
    nc <- length(vi)
    cmpid <- paste(vi, "-", vj, sep = "")
    gn <- n[vi] + n[vj]
    cmpid <- paste("p(",fl[vi],",", fl[vj],")")
    intRanks <- lapply(samples, rank)
    pairRanks <- lapply(1:nc, function(arg) {rank(c(samples[[vi[arg]]], samples[[vj[arg]]]))})
    pd <- sapply(1:nc, function(arg) {
        i <- vi[arg]
        j <- vj[arg]
        (sum(pairRanks[[arg]][(n[i] + 1):gn[arg]])/n[j] - (n[j] +1)/2)/n[i]})
        dij <- dji <- list(0)
        lambda <- sqrt(n[vj]/(gn + 1))
        vd.st <- sapply(1:nc, function(arg) ssq(pairRanks[[arg]] -
        (gn[arg] + 1)/2))/(n[vi] * n[vj] * (gn - 1))
  std.st <- sqrt(vd.st/gn)
  t.st <- (pd - 0.5) * sqrt(gn/vd.st)
   rho.bf <- rho.st <- diag(nc)
   for (x in 1:(nc - 1)) {
        for (y in (x + 1):nc) {
            i <- vi[x]
            j <- vj[x]
            v <- vi[y]
            w <- vj[y]
            p <- c(i == v, j == w, i == w, j == v)
            if (sum(p) == 1) {
                cl <- list(function() sqrt(n[j]/(n[i]+n[j]+1))*sqrt(n[w]/(n[i]+n[w]+1))
                 , function() sqrt(n[i]/(n[i]+n[j]+1))*sqrt(n[v]/(n[v]+n[j]+1)) ,
                  function() -(sqrt(n[j]/(n[i]+n[j]+1))*sqrt(n[v]/(n[v]+n[i]+1))),
                  function() -(sqrt(n[w]/(n[w]+n[j]+1))*sqrt(n[i]/(n[j]+n[i]+1))))
                case <- (1:4)[p]
                rho.st[x, y] <- rho.st[y, x] <- 1 * cl[[case]]()

            }
        }
    }

    switch(alternative,two.sided={
    text.alternative = paste("Relative effect is unequal to 1/2 (distributions differ)")
    p.st<-sapply(1:nc,function(arg){
    1-pmvnorm(lower=-abs(t.st[arg]), abs(t.st[arg]),corr=rho.st,mean=rep(0,nc))})},
     less={text.alternative=paste("Relative effect p is less than 1/2 (distribution 2 is smaller)")
       p.st<-sapply(1:nc,function(arg){
    pmvnorm(lower=-Inf, t.st[arg],corr=rho.st,mean=rep(0,nc))})},
    greater={text.alternative=paste("Relative effect p is greater than 1/2 (distribution 2 is larger)")
       p.st<-sapply(1:nc,function(arg){
    pmvnorm(lower=t.st[arg],Inf,corr=rho.st,mean=rep(0,nc))})})
    dataStructure <- data.frame("group index" = 1:a, "grp level" = fl,nobs = n)
    data.info <- data.frame(row.names = 1:a, Sample = fl, Size = n)
    test.st <- data.frame(Comparison = cmpid, Rel.Effect = pd, Variance=vd.st, Std.Error = std.st, T = t.st, p.value = p.st)
    rownames(test.st)<-1:nc
 if(correlation==TRUE){R=rho.st}
 if(correlation==FALSE){R=NULL}
 result <- list(Data.Info = data.info, Analysis = test.st, Correlation=R)
    if (info == TRUE) {
    cat("\n", "#--------------Nonparametric Multiple Steel Type Tests-----------------------#", 
    "\n", "\n", "-", "Alternative Hypothesis:  ",text.alternative, 
    "\n", "\n", "#---------------------------Interpretation----------------------------------#", 
    "\n", "p(a,b)", ">", "1/2", 
    ":", "b tends to be larger than a", "\n", 
    "#---------------------------------------------------------------------------#", 
    "\n", "\n")}
    #class(result)<-"rankFD"
    return(result)}













