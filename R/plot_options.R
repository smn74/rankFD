#################################################################
#Grafik Option
#################################################################
plotting <- function(nf, fac_names, n.hypotheses,Descriptive.Factors, CI.method, alpha, cex=1.3, cex.lab=1.3,cex.axis=1.3,col=1,pch=19,lwd=4,cex.ci=4,main=NULL,xlab=NULL,ylab=NULL, ...){

  if(nf ==1){Faktor = fac_names[1]}
  if(nf > 1){print("Please choose the factor you wish to plot (for interaction type something like group1:group2) and confirm by pressing 'Enter'")
    Faktor <- scan("", what = "character")}
  Fak.split <- strsplit(Faktor,":")[[1]]
  l.Fak.split <- length(Fak.split)
  if (!(Faktor %in% fac_names)) {stop("Please enter a valid factor name!")}
  for(i in 1:n.hypotheses){
    
    if(names(Descriptive.Factors)[i] == Faktor){
      posP <- which(names(Descriptive.Factors)[i] == Faktor)
      DatenPlot <- data.frame(Descriptive.Factors[[i]])
       lower=DatenPlot$Lower
       upper=DatenPlot$Upper
       DatenPlot$pd<-DatenPlot$Rel.Effect
      if(is.null(xlab)){
      text.X<-paste(names(DatenPlot[1]))}
       if(!is.null(xlab)){text.X<-xlab}
      if(is.null(ylab)){
        text.Ci <- paste((1 - alpha) * 100, "%", "Confidence Intervals for Relative Effects")}
        if(!is.null(ylab)){text.Ci<-ylab}
      if (l.Fak.split==1){
        print(xyplot(pd ~ DatenPlot[,1], group=DatenPlot[,1],data = DatenPlot, 
                     type = 'p',  main=main,
                     xlab=list(text.X,cex=cex.lab),
                     ylab=list(text.Ci,cex=cex.lab),
                      scales=list(x=list(cex=cex.axis),y=list(cex=cex.axis)),
                     col = col, pch = pch, cex=cex, lwd=lwd,  ylim = c(0, 1), 
                     upper = upper,
                     lower = lower,
                      panel = function(x, y, ...){
                     panel.superpose(x, y,panel.groups = function(x, y, upper, lower,subscripts, ..., font, fontface) {
                       upper <- upper[subscripts]
                       lower <- lower[subscripts]
                       panel.arrows(x, lower, x, upper,code=4,lwd=lwd,col=col)   
                       panel.points(x, lower,pch="_",cex=cex.ci,lwd=lwd,col=col)
                       panel.points(x, upper,pch="_",cex=cex.ci,lwd=lwd,col=col)
                     } , ...)
                       panel.xyplot(x, y, ...)}
                    ))
                    
                    }
      
      if (l.Fak.split==2){
        print(xyplot(pd ~ DatenPlot[,1]|DatenPlot[,2], 
                     group=DatenPlot[,1],data = DatenPlot, type = 'p', main=main,
                     xlab=list(text.X,cex=cex.lab),
                     ylab=list(text.Ci,cex=cex.lab),
                     scales=list(x=list(cex=cex.axis),y=list(cex=cex.axis)),
                      col = col, pch = pch, cex=cex, lwd=lwd,  ylim = c(0, 1), 
                     upper = upper,
                     lower = lower,
                     par.strip.text =list(cex=cex),
                     panel = function(x, y, ...){
                       panel.superpose(x, y,panel.groups = function(x, y, upper, lower,subscripts, ..., font, fontface) {
                         upper <- upper[subscripts]
                         lower <- lower[subscripts]
                         panel.arrows(x, lower, x, upper,code=4,lwd=lwd,col=col)   
                       panel.points(x, lower,pch="_",cex=cex.ci,lwd=lwd,col=col)
                       panel.points(x, upper,pch="_",cex=cex.ci,lwd=lwd,col=col)
                       }, ...)
                       panel.xyplot(x, y, ...)
                     }))}
      
      if (l.Fak.split==3){
        print(xyplot(pd ~ DatenPlot[,1]|DatenPlot[,2]*DatenPlot[,3], 
                     group=DatenPlot[,1],data = DatenPlot, type = 'p', main=main,
                    xlab=list(text.X,cex=cex.lab),
                     ylab=list(text.Ci,cex=cex.lab),
                     scales=list(x=list(cex=cex.axis),y=list(cex=cex.axis)),
                      col = col, pch = pch, cex=cex, lwd=lwd,  ylim = c(0, 1), 
                     upper = upper,
                     lower = lower,
                     par.strip.text =list(cex=cex),
                     panel = function(x, y, ...){
                       panel.superpose(x, y,panel.groups = function(x, y, upper, lower,subscripts, ..., font, fontface) {
                         upper <- upper[subscripts]
                         lower <- lower[subscripts]
                          panel.arrows(x, lower, x, upper,code=4,lwd=lwd,col=col)   
                          panel.points(x, lower,pch="_",cex=cex.ci,lwd=lwd,col=col)
                          panel.points(x, upper,pch="_",cex=cex.ci,lwd=lwd,col=col)
                       }, ...)
                       panel.xyplot(x, y, ...)
                     })) }
      
      if (l.Fak.split>=4){
        stop("4 and higher way interactions cannot be plotted!")
      }
    }}
}

plotting.sci <- function(resmctp, hypothesis=hypothesis, sci.method=sci.method,alpha=alpha,cex=1.3,cex.lab=1.3,cex.axis=1.3,col=1,pch=19,lwd=4,cex.ci=4,main=NULL, xlab=NULL,ylab=NULL,...){
if(is.null(resmctp)){
  stop("No contrast has been specified in the rankFD call.")
}
  if(hypothesis=="H0F"){stop("Simultaneous Confidence Intervals can only be computed when hypothesis='H0p'. Please change the argument.")}
if(hypothesis=="H0p"){
DatenPlot <- data.frame(resmctp[[2]])
nc <- nrow(DatenPlot)
comp <- NULL
DatenPlot$comp <- 1:nc
lower=DatenPlot[,4]
upper=DatenPlot[,5]

if(is.null(xlab)){
      text.X<-paste(names(DatenPlot[1]))}
       if(!is.null(xlab)){text.X<-xlab}
      if(is.null(ylab)){
        text.Ci <- paste((1 - alpha) * 100, "%", "Simultaneous Confidence Intervals")}
        if(!is.null(ylab)){text.Ci<-ylab}
   print(xyplot(DatenPlot[,1] ~ comp, group=comp,data = DatenPlot, 
                     type = 'p',  main=main,
                      xlab=list(text.X,cex=cex.lab),
                     ylab=list(text.Ci,cex=cex.lab),
                      scales=list(x=list(at=1:nc, labels=c(paste("C",1:nc,sep="")),cex=cex.axis),y=list(cex=cex.axis)),
                     col = col, pch = pch, cex=cex, lwd=lwd,  ylim = c(-1, 1), 
                     upper = upper,
                     lower = lower,
                      panel = function(x, y, ...){
                     panel.superpose(x, y,panel.groups = function(x, y, upper, lower,subscripts, ..., font, fontface) {
                       upper <- upper[subscripts]
                       lower <- lower[subscripts]
                       panel.arrows(x, lower, x, upper,code=4,lwd=lwd,col=col)   
                       panel.points(x, lower,pch="_",cex=cex.ci,lwd=lwd,col=col)
                       panel.points(x, upper,pch="_",cex=cex.ci,lwd=lwd,col=col)
                       panel.abline(h=0, col="red", lty=1, lwd=2)
                     } , ...)
                       panel.xyplot(x, y, ...)}
                    ))} }

plot2samples <- function(method, Effects, alpha,cex.lab=1.3, cex=1.3,cex.axis=1.3,col=1,pch=19,lwd=4,cex.ci=4,main=NULL,xlab=NULL,ylab=NULL, ...){
         
        DatenPlot <-Effects
        if(is.null(xlab)){
        text.Ci <- paste((1 - alpha) * 100, "%", "Confidence Interval for Relative Effect")}
        if(!is.null(xlab)){text.Ci<-xlab}
        if(is.null(ylab)){text.Y <- ""}
        if(!is.null(ylab)){text.Y <- ylab}
        Lowerp <- "|"
        pd <- DatenPlot[,2]
        Lower <-DatenPlot[,5]
        Upper <-DatenPlot[,6]
        plot(rep(pd, 1), 1:1, xlim = c(0, 1), pch = pch, 
            axes = FALSE, xlab = text.Ci,cex.lab=cex.lab, ylab = text.Y,cex=cex,main=main)
        points(Lower, 1:1, pch = Lowerp, font = 2, cex = cex.ci)
        points(Upper, 1:1, pch = Lowerp, font = 2, cex = cex.ci)
        abline(v = 0.5, lty = 3, lwd = 2,col="red")
        polygon(x = c(Lower, Upper), y = c(1, 1), lwd = lwd)   
        axis(1, at = seq(0, 1, 0.1),cex.axis=cex.axis)
        box()
       
   }
   



