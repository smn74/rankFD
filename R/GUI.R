#' A graphical user interface for the package rankFD
#' 
#' This function provides a graphical user interface for calculating rank-based 
#' statistical tests in general factorial designs.
#' 
#' The function produces a GUI for the calculation of the test statistics and 
#' for plotting. Data can be loaded via the "load data" button. The formula 
#' and the significance level alpha (default: 0.05) need to be specified.
#' One can choose between two different null hypotheses (\eqn{H_0^F} and \eqn{H_0^p}, see the details
#' to \code{\link{rankFD}}) to be tested as well as
#' weighted or unweighted effects as discussed in Brunner et al. (2016) (\eqn{r_i} and \eqn{p_i} in
#' their notation).
#' If the plot option is chosen, an additional window opens containing
#' information on the plots.
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
#' @export

# GUI
rankFD_GUI <- function() {
  ## Run on "Load"
  requireNamespace("RGtk2", quietly=TRUE)
  if(!("package:RGtk2" %in% search())){attachNamespace("RGtk2")}
  getDirectory <- function(button, user.data){
    directory = file.choose()
    RGtk2::gtkEntrySetText(filename,directory)
  }  
  ## Run on "OK"
  performStatistics <- function(button, user.data) {
    res <- NULL
    d <- NULL
    error <- NULL
    warning <- NULL
    # Get the information about data and the file
    the.file <- filename$getText()
    the.formula <- formula(filename1$getText())
    the.alpha <- as.numeric(filename3$getText())
    the.rounds <- as.numeric(filenameround$getText())
    the.plot <- toPlot$active
    the.output <-toOutput$active
    the.effect <-c("weighted","unweighted")[comboboxeffect$active+1]
    the.hypothesis <-c("H0F","H0p")[comboboxhypothesis$active+1]
    the.CImethod <-c("logit","normal")[comboboxCI.method$active+1]
    the.scimethod <- c("fisher","multi.t")[comboboxscimethod.list$active+1]
    the.SCIplot <- toSCIPlot$active
    the.contrast <-filenamecontrast$getText()
    the.contrasttype <-c("","Dunnett","Tukey", "Sequen", "AVE", "Changepoint", 
                    "Williams", "Marcus", "McDermott", "UmbrellaWilliams", "GrandMean")[comboboxcontrast.list$active+1]
    the.sep <- sepEntry$getText()
    the.headers <- headersEntry$active
    the.dec <- decEntry$getText()
    the.factors = toFactors$active
    d <- read.table(the.file, sep = the.sep, header = the.headers,
                    dec = the.dec)
    
    dat.ModelTEST <- model.frame(the.formula,d)
    #--------Numbers of factors---------#
    
    for(i in length(dat.ModelTEST):2) dat.ModelTEST[,i] <-as.factor(dat.ModelTEST[,i])
    nfTEST = ncol(dat.ModelTEST)-1
    
    n.levelsTEST = c()
    for (i in 2:(nfTEST+1)){
      n.levelsTEST[i-1]= nlevels(dat.ModelTEST[,i])
    }
    if(!(nfTEST==1 && n.levelsTEST ==2)){
    if(nchar(the.contrast)==0){
    res <- rangtests(the.formula, d,plot_CI = the.plot,alpha=the.alpha,effect=the.effect,
                       hypothesis=the.hypothesis, Factor.Information=the.factors,CI.method=the.CImethod, info=the.output,rounds=the.rounds)
      print(res)                 
    }
    if(nchar(the.contrast)!=0 && nchar(the.contrasttype)!=0){
      res <- rangtests(the.formula, d,plot_CI = the.plot,alpha=the.alpha,effect=the.effect,
      hypothesis=the.hypothesis, Factor.Information=the.factors,CI.method=the.CImethod,contrast=list(the.contrast,the.contrasttype),sci.method=the.scimethod,
      info=the.output, plot_SCI= the.SCIplot,rounds=the.rounds)
       print(res) 
       
      }
      
      if(nchar(the.contrast)!=0 && nchar(the.contrasttype)==0){
      res <- rangtests(the.formula, d,plot_CI = the.plot,alpha=the.alpha,effect=the.effect, plot_SCI= the.SCIplot,
      hypothesis=the.hypothesis, Factor.Information=the.factors,CI.method=the.CImethod,contrast=list(the.contrast),sci.method=the.scimethod,info=the.output,rounds=the.rounds)
       print(res)  
      }
     
    }
    
    if(nfTEST == 1 && n.levelsTEST == 2){
      ###########################################################
      
      #-------------GUI Function!!--------#
      calculateGUItwosamples <- function() {
        twosamples <- function(button, user.data) {
          
          error <- NULL
                    
          if (!is.null(error)) {
            hbox <- RGtk2::gtkHBoxNew()
            vbox$packStart(hbox,FALSE,FALSE,0)
            label <- RGtk2::gtkLabel(error)
            hbox$packStart(label,FALSE,FALSE,0)
          }
          
          
          the.alternative <-c("two.sided","less","greater")[comboboxalternative$active+1]
          the.nperm <-as.numeric(filename14$getText())
          the.alpha <- as.numeric(filename16$getText())
          the.rounds <- as.numeric(filenameround$getText())
          the.shift <- toPlot$active
          the.Permutation<-toPermutation$active
          the.output <-toOutput$active
          the.wilcoxon <-c("asymptotic","exact")[comboboxwilcoxon$active+1]
          the.method = c("normal","t.app","logit","probit")[combobox$active+1]
          res<-rank.two.samples(the.formula, d, conf.level = 1-the.alpha, #plot.simci=the.plot,
                                alternative = the.alternative, method =the.method, permu=the.Permutation,  info = the.output, 
                                shift.int=the.shift,nperm = the.nperm,wilcoxon=the.wilcoxon,rounds=the.rounds) 
                
          print(res)
          ## Run on "OK" 
        }
        
        ##################################################
        # Create window
        window <- RGtk2::gtkWindow()
        # Add title
        window["title"] <- "Two-sample Rank Tests"
        
              # Add a frame
      frame <- RGtk2::gtkFrameNew("")
      window$add(frame)
        
        # Create vertical container for file name entry
        vbox <- RGtk2::gtkVBoxNew(FALSE, 8)
        vbox$setBorderWidth(24)
        frame$add(vbox)
        # Add horizontal container for every widget line
        hbox <- RGtk2::gtkHBoxNew(FALSE, 8)
        vbox$packStart(hbox, FALSE, FALSE, 0)
        
        labeltest <- RGtk2::gtkLabelNewWithMnemonic("#------------------------------------------Tests and Confidence Intervals for Relative Effects--------------------------------#")
        hbox$packStart(labeltest,FALSE,FALSE,0)    
        
        ############################################################
        
        # Add an horizontal container to specify parameters
        
        hbox <- RGtk2::gtkHBoxNew(FALSE,8)
        vbox$packStart(hbox, FALSE, FALSE, 0)
        labeltest <- RGtk2::gtkLabelNewWithMnemonic("Method")
        hbox$packStart(labeltest,FALSE,FALSE,0)
        model<-RGtk2::rGtkDataFrame(c("normal","t.app","logit","probit"))
        combobox <- RGtk2::gtkComboBox(model)
        #combobox allowing to decide whether we want result as integer or double
        crt <- RGtk2::gtkCellRendererText()
        combobox$packStart(crt)
        combobox$addAttribute(crt, "text", 0)
        
        RGtk2::gtkComboBoxSetActive(combobox,0)
        hbox$packStart(combobox)
        
        labelAlternative <- RGtk2::gtkLabelNewWithMnemonic("Alternative")
        hbox$packStart(labelAlternative,FALSE,FALSE,0)
        alternative<-RGtk2::rGtkDataFrame(c("two.sided","less","greater"))
        comboboxalternative <- RGtk2::gtkComboBox(alternative)
        #combobox allowing to decide whether we want result as integer or double
        crtalternative <- RGtk2::gtkCellRendererText()
        comboboxalternative$packStart(crtalternative)
        comboboxalternative$addAttribute(crtalternative, "text", 0)
        RGtk2::gtkComboBoxSetActive(comboboxalternative,0)
        hbox$packStart(comboboxalternative)
        
           # Add Shift Effect Option
        #hbox <- RGtk2::gtkHBoxNew(FALSE,8)
        #vbox$packStart(hbox, FALSE, FALSE, 0)
        label <- RGtk2::gtkLabelNewWithMnemonic("Permutation Test")
        hbox$packStart(label,FALSE,FALSE,0)
        toPermutation <- RGtk2::gtkCheckButton()
        hbox$packStart(toPermutation,FALSE,FALSE,0)
        
        label14 <- RGtk2::gtkLabelNewWithMnemonic("nperm")
        hbox$packStart(label14,FALSE,FALSE,0)
        # Add entry in the second column; named "filename14"
        filename14 <- RGtk2::gtkEntryNew()
        filename14$setWidthChars(10)
        filename14$setText(10000)
        label14$setMnemonicWidget(filename14)
        hbox$packStart(filename14,FALSE,FALSE,0)
        
          label16 <- RGtk2::gtkLabelNewWithMnemonic("Alpha")
        hbox$packStart(label16,FALSE,FALSE,0)
        # Add entry in the second column; named "filename16"
        filename16 <- RGtk2::gtkEntryNew()
        filename16$setWidthChars(10)
        filename16$setText(0.05)
        label16$setMnemonicWidget(filename16)
        hbox$packStart(filename16,FALSE,FALSE,0)
        
        #hbox <- RGtk2::gtkHBoxNew(FALSE,8)
        #vbox$packStart(hbox, FALSE, FALSE, 0)
        
        # Add separator
        vbox$packStart(RGtk2::gtkHSeparatorNew(), FALSE, FALSE, 0)
        hbox <- RGtk2::gtkHBoxNew(FALSE,8)
        vbox$packStart(hbox, FALSE, FALSE, 0)
        labeltest <- RGtk2::gtkLabelNewWithMnemonic("#-----------------------------------------------Wilcoxon-Mann-Whitney Test---------------------------------------------#")
        hbox$packStart(labeltest,FALSE,FALSE,0) 
        vbox$packStart(RGtk2::gtkHSeparatorNew(), FALSE, FALSE, 0)
        hbox <- RGtk2::gtkHBoxNew(FALSE,8)
        vbox$packStart(hbox, FALSE, FALSE, 0)
        labelWilcoxon <- RGtk2::gtkLabelNewWithMnemonic("Wilcoxon Test")
        hbox$packStart(labelWilcoxon,FALSE,FALSE,0)
        wilcoxon<-RGtk2::rGtkDataFrame(c("asymptotic","exact"))
        comboboxwilcoxon <- RGtk2::gtkComboBox(wilcoxon)
        crtwilcoxon <- RGtk2::gtkCellRendererText()
        comboboxwilcoxon$packStart(crtwilcoxon)
        comboboxwilcoxon$addAttribute(crtwilcoxon, "text", 0)
        
        RGtk2::gtkComboBoxSetActive(comboboxwilcoxon,0)
        hbox$packStart(comboboxwilcoxon)  
        
         RGtk2::gtkComboBoxSetActive(comboboxwilcoxon,0)
        hbox$packStart(comboboxwilcoxon)        
        
      
        
        # Add Shift Effect Option
        #hbox <- RGtk2::gtkHBoxNew(FALSE,8)
        #vbox$packStart(hbox, FALSE, FALSE, 0)
        label <- RGtk2::gtkLabelNewWithMnemonic("Shift Effects?")
        hbox$packStart(label,FALSE,FALSE,0)
        toPlot <- RGtk2::gtkCheckButton()
        hbox$packStart(toPlot,FALSE,FALSE,0)
        
        # Add Output Information Option
        #hbox <- RGtk2::gtkHBoxNew(FALSE,8)
        #vbox$packStart(hbox, FALSE, FALSE, 0)
        label <- RGtk2::gtkLabelNewWithMnemonic("Output Information?")
        hbox$packStart(label,FALSE,FALSE,0)
        toOutput <- RGtk2::gtkCheckButton()
        hbox$packStart(toOutput,FALSE,FALSE,0)
        
        labelround <- RGtk2::gtkLabelNewWithMnemonic("_Output Decimals")
        hbox$packStart(labelround,FALSE,FALSE,0)
        filenameround <- RGtk2::gtkEntryNew()
        filenameround$setWidthChars(10)
        filenameround$setText(4)
        labelround$setMnemonicWidget(filenameround)
        hbox$packStart(filenameround,FALSE,FALSE,0)
        
        ############################################################
        
        # Add button
        the.buttons <- RGtk2::gtkHButtonBoxNew()
        the.buttons$setBorderWidth(5)
        vbox$add(the.buttons)
        the.buttons$setLayout("spread")
        the.buttons$setSpacing(40)
        buttonOK <- RGtk2::gtkButtonNewFromStock("gtk-ok")
        RGtk2::gSignalConnect(buttonOK, "clicked", twosamples)
        the.buttons$packStart(buttonOK,fill=F)
        
        ##################################################
        
        buttonCancel <- RGtk2::gtkButtonNewFromStock("gtk-close")
        RGtk2::gSignalConnect(buttonCancel, "clicked", window$destroy)
        the.buttons$packStart(buttonCancel,fill=F)
        
      }
      ###################################
      
      #end of two samples!!
      calculateGUItwosamples()
            
    }
  }
  
  # Create window
  window <- RGtk2::gtkWindow()
  # Add title
  window["title"] <- "Rank Methods for Factorial Designs"
  
  # Add a frame
  frame <- RGtk2::gtkFrameNew("") #add title here
  window$add(frame)
  
  # Create vertical container for file name entry
  vbox <- RGtk2::gtkVBoxNew(FALSE, 8)
  vbox$setBorderWidth(24)
  frame$add(vbox)
  # Add horizontal container for every widget line
  hbox <- RGtk2::gtkHBoxNew(FALSE, 8)
  vbox$packStart(hbox, FALSE, FALSE, 0)
  
  #-------------------Data Set Location and its Characteristics----------------#
  # Add an horizontal container to specify input file options
  # are headers included in the file?
 
    label <- RGtk2::gtkLabelNewWithMnemonic("_File name")
  hbox$packStart(label,FALSE,FALSE,0)
  # Add entry in the second column; named "filename"
  filename <- RGtk2::gtkEntryNew()
  filename$setWidthChars(50)
  label$setMnemonicWidget(filename)
  hbox$packStart(filename,FALSE,FALSE,0)
  
  #hbox <- RGtk2::gtkHBoxNew(FALSE,8)
  vbox$packStart(hbox, FALSE, FALSE, 0)
  label <- RGtk2::gtkLabelNewWithMnemonic("_Headers?")
  hbox$packStart(label,FALSE,FALSE,0)
  headersEntry <- RGtk2::gtkCheckButton()
  headersEntry$active <- TRUE
  hbox$packStart(headersEntry,FALSE,FALSE,0)
  label$setMnemonicWidget(headersEntry)
  
  # are headers included in the file?
  label <- RGtk2::gtkLabelNewWithMnemonic("Col. _Separator?")
  hbox$packStart(label,FALSE,FALSE,0)
  sepEntry <- RGtk2::gtkEntryNew()
  sepEntry$setWidthChars(1)
  sepEntry$setText("")
  hbox$packStart(sepEntry,FALSE,FALSE,0)
  label$setMnemonicWidget(sepEntry)
  
  # what's the character used for decimal points?
  label <- RGtk2::gtkLabelNewWithMnemonic("_Dec. character?")
  hbox$packStart(label,FALSE,FALSE,0)
  decEntry <- RGtk2::gtkEntryNew()
  decEntry$setWidthChars(1)
  decEntry$setText(".")
  hbox$packStart(decEntry,FALSE,FALSE,0)
  label$setMnemonicWidget(decEntry)
  

  
  #----------------------------------------------------------------------------#
      # Add separator
  hbox <- RGtk2::gtkHBoxNew(FALSE,8)
  vbox$packStart(hbox, FALSE, FALSE, 0)
  vbox$packStart(RGtk2::gtkHSeparatorNew(), FALSE, FALSE, 0)
  # Add label in first column
  label1 <- RGtk2::gtkLabelNewWithMnemonic("_Formula")
  hbox$packStart(label1,FALSE,FALSE,0)
  # Add entry in the second column; named "filename1"
  filename1 <- RGtk2::gtkEntryNew()
  filename1$setWidthChars(50)
  label1$setMnemonicWidget(filename1)
  hbox$packStart(filename1,FALSE,FALSE,0)
  
  
  
 
 #---------------------Hypothesis and Ranking Method---------------------------#
  
  
  labelhypothesis <- RGtk2::gtkLabelNewWithMnemonic("Hypothesis")
  hbox$packStart(labelhypothesis,FALSE,FALSE,0)
  hypothesis<-RGtk2::rGtkDataFrame(c("Distribution Functions H_0^F","Relative Effects H_0^p"))
  comboboxhypothesis <- RGtk2::gtkComboBox(hypothesis)
  crthypothesis <- RGtk2::gtkCellRendererText()
  comboboxhypothesis$packStart(crthypothesis)
  comboboxhypothesis$addAttribute(crthypothesis, "text", 0)
  RGtk2::gtkComboBoxSetActive(comboboxhypothesis,0)
  hbox$packStart(comboboxhypothesis)
  
  labeleffect <- RGtk2::gtkLabelNewWithMnemonic("Effects")
  hbox$packStart(labeleffect,FALSE,FALSE,0)
  effect<-RGtk2::rGtkDataFrame(c("weighted (n_i/N)","unweighted (1/d)"))
  comboboxeffect <- RGtk2::gtkComboBox(effect)
  crteffect <- RGtk2::gtkCellRendererText()
  comboboxeffect$packStart(crteffect)
  comboboxeffect$addAttribute(crteffect, "text", 0)
  RGtk2::gtkComboBoxSetActive(comboboxeffect,0)
  hbox$packStart(comboboxeffect)
  
  labelCI.method <- RGtk2::gtkLabelNewWithMnemonic("CI Method")
  hbox$packStart(labelCI.method,FALSE,FALSE,0)
  CI.method<-RGtk2::rGtkDataFrame(c("logit","normal"))
  comboboxCI.method <- RGtk2::gtkComboBox(CI.method)
  crtCI.method <- RGtk2::gtkCellRendererText()
  comboboxCI.method$packStart(crtCI.method)
  comboboxCI.method$addAttribute(crtCI.method, "text", 0)
  RGtk2::gtkComboBoxSetActive(comboboxCI.method,0)
  hbox$packStart(comboboxCI.method)
  
    #####################################################
  # Add separator
  hbox <- RGtk2::gtkHBoxNew(FALSE,8)
  vbox$packStart(hbox, FALSE, FALSE, 0)
  vbox$packStart(RGtk2::gtkHSeparatorNew(), FALSE, FALSE, 0)
  
  
 #------------------------------Contrast Statement-----------------------------#

  
  
   label <- RGtk2::gtkLabelNewWithMnemonic("Multiple Contrast Tests")
   hbox$packStart(label,FALSE,FALSE,0)
    hbox <- RGtk2::gtkHBoxNew(FALSE,8)
  vbox$packStart(hbox, FALSE, FALSE, 0)
  vbox$packStart(RGtk2::gtkHSeparatorNew(), FALSE, FALSE, 0)
  
  
    labelcontrast <- RGtk2::gtkLabelNewWithMnemonic("_Contrast Factor")
        hbox$packStart(labelcontrast,FALSE,FALSE,0)
        # Add entry in the second column; named "filename"
        filenamecontrast <- RGtk2::gtkEntryNew()
        filenamecontrast$setWidthChars(50)
        label$setMnemonicWidget(filenamecontrast)
        hbox$packStart(filenamecontrast,FALSE,FALSE,0)
        
 labelcontrastlist <- RGtk2::gtkLabelNewWithMnemonic("Contrast")
  hbox$packStart(labelcontrastlist,FALSE,FALSE,0)
  contrast.list<-RGtk2::rGtkDataFrame(c("","Dunnett","Tukey", "Sequen", "AVE", "Changepoint", 
                    "Williams", "Marcus", "McDermott", "UmbrellaWilliams", "GrandMean"))
  comboboxcontrast.list <- RGtk2::gtkComboBox(contrast.list)
  crtcontrast.list <- RGtk2::gtkCellRendererText()
  comboboxcontrast.list$packStart(crtcontrast.list)
  comboboxcontrast.list$addAttribute(crtcontrast.list, "text", 0)
  RGtk2::gtkComboBoxSetActive(comboboxcontrast.list,0)
  hbox$packStart(comboboxcontrast.list)
  
  labelscimethod<- RGtk2::gtkLabelNewWithMnemonic("SCI Method")
  hbox$packStart(labelscimethod,FALSE,FALSE,0)
  scimethod.list<-RGtk2::rGtkDataFrame(c("fisher","multi.t"))
  comboboxscimethod.list <- RGtk2::gtkComboBox(scimethod.list)
  crtscimethod.list <- RGtk2::gtkCellRendererText()
  comboboxscimethod.list$packStart(crtscimethod.list)
  comboboxscimethod.list$addAttribute(crtscimethod.list, "text", 0)
  RGtk2::gtkComboBoxSetActive(comboboxscimethod.list,0)
  hbox$packStart(comboboxscimethod.list)
  
  
  label3 <- RGtk2::gtkLabelNewWithMnemonic("_alpha")
  hbox$packStart(label3,FALSE,FALSE,0)
  # Add entry in the second column; named "filename3"
  filename3 <- RGtk2::gtkEntryNew()
  filename3$setWidthChars(10)
  filename3$setText(0.05)
  label3$setMnemonicWidget(filename3)
  hbox$packStart(filename3,FALSE,FALSE,0)
  
   label <- RGtk2::gtkLabelNewWithMnemonic("SCI Plot")
  hbox$packStart(label,FALSE,FALSE,0)
  toSCIPlot <- RGtk2::gtkCheckButton()
  hbox$packStart(toSCIPlot,FALSE,FALSE,0)

  
  #####################################################
  
  ############################################################
  
  
  
  # Add plot-option
  hbox <- RGtk2::gtkHBoxNew(FALSE,8)
  vbox$packStart(hbox, FALSE, FALSE, 0)
  label <- RGtk2::gtkLabelNewWithMnemonic("CI Plot")
  hbox$packStart(label,FALSE,FALSE,0)
  toPlot <- RGtk2::gtkCheckButton()
  hbox$packStart(toPlot,FALSE,FALSE,0)
  
  # Add Factor information
  #hbox <- RGtk2::gtkHBoxNew(FALSE,8)
  #vbox$packStart(hbox, FALSE, FALSE, 0)
  label <- RGtk2::gtkLabelNewWithMnemonic("Individual Factor Results")
  hbox$packStart(label,FALSE,FALSE,0)
  toFactors <- RGtk2::gtkCheckButton()
  hbox$packStart(toFactors,FALSE,FALSE,0)
  
   # Output Info
  label <- RGtk2::gtkLabelNewWithMnemonic("Output Information")
  hbox$packStart(label,FALSE,FALSE,0)
  toOutput <- RGtk2::gtkCheckButton()
   toOutput$active <- TRUE
  hbox$packStart( toOutput,FALSE,FALSE,0)
  label$setMnemonicWidget( toOutput)
  
  labelround <- RGtk2::gtkLabelNewWithMnemonic("_Output Decimals")
  hbox$packStart(labelround,FALSE,FALSE,0)
  # Add entry in the second column; named "filename3"
  filenameround <- RGtk2::gtkEntryNew()
  filenameround$setWidthChars(10)
  filenameround$setText(4)
  labelround$setMnemonicWidget(filenameround)
  hbox$packStart(filenameround,FALSE,FALSE,0)
  
  
   


 
   
  

   
    #####################################################
  # Add separator
  hbox <- RGtk2::gtkHBoxNew(FALSE,8)
  vbox$packStart(hbox, FALSE, FALSE, 0)
  vbox$packStart(RGtk2::gtkHSeparatorNew(), FALSE, FALSE, 0)
  
  # Add button
  the.buttons <- RGtk2::gtkHButtonBoxNew()
  the.buttons$setBorderWidth(5)
  vbox$add(the.buttons)
  the.buttons$setLayout("spread")
  the.buttons$setSpacing(40)
  buttonOK <- RGtk2::gtkButtonNewFromStock("gtk-ok")
  buttonLoad <- RGtk2::gtkButtonNewFromStock("Load Data")
  RGtk2::gSignalConnect(buttonOK, "clicked", performStatistics)
  RGtk2::gSignalConnect(buttonLoad, "clicked", getDirectory)
  the.buttons$packStart(buttonOK,fill=F)
  the.buttons$packStart(buttonLoad,fill=F)
  
  ##################################################
  
  buttonCancel <- RGtk2::gtkButtonNewFromStock("gtk-close")
  RGtk2::gSignalConnect(buttonCancel, "clicked", window$destroy)
  the.buttons$packStart(buttonCancel,fill=F)
}



rangtests <- function(formula, data,alpha=0.05, CI.method,
                      plot_CI=FALSE,effect,hypothesis,contrast = NULL, sci.method=c("fisher","multi.t"), 
                      Factor.Information=FALSE, info = TRUE, plot_SCI=FALSE,rounds){
  
  requireNamespace("RGtk2", quietly=TRUE)
  
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
 
  #n$pd <- c(H0pW$pd)
  #n$Var <- c(diag(H0pW$VBF))
  
  #CI <- Limits(c(H0pW$pd),H0pW$VBF,alpha,c(H0pW$N),CI.method)
  #n$Lower<- CI[,1]
  #n$Upper <- CI[,2]
  
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
    Descriptives <-round(data.frame(Rel.Effect=CI.Matrices[[i]]%*%H0pW$pd,
                              Std.Error= sqrt(c(diag(CI.Matrices[[i]]%*%H0pW$VBF%*%t(CI.Matrices[[i]])))/H0pW$N),
                              Lower=CILimits[,1],
                              Upper = CILimits[,2]),rounds)
                              
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
   
  
  if(plot_CI==TRUE){
    
    calculateGUIplot <- function() {
      plotting <- function(button, user.data) {
        
        ######################################################################
        # PLOT GUI for One way!!!
        if(nf ==1){
          Faktor = fac_names[1]}
        
        if(nf > 1){
          Faktor <- filename$getText()}
        
        Fak.split <- strsplit(Faktor,":")[[1]]
        l.Fak.split <- length(Fak.split)
        error <- NULL
        
        Title <- filename2$getText()
        lwd.gui <- as.numeric(filename3$getText())
        cex.gui <- as.numeric(filename4$getText())
        cex.lab.gui <- as.numeric(filename5$getText())
        cex.axis.gui <- as.numeric(filename6$getText())
        if(is.character(filename7$getText())){
        col.gui<- filename7$getText()}
        if(!is.character(filename7$getText())){
        col.gui <-as.numeric(filename7$getText())} 
  
        pch.gui <- as.numeric(filename8$getText())
       cex.ci.gui <- as.numeric(filename9$getText())
       xlab.gui <- filename10$getText()
       ylab.gui <- filename11$getText()
       
        if (!(Faktor %in% fac_names)) {
          error <- "Please enter a valid factor name"
        }
        
        if (!is.null(error)) {
          hbox <- RGtk2::gtkHBoxNew()
          vbox$packStart(hbox,FALSE,FALSE,0)
          label <- RGtk2::gtkLabel(error)
          hbox$packStart(label,FALSE,FALSE,0)
        }
        
        for(i in 1:n.hypotheses){
          
          if(names(Descriptive.Factors)[i] == Faktor){
            posP <- which(names(Descriptive.Factors)[i] == Faktor)
            DatenPlot <- data.frame(Descriptive.Factors[[i]])
            
            
           lower=DatenPlot$Lower
           upper=DatenPlot$Upper
           DatenPlot$pd<-DatenPlot$Rel.Effect
           if(nchar(xlab.gui)==0){
      text.X<-paste(names(DatenPlot[1]))}
       if(nchar(xlab.gui)>0){text.X<-xlab.gui}
      if(nchar(ylab.gui)==0){
        text.Ci <- paste((1 - alpha) * 100, "%", "Confidence Intervals for Relative Effects")}
        if(nchar(ylab.gui)>0){text.Ci<-ylab.gui}
            
            if (l.Fak.split==1){
               print(xyplot(pd ~ DatenPlot[,1], group=DatenPlot[,1],data = DatenPlot, 
                     type = 'p',  main=Title,
                     xlab=list(text.X,cex=cex.lab.gui),
                     ylab=list(text.Ci,cex=cex.lab.gui),
                      scales=list(x=list(cex=cex.axis.gui),y=list(cex=cex.axis.gui)),
                     col = col.gui, pch = pch.gui, cex=cex.gui, lwd=lwd.gui,  ylim = c(0, 1), 
                     upper = upper,
                     lower = lower,
                      panel = function(x, y, ...){
                     panel.superpose(x, y,panel.groups = function(x, y, upper, lower,subscripts, ..., font, fontface) {
                       upper <- upper[subscripts]
                       lower <- lower[subscripts]
                       panel.arrows(x, lower, x, upper,code=4,lwd=lwd.gui,col=col.gui)   
                       panel.points(x, lower,pch="_",cex=cex.ci.gui,lwd=lwd.gui,col=col.gui)
                       panel.points(x, upper,pch="_",cex=cex.ci.gui,lwd=lwd.gui,col=col.gui)
                     } , ...)
                       panel.xyplot(x, y, ...)}
                    ))
                             
                             }
            
            if (l.Fak.split==2){
               print(xyplot(pd ~ DatenPlot[,1]|DatenPlot[,2], 
                     group=DatenPlot[,1],data = DatenPlot, type = 'p', main=Title,
                     xlab=list(text.X,cex=cex.lab.gui),
                     ylab=list(text.Ci,cex=cex.lab.gui),
                     scales=list(x=list(cex=cex.axis.gui),y=list(cex=cex.axis.gui)),
                      col = col.gui, pch = pch.gui, cex=cex.gui, lwd=lwd.gui,  ylim = c(0, 1), 
                     upper = upper,
                     lower = lower,
                     par.strip.text =list(cex=cex.gui),
                     panel = function(x, y, ...){
                       panel.superpose(x, y,panel.groups = function(x, y, upper, lower,subscripts, ..., font, fontface) {
                         upper <- upper[subscripts]
                         lower <- lower[subscripts]
                         panel.arrows(x, lower, x, upper,code=4,lwd=lwd.gui,col=col.gui)   
                       panel.points(x, lower,pch="_",cex=cex.ci.gui,lwd=lwd.gui,col=col.gui)
                       panel.points(x, upper,pch="_",cex=cex.ci.gui,lwd=lwd.gui,col=col.gui)
                       }, ...)
                       panel.xyplot(x, y, ...)
                     }))
              
            }
            
            
            if (l.Fak.split==3){
              print(xyplot(pd ~ DatenPlot[,1]|DatenPlot[,2]*DatenPlot[,3], 
                     group=DatenPlot[,1],data = DatenPlot, type = 'p', main=Title,
                    xlab=list(text.X,cex=cex.lab.gui),
                     ylab=list(text.Ci,cex=cex.lab.gui),
                     scales=list(x=list(cex=cex.axis.gui),y=list(cex=cex.axis.gui)),
                      col = col.gui, pch = pch.gui, cex=cex.gui, lwd=lwd.gui,  ylim = c(0, 1), 
                     upper = upper,
                     lower = lower,
                     par.strip.text =list(cex=cex.gui),
                     panel = function(x, y, ...){
                       panel.superpose(x, y,panel.groups = function(x, y, upper, lower,subscripts, ..., font, fontface) {
                         upper <- upper[subscripts]
                         lower <- lower[subscripts]
                          panel.arrows(x, lower, x, upper,code=4,lwd=lwd.gui,col=col.gui)   
                          panel.points(x, lower,pch="_",cex=cex.ci.gui,lwd=lwd.gui,col=col.gui)
                          panel.points(x, upper,pch="_",cex=cex.ci.gui,lwd=lwd.gui,col=col.gui)
                       }, ...)
                       panel.xyplot(x, y, ...)
                     }))
              
            }
            
            if (l.Fak.split>=4){
              stop("4 and higher way interactions cannot be plotted!")
            }
            
          }
          
        }
        
      }
      
      # Create window
      window <- RGtk2::gtkWindow()
      # Add title
      window["title"] <- "Plot"
      
      # Add a frame
      frame <- RGtk2::gtkFrameNew("Please choose the factor you wish to plot (for interaction type something like group1:group2).")
      window$add(frame)
      
      # Create vertical container for file name entry
      vbox <- RGtk2::gtkVBoxNew(FALSE, 8)
      vbox$setBorderWidth(24)
      frame$add(vbox)
      # Add horizontal container for every widget line
      hbox <- RGtk2::gtkHBoxNew(FALSE, 8)
      vbox$packStart(hbox, FALSE, FALSE, 0)
      
      # Add label in first column
      if(nf >1){
        label <- RGtk2::gtkLabelNewWithMnemonic("_Factor")
        hbox$packStart(label,FALSE,FALSE,0)
        # Add entry in the second column; named "filename"
        filename <- RGtk2::gtkEntryNew()
        filename$setWidthChars(50)
        label$setMnemonicWidget(filename)
        hbox$packStart(filename,FALSE,FALSE,0)
      } 
      ############################################################
      
      # Add an horizontal container to specify parameters
      hbox <- RGtk2::gtkHBoxNew(FALSE,8)
      vbox$packStart(hbox, FALSE, FALSE, 0)
    
      
      label2 <- RGtk2::gtkLabelNewWithMnemonic("_Title")
      hbox$packStart(label2,FALSE,FALSE,0)
      filename2 <- RGtk2::gtkEntryNew()
      filename2$setWidthChars(10)
      label2$setMnemonicWidget(filename2)
      hbox$packStart(filename2,FALSE,FALSE,0)
      
          label10 <- RGtk2::gtkLabelNewWithMnemonic("_xlab")
       hbox$packStart(label10,FALSE,FALSE,0)
      filename10 <- RGtk2::gtkEntryNew()
      filename10$setWidthChars(10)
      label10$setMnemonicWidget(filename10)
      hbox$packStart(filename10,FALSE,FALSE,0)
      
      label11 <- RGtk2::gtkLabelNewWithMnemonic("_ylab")
       hbox$packStart(label11,FALSE,FALSE,0)
      filename11 <- RGtk2::gtkEntryNew()
      filename11$setWidthChars(10)
      label11$setMnemonicWidget(filename11)
      hbox$packStart(filename11,FALSE,FALSE,0)
      
        label5 <- RGtk2::gtkLabelNewWithMnemonic("_cex.lab")
      hbox$packStart(label5,FALSE,FALSE,0)
      filename5 <- RGtk2::gtkEntryNew()
      filename5$setWidthChars(10)
      filename5$setText(1.3)
      label5$setMnemonicWidget(filename5)
      hbox$packStart(filename5,FALSE,FALSE,0)
      
       label6 <- RGtk2::gtkLabelNewWithMnemonic("_cex.axis")
      hbox$packStart(label6,FALSE,FALSE,0)
      filename6 <- RGtk2::gtkEntryNew()
      filename6$setWidthChars(10)
      filename6$setText(1.3)
      label6$setMnemonicWidget(filename6)
      hbox$packStart(filename6,FALSE,FALSE,0)
      
          #####################################################
  # Add separator
  hbox <- RGtk2::gtkHBoxNew(FALSE,8)
  vbox$packStart(hbox, FALSE, FALSE, 0)
  vbox$packStart(RGtk2::gtkHSeparatorNew(), FALSE, FALSE, 0)
  
  
      
           label7 <- RGtk2::gtkLabelNewWithMnemonic("_color")
      hbox$packStart(label7,FALSE,FALSE,0)
      filename7 <- RGtk2::gtkEntryNew()
      filename7$setWidthChars(10)
      filename7$setText(1)
      label7$setMnemonicWidget(filename7)
      hbox$packStart(filename7,FALSE,FALSE,0)
      
      label3 <- RGtk2::gtkLabelNewWithMnemonic("_lwd")
      hbox$packStart(label3,FALSE,FALSE,0)
      filename3 <- RGtk2::gtkEntryNew()
      filename3$setWidthChars(10)
      filename3$setText(2)
      label3$setMnemonicWidget(filename3)
      hbox$packStart(filename3,FALSE,FALSE,0)
      
      label4 <- RGtk2::gtkLabelNewWithMnemonic("_cex")
      hbox$packStart(label4,FALSE,FALSE,0)
      filename4 <- RGtk2::gtkEntryNew()
      filename4$setWidthChars(10)
      filename4$setText(1.3)
      label4$setMnemonicWidget(filename4)
      hbox$packStart(filename4,FALSE,FALSE,0)
      
      label8 <- RGtk2::gtkLabelNewWithMnemonic("_pch")
      hbox$packStart(label8,FALSE,FALSE,0)
      filename8 <- RGtk2::gtkEntryNew()
      filename8$setWidthChars(10)
      filename8$setText(19)
      label8$setMnemonicWidget(filename8)
      hbox$packStart(filename8,FALSE,FALSE,0)
      
      label9 <- RGtk2::gtkLabelNewWithMnemonic("_cex.ci")
      hbox$packStart(label9,FALSE,FALSE,0)
      filename9 <- RGtk2::gtkEntryNew()
      filename9$setWidthChars(10)
      filename9$setText(4)
      label9$setMnemonicWidget(filename9)
      hbox$packStart(filename9,FALSE,FALSE,0)
      
  
      
      #labelcimethod <- RGtk2::gtkLabelNewWithMnemonic("CI Method")
      #hbox$packStart(labelcimethod,FALSE,FALSE,0)
      #cimethod<-RGtk2::rGtkDataFrame(c("logit","normal"))
      #comboboxcimethod <- RGtk2::gtkComboBox(cimethod)
      #combobox allowing to decide whether we want result as integer or double
      #crtcimethod <- RGtk2::gtkCellRendererText()
      #comboboxcimethod$packStart(crtcimethod)
      #comboboxcimethod$addAttribute(crtcimethod, "text", 0)
      #RGtk2::gtkComboBoxSetActive(comboboxcimethod,0)
      #hbox$packStart(comboboxcimethod)
      
      ############################################################
      
      # Add button
      the.buttons <- RGtk2::gtkHButtonBoxNew()
      the.buttons$setBorderWidth(5)
      vbox$add(the.buttons)
      the.buttons$setLayout("spread")
      the.buttons$setSpacing(40)
      buttonOK <- RGtk2::gtkButtonNewFromStock("gtk-ok")
      RGtk2::gSignalConnect(buttonOK, "clicked", plotting)
      the.buttons$packStart(buttonOK,fill=F)
      
      ##################################################
      
      buttonCancel <- RGtk2::gtkButtonNewFromStock("gtk-close")
      RGtk2::gSignalConnect(buttonCancel, "clicked", window$destroy)
      the.buttons$packStart(buttonCancel,fill=F)
    }
    
    
    #end of plot_ci==TRUE
    calculateGUIplot()
    
  }
  
  #-------------------------End of Plot GUI------------------------------------#
  
  #---------------------Compute MCTP for Main Effects--------------------------#
 if(is.null(contrast)){res.mctp <- NULL }
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
if(class(contrast[[2]])=="character"){
 CC <- contrMat(n=rep(10,nrow(CI.Matrices[[pos.contr]])),contrast[[2]]) 
 switch(hypothesis,H0p={
 res.mctp <- MCTP(H0pW$pd,CI.Matrices[[pos.contr]],H0pW$VBF,CC,n$Size,H0pW$placements,alpha,sci.method)},
 H0F={res.mctp <- MCTP.H0F(H0pW$pd,CI.Matrices[[pos.contr]],H0pW$VH0F,CC,n$Size,alpha,sci.method)})}
 if(class(contrast[[2]])!="character"){
 CC <- matrix(contrast[[2]],ncol=nrow(CI.Matrices[[pos.contr]]))
 switch(hypothesis,H0p={
 res.mctp <- MCTP(H0pW$pd,CI.Matrices[[pos.contr]],H0pW$VBF,CC,n$Size,H0pW$placements,alpha,sci.method) },
 H0F={res.mctp <- MCTP.H0F(H0pW$pd,CI.Matrices[[pos.contr]],H0pW$VH0F,CC,n$Size,alpha,sci.method) })}}
    colnames(res.mctp[[1]]) <- 1:ncol(res.mctp[[1]])
   rownames(res.mctp[[1]]) <- paste("C",1:nrow(CC),sep="")
   rownames(res.mctp[[2]]) <- paste("C",1:nrow(CC),sep="")
   res.mctp[[2]]<-round(res.mctp[[2]],rounds)
   res.mctp[[3]]<-round(res.mctp[[3]],rounds)
   res.mctp[[5]]<-res.mctp[[5]]}
  #-----------------------------Descriptive Output------------------------------#
  n$Rel.Effect <- c(H0pW$pd)
  n$Std.Error <- sqrt(c(diag(H0pW$VBF))/sum(n$Size))
  CI <- Limits(c(H0pW$pd),H0pW$VBF,alpha,c(H0pW$N),CI.method)
  n$Lower<- CI[,1]
  n$Upper <- CI[,2]
  rownames(n) <- 1:nrow(n)
  n[,(ncol(n)-3):ncol(n)]<- round(n[,(ncol(n)-3):ncol(n)],rounds)
  if (Factor.Information){res.factor.information <-Descriptive.Factors}
  if (!Factor.Information){res.factor.information <-NULL}
  
  if(hypothesis=="H0F"){
  if (nf ==1){
  colnames(KW) <- c("Statistic", "df", "p-Value")
  rownames(KW) <- Output.names
 result <- list(Call=formula,Descriptive=n, Wald.Type.Statistic = round(WTS,rounds), ANOVA.Type.Statistic=round(ATS,rounds), 
 Kruskal.Wallis.Test = round(KW,rounds),  MCTP=res.mctp, Factor.Information=res.factor.information)}
 if(nf!=1){result <- list(Call=formula,Descriptive=n, Wald.Type.Statistic =round(WTS,rounds), ANOVA.Type.Statistic=round(ATS,rounds), MCTP=res.mctp,Factor.Information=res.factor.information)}}
  if(hypothesis=="H0p"){
    result <- list(Call=formula,Descriptive=n, Wald.Type.Statistic=round(WTSp,rounds), ANOVA.Type.Statistic=round(ATSp,rounds), MCTP=res.mctp, Factor.Information=res.factor.information)} 
  
  
  #------------------------GUI PLOT SCI----------------------------------------#
  
  
  
  if(plot_SCI==TRUE){
    
    calculateGUISCIplot <- function() {
      SCIplotting <- function(button, user.data) {
        
                  ######################################################################
      
        error <- NULL
        
        Title <- filename2$getText()
        lwd.gui <- as.numeric(filename3$getText())
        cex.gui <- as.numeric(filename4$getText())
        cex.lab.gui <- as.numeric(filename5$getText())
        cex.axis.gui <- as.numeric(filename6$getText())
        if(is.character(filename7$getText())){
        col.gui<- filename7$getText()}
        if(!is.character(filename7$getText())){
        col.gui <-as.numeric(filename7$getText())} 
  
        pch.gui <- as.numeric(filename8$getText())
       cex.ci.gui <- as.numeric(filename9$getText())
       xlab.gui <- filename10$getText()
       ylab.gui <- filename11$getText()
       
        if (hypothesis=="H0F") {
          error <- "Simultaneous confidence intervals cannot be computed under H0F. Please choose hypothesis=H0p"
        }
        
        if (!is.null(error)) {
          hbox <- RGtk2::gtkHBoxNew()
          vbox$packStart(hbox,FALSE,FALSE,0)
          label <- RGtk2::gtkLabel(error)
          hbox$packStart(label,FALSE,FALSE,0)
        }
        
        if(hypothesis=="H0p"){
DatenPlot <- data.frame(res.mctp[[2]])
nc <- nrow(DatenPlot)
comp <- NULL
DatenPlot$comp <- 1:nc
lower=DatenPlot[,4]
upper=DatenPlot[,5]

if(nchar(xlab.gui)==0){
      text.X<-paste(names(DatenPlot[1]))}
       if(nchar(xlab.gui)!=0){text.X<-xlab.gui}
      if(nchar(ylab.gui)==0){
        text.Ci <- paste((1 - alpha) * 100, "%", "Simultaneous Confidence Intervals")}
        if(nchar(ylab.gui)!=0){text.Ci<-ylab.gui}
   print(xyplot(DatenPlot[,1] ~ comp, group=comp,data = DatenPlot, 
                     type = 'p',  main=Title,
                      xlab=list(text.X,cex=cex.lab.gui),
                     ylab=list(text.Ci,cex=cex.lab.gui),
                      scales=list(x=list(at=1:nc, labels=c(paste("C",1:nc,sep="")),cex=cex.axis.gui),y=list(cex=cex.axis.gui)),
                     col = col.gui, pch = pch.gui, cex=cex.gui, lwd=lwd.gui,  ylim = c(-1, 1), 
                     ,upper = upper,
                     lower = lower,
                      panel = function(x, y, ...){
                     panel.superpose(x, y,panel.groups = function(x, y, upper, lower,subscripts, ..., font, fontface) {
                       upper <- upper[subscripts]
                       lower <- lower[subscripts]
                       panel.arrows(x, lower, x, upper,code=4,lwd=lwd.gui,col=col.gui)   
                       panel.points(x, lower,pch="_",cex=cex.ci.gui,lwd=lwd.gui,col=col.gui)
                       panel.points(x, upper,pch="_",cex=cex.ci.gui,lwd=lwd.gui,col=col.gui)
                       panel.abline(h=0, col="red", lty=1, lwd=2)
                     } , ...)
                       panel.xyplot(x, y, ...)}
                    ))}
        
      }
      
      # Create window
      window <- RGtk2::gtkWindow()
      # Add title
      window["title"] <- "Plot of Simultaneous Confidence Intervals"
      
      # Add a frame
      frame <- RGtk2::gtkFrameNew("Plot Parameter Settings.")
      window$add(frame)
      
      # Create vertical container for file name entry
      vbox <- RGtk2::gtkVBoxNew(FALSE, 8)
      vbox$setBorderWidth(24)
      frame$add(vbox)
      # Add horizontal container for every widget line
      hbox <- RGtk2::gtkHBoxNew(FALSE, 8)
      vbox$packStart(hbox, FALSE, FALSE, 0)
      
  
      ############################################################
      
      # Add an horizontal container to specify parameters
      hbox <- RGtk2::gtkHBoxNew(FALSE,8)
      vbox$packStart(hbox, FALSE, FALSE, 0)
    
      
      label2 <- RGtk2::gtkLabelNewWithMnemonic("_Title")
      hbox$packStart(label2,FALSE,FALSE,0)
      filename2 <- RGtk2::gtkEntryNew()
      filename2$setWidthChars(10)
      label2$setMnemonicWidget(filename2)
      hbox$packStart(filename2,FALSE,FALSE,0)
      
          label10 <- RGtk2::gtkLabelNewWithMnemonic("_xlab")
       hbox$packStart(label10,FALSE,FALSE,0)
      filename10 <- RGtk2::gtkEntryNew()
      filename10$setWidthChars(10)
      label10$setMnemonicWidget(filename10)
      hbox$packStart(filename10,FALSE,FALSE,0)
      
      label11 <- RGtk2::gtkLabelNewWithMnemonic("_ylab")
       hbox$packStart(label11,FALSE,FALSE,0)
      filename11 <- RGtk2::gtkEntryNew()
      filename11$setWidthChars(10)
      label11$setMnemonicWidget(filename11)
      hbox$packStart(filename11,FALSE,FALSE,0)
      
        label5 <- RGtk2::gtkLabelNewWithMnemonic("_cex.lab")
      hbox$packStart(label5,FALSE,FALSE,0)
      filename5 <- RGtk2::gtkEntryNew()
      filename5$setWidthChars(10)
      filename5$setText(1.3)
      label5$setMnemonicWidget(filename5)
      hbox$packStart(filename5,FALSE,FALSE,0)
      
       label6 <- RGtk2::gtkLabelNewWithMnemonic("_cex.axis")
      hbox$packStart(label6,FALSE,FALSE,0)
      filename6 <- RGtk2::gtkEntryNew()
      filename6$setWidthChars(10)
      filename6$setText(1.3)
      label6$setMnemonicWidget(filename6)
      hbox$packStart(filename6,FALSE,FALSE,0)
      
          #####################################################
  # Add separator
    hbox <- RGtk2::gtkHBoxNew(FALSE,8)
    vbox$packStart(hbox, FALSE, FALSE, 0)
    vbox$packStart(RGtk2::gtkHSeparatorNew(), FALSE, FALSE, 0)
  
  
      
           label7 <- RGtk2::gtkLabelNewWithMnemonic("_color")
      hbox$packStart(label7,FALSE,FALSE,0)
      filename7 <- RGtk2::gtkEntryNew()
      filename7$setWidthChars(10)
      filename7$setText(1)
      label7$setMnemonicWidget(filename7)
      hbox$packStart(filename7,FALSE,FALSE,0)
      
      label3 <- RGtk2::gtkLabelNewWithMnemonic("_lwd")
      hbox$packStart(label3,FALSE,FALSE,0)
      filename3 <- RGtk2::gtkEntryNew()
      filename3$setWidthChars(10)
      filename3$setText(2)
      label3$setMnemonicWidget(filename3)
      hbox$packStart(filename3,FALSE,FALSE,0)
      
      label4 <- RGtk2::gtkLabelNewWithMnemonic("_cex")
      hbox$packStart(label4,FALSE,FALSE,0)
      filename4 <- RGtk2::gtkEntryNew()
      filename4$setWidthChars(10)
      filename4$setText(1.3)
      label4$setMnemonicWidget(filename4)
      hbox$packStart(filename4,FALSE,FALSE,0)
      
      label8 <- RGtk2::gtkLabelNewWithMnemonic("_pch")
      hbox$packStart(label8,FALSE,FALSE,0)
      filename8 <- RGtk2::gtkEntryNew()
      filename8$setWidthChars(10)
      filename8$setText(19)
      label8$setMnemonicWidget(filename8)
      hbox$packStart(filename8,FALSE,FALSE,0)
      
      label9 <- RGtk2::gtkLabelNewWithMnemonic("_cex.ci")
      hbox$packStart(label9,FALSE,FALSE,0)
      filename9 <- RGtk2::gtkEntryNew()
      filename9$setWidthChars(10)
      filename9$setText(4)
      label9$setMnemonicWidget(filename9)
      hbox$packStart(filename9,FALSE,FALSE,0)
      
      
      ############################################################
      
      # Add button
      the.buttons <- RGtk2::gtkHButtonBoxNew()
      the.buttons$setBorderWidth(5)
      vbox$add(the.buttons)
      the.buttons$setLayout("spread")
      the.buttons$setSpacing(40)
      buttonOK <- RGtk2::gtkButtonNewFromStock("gtk-ok")
      RGtk2::gSignalConnect(buttonOK, "clicked", SCIplotting)
      the.buttons$packStart(buttonOK,fill=F)
      
      ##################################################
      
      buttonCancel <- RGtk2::gtkButtonNewFromStock("gtk-close")
      RGtk2::gSignalConnect(buttonCancel, "clicked", window$destroy)
      the.buttons$packStart(buttonCancel,fill=F)
    }
    
    
    #end of plot_SCI==TRUE
    calculateGUISCIplot()
    
  }
  
  #-------------------------End of Plot GUI------------------------------------#
  
  #---------------------------------RESULT and Output--------------------------#
 output.hypothesis=switch(hypothesis, H0F={"Distribution Functions"}, H0p={"Relative Effects"})
 output.ranks=switch(effect, weighted={"Global Ranks"}, unweighted={"Pseudo-Ranks"})
 output.confidence=switch(CI.method, logit={"Logit-Transformation"}, normal={"Normal Approximation"})

 
 result$output <- list(output.hypothesis=output.hypothesis,output.ranks=output.ranks,output.confidence=output.confidence,
 output.info=info,output.contrast=contrast,output.sci.method=sci.method,output.alpha=alpha)
  
     class(result) <- "rankFD"

  return(result)
}

