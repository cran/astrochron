### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2017 Stephen R. Meyers
###
###########################################################################
### mwCor: moving window correlation, allowing for dynamic
###       adustment of window so it has a constant duration
###       in time or space (SRM: March 16, 2013; May 20, 2013; 
###       June 5-6, 2013; February 13, 2014; May 25, 2017; May 20, 2017
###       May 31, 2017; June 13-14, 2017; July 18, 2017; July 29, 2017)
###########################################################################

mwCor <- function (dat,cols=NULL,win=NULL,conv=1,cormethod=1,output=T,pl=1,genplot=T,verbose=T)

{

  if(verbose) cat("\n----- CALCULATING MOVING WINDOW CORRELATION WITH DYNAMIC WINDOW-----\n")

# ensure we have a data frame
  dat=data.frame(dat)
  ipts <- length(dat[,1]) 
### sorting will be done in mwin
   if(verbose) cat(" * Number of data points=", ipts,"\n")  
  
  if(is.null(cols) && length(dat) > 3)
    {
      if(verbose) cat("\n**** ERROR: number of variables > 2.  Please specify which variables to analyze.\n")  
      stop("**** TERMINATING NOW")
    }
 
  if(!is.null(cols) && length(cols) == 2) 
    {
# reassign dat
      dat <- data.frame( cbind(dat[,1],dat[,cols[1]],dat[,cols[2]]) )
    }
    
  if(!is.null(cols) && length(cols) != 2) 
    {
      if(verbose) cat("\n**** ERROR: you must specify two variables to examine.\n")  
      stop("**** TERMINATING NOW")
    }  

   if(is.null(win)) 
    {
      win=(dat[ipts,1]-dat[1,1])/10
      if(verbose) cat(" * Setting window size to default value of", win,"\n")  
    }  
  
   if(win >= (dat[ipts,1]-dat[1,1])) 
    {
      if(verbose) cat("\n**** ERROR: window length is >= data series duration.\n")  
      stop("**** TERMINATING NOW")    
    } 
  
   if(win < 0) 
    {
      if(verbose) cat("\n**** ERROR: window length is a negative number, must be positive!\n")  
      stop("**** TERMINATING NOW")    
    } 
  
# determine location of each window
  loc <- mwin(dat,win,conv=conv,verbose=F)
  
  npts=length(loc$n1)
  r=double(npts)
  n=double(npts)
  
# now loop over all windows
  for (i in 1:npts)
   {
# shuffle in data
     if(loc$n1[i]==loc$n2[i]) stop("**** ERROR: Only one data point in window! Terminating!")
     x <- dat[loc$n1[i]:loc$n2[i],1]
     y <- dat[loc$n1[i]:loc$n2[i],2]
     z <- dat[loc$n1[i]:loc$n2[i],3]
# number of data points in window     
     n[i]=loc$n2[i]-loc$n1[i]+1
     if(cormethod == 1) r[i] = cor(y,z,method=c("pearson"))
     if(cormethod == 2) r[i] = cor(y,z,method=c("spearman"))
     if(cormethod == 3) r[i] = cor(y,z,method=c("kendall"))
    }

 if(genplot)
   {
     dev.new(height = 7.7, width = 5.7, title = "mwCor Results")
     par(mfrow=c(4,1),mar=c(3.1, 4.1, 4.1, 2.1))
     
     xmin=min(dat[,1])
     xmax=max(dat[,1])
     plot(dat[,1],dat[,2],type="l",xlim=c(xmin,xmax),xlab="",ylab="",main="Variable A",bty="n")
     mtext("Location",side=1,line=2,cex=0.7)
     mtext("Value",side=2,line=2,cex=0.7)

     par(mar=c(4.1, 4.1, 3.1, 2.1))
     plot(dat[,1],dat[,3],xlim=c(xmin,xmax),type="l",xlab="",ylab="",main="Variable B",bty="n")
     mtext("Location",side=1,line=2,cex=0.7)
     mtext("Value",side=2,line=2,cex=0.7)

    if(pl==1)
     {
      plot(loc$center,r,xlim=c(xmin,xmax),type="l",xlab="",ylab="",main="Moving Window Correlation Coefficient",col="red",lwd=2,bty="n")
      mtext("Center of window",side=1,line=2,cex=0.7)
      mtext("Correlation",side=2,line=3,cex=0.7)
      mtext("Coefficient",side=2,line=2,cex=0.7)
     }
     
     if(pl==2)
      {
       plot(0,0,xlim=c(xmin,xmax),ylim=c(min(r),max(r)),type="l",xlab="",ylab="",main="Moving Window Correlation Coefficient",bty="n")
       for(i in 1:npts)
        {
          lines(c(dat[loc$n1[i],1],dat[loc$n2[i],1]),c(r[i],r[i]),col="red")
          points(c(dat[loc$n1[i]:loc$n2[i],1]),c(rep(r[i],loc$n2[i]-loc$n1[i]+1)))
        }  
          
       mtext("Location",side=1,line=2,cex=0.7)
       mtext("Correlation",side=2,line=3,cex=0.7)
       mtext("Coefficient",side=2,line=2,cex=0.7)
      }
     
     plot(loc$center,n,xlim=c(xmin,xmax),type="l",xlab="",ylab="",main="Number of Data Points in Window",lwd=2,bty="n")
     mtext("Center of window",side=1,line=2,cex=0.7)
     mtext("# Points",side=2,line=2,cex=0.7)
   }

    if(output)
     {
       out = data.frame (cbind (loc$center,r,n) ) 
       colnames(out)[1] <- 'Location'
       colnames(out)[2] <- 'Correlation'
       colnames(out)[3] <- 'Points'
       return(out)
     }  

### END function mwCor
}