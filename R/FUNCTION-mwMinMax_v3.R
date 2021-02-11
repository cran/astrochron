### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2021 Stephen R. Meyers
###
###########################################################################
### mwMinMax : Identify minimum and maximum values in dynamic moving window 
###                    (SRM: May 27, 2020; Sept. 28, 2020; Jan. 14, 2021)
###
###########################################################################

mwMinMax <- function (dat,cols=NULL,win=NULL,conv=1,output=T,genplot=T,verbose=T)

{
  if(verbose) cat("\n----- PERFORMING Min/Max EVALUATION USING DYNAMIC WINDOW-----\n")

# ensure we have a data frame
  dat=data.frame(dat)
  
  if(is.null(cols) && length(dat) != 2)
    {
      if(verbose) cat("\n**** ERROR: number of variables > 1.  Please specify which variable to analyze.\n")  
      stop("**** TERMINATING NOW")
    }
# no need to do anything if is.null(cols) && length(dat) == 2
  if(!is.null(cols) && length(dat) == 2) 
    {
     if(verbose) cat("\n**** WARNING: You only have two columns, will ignore cols.\n") 
    } 
  if(!is.null(cols) && length(cols) != 1 && length(dat) > 2) 
    {
      if(verbose) cat("\n**** WARNING: cols not specified correctly, will use data from column 2.\n")  
# reassign dat
      dat <- data.frame( cbind(dat[,1],dat[,2]) )
    }  
  if(!is.null(cols) && length(cols) == 1 && length(dat) > 2) 
    {
# reassign dat
      dat <- data.frame( cbind(dat[,1],dat[,cols[1]]) )
    }  

### sort to ensure increasing depth/height/time
   if(verbose) cat(" * Sorting into increasing order, removing empty entries\n")
   dat <- dat[order(dat[,1],na.last=NA,decreasing=F),]

   ipts <- length(dat[,1]) 
   if(verbose) cat(" * Number of data points=", ipts,"\n")  

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
   
# determine location of each window
  loc <- mwin(dat,win,conv=conv,verbose=F)
  
  npts=length(loc$n1)
  resMax=double(npts)
  resMin=double(npts)
  n=double(npts)
  
# now loop over all windows
  for (i in 1:npts)
   {
# shuffle in data
     if(loc$n1[i]==loc$n2[i]) 
       {
         cat("**** WARNING: Only one data point in window",i,":",loc$center[i],"\n")
         y <- dat[loc$n1[i],2]
# number of data points in window
         n[i]=1
         resMax[i] = y
         resMin[i] = y
       }
     if(loc$n1[i]!=loc$n2[i])
       {
         y <- dat[loc$n1[i]:loc$n2[i],2]
# number of data points in window
         n[i]=loc$n2[i]-loc$n1[i]+1
         resMax[i] = max(y)
         resMin[i] = min(y)        
       }
    }

 if(genplot)
   {
     dev.new(height = 7.7, width = 5.7, title = "mwMinMax Results")
     par(mfrow=c(4,1),mar=c(3.1, 4.1, 4.1, 2.1))
     
     ymin=min(resMin)
     ymax=max(resMax)
     xmin=min(dat[,1])
     xmax=max(dat[,1])
     plot(dat[,1],dat[,2],type="l",xlab="",ylab="",main="Stratigraphic Series",bty="n",xlim=c(xmin,xmax))
     mtext("Location",side=1,line=2,cex=0.7)
     mtext("Value",side=2,line=2,cex=0.7)

     par(mar=c(4.1, 4.1, 3.1, 2.1))
     plot(loc$center,resMin,ylim=c(ymin,ymax),xlim=c(xmin,xmax),type="l",lwd=2,xlab="",ylab="",main="Maximum and Minimum Values",bty="n")
     lines(loc$center,resMax)
     mtext("Center of window",side=1,line=2,cex=0.7)
     mtext("Value",side=2,line=2,cex=0.7)
 
     plot(loc$center,resMax-resMin,xlim=c(xmin,xmax),type="l",lwd=2,xlab="",ylab="",main="Max-Min Values",bty="n")
     
     plot(loc$center,n,xlim=c(xmin,xmax),type="l",xlab="",ylab="",main="Number of Data Points in Window",bty="n")
     mtext("Center of window",side=1,line=2,cex=0.7)
     mtext("# Points",side=2,line=2,cex=0.7)
   }

    
    if(output)
     {
       out = data.frame (cbind (loc$center,resMin,resMax,resMax-resMin,n) ) 
       colnames(out)[1] <- 'Center_win'
       colnames(out)[2] <- 'Min'
       colnames(out)[3] <- 'Max'
       colnames(out)[4] <- 'Max-Min'
       colnames(out)[5] <- 'Points'
       return(out)
     }  

### END function mwMinMax
}
