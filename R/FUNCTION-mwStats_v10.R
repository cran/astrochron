### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2021 Stephen R. Meyers
###
###########################################################################
### mwStats : moving window average, median, and standard deviation, 
###               allowing for dynamic adustment of window so that it has a 
###               constant duration in time or space (SRM: April 28, 2014; 
###               April 29, 2014; June 10, 2016; May 24-25, 2017; 
###               May 31, 2017; June 13-14, 2017; August 23, 2020; 
###               Sept. 27-29, 2020; Nov. 26, 2020; Jan. 14, 2021)
###
###########################################################################

mwStats <- function (dat,cols=NULL,win=NULL,conv=1,ends=F,CI=0,output=T,genplot=T,verbose=T)

{

  if(verbose) cat("\n----- CALCULATING MOVING AVERAGE, MEDIAN AND VARIANCE USING DYNAMIC WINDOW-----\n")

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
  res=double(npts)
  res2=double(npts)
  res3=double(npts)
  if(CI>0) 
    {  
      res4=double(npts)
      alpha=(1 - (0.01*CI) )/2
      conf=1-alpha
    } 
   n=double(npts)
   locCenter=loc$center
  
# now loop over all windows
  for (i in 1:npts)
   {
# shuffle in data
     if(loc$n1[i]==loc$n2[i]) cat("**** WARNING: Only one data point in window",i,":",loc$center[i],"\n")
      
     y <- dat[loc$n1[i]:loc$n2[i],2]
# number of data points in window
     n[i]=loc$n2[i]-loc$n1[i]+1
     res[i] = mean(y)
     res2[i] = median(y)
     res3[i] = var(y)
     if(CI>0 && n[i]>1) res4[i] = qt(p=conf,df=n[i]-1) * sqrt(res3[i])/sqrt(n[i])      
     if(CI>0 && n[i]==1) res4[i] = NA
    }

# if conv=1 and ends=T calculate estimates for endpoints <0.5*win that were not evaluated
#  these will be assigned as a constant value
    if(conv==1 && ends)
      {
        cat("**** WARNING: The endpoints use a smaller window; points=", 1:(loc$n1[1]-1),(loc$n2[npts]+1):ipts,"\n")
# note: no warning is presently given if there is only one data point in the window
        nB=max(which(dat[,1]<loc$center[1])) 
        nE=ipts-min(which(dat[,1]>loc$center[npts]))+1       
        n=c(rep(nB,nB),n,rep(nE,nE))
        y <- dat[1:nB,2]
        resB = mean(y)
        resB2 = median(y)
        resB3 = var(y)
        if(CI>0) 
          {
            if(n[1]==1) resB4 = NA
            if(n[1]>1) resB4 = qt(p=conf,df=n[1]-1) * sqrt(resB3)/sqrt(n[1])
          }  
        y <- dat[(ipts-nE+1):ipts,2]
        resE = mean(y)
        resE2 = median(y)
        resE3 = var(y)
        if(CI>0) 
          {
            if(n[ipts]==1) resE4 = NA
            if(n[ipts]>1) resE4 = qt(p=conf,df=n[ipts]-1) * sqrt(resE3)/sqrt(n[ipts])
          }  
         
         res=c(rep(resB,n[1]),res,rep(resE,n[ipts]))
         res2=c(rep(resB2,n[1]),res2,rep(resE2,n[ipts]))
         res3=c(rep(resB3,n[1]),res3,rep(resE3,n[ipts])) 
         if(CI>0) res4=c(rep(resB4,n[1]),res4,rep(resE4,n[ipts]))
         locCenter=c(dat[1:n[1],1],loc$center,dat[(ipts-n[ipts]+1):ipts,1])
      }


 if(genplot)
   {
     dev.new(height = 7.7, width = 5.7, title = "mwStats Results")
     par(mfrow=c(4,1),mar=c(3.1, 4.1, 4.1, 2.1))
     
     xmin=min(dat[,1])
     xmax=max(dat[,1])
     plot(dat[,1],dat[,2],type="l",xlab="",ylab="",main="Stratigraphic Series",bty="n",xlim=c(xmin,xmax))
     mtext("Location",side=1,line=2,cex=0.7)
     mtext("Value",side=2,line=2,cex=0.7)

     par(mar=c(4.1, 4.1, 3.1, 2.1))
     if(CI==0) 
       {
         ymin=min(res,res2)
         ymax=max(res,res2)
         plot(locCenter,res,ylim=c(ymin,ymax),xlim=c(xmin,xmax),type="l",lwd=2,xlab="",ylab="",main="Moving Window Average (black) and Median (red)",bty="n")
         lines(locCenter,res2,col="red",lwd=2)
       }  
     if(CI>0)
       {
         ymin=min(res-res4,na.rm=TRUE)
         ymax=max(res+res4,na.rm=TRUE)         
         plot(locCenter,res,ylim=c(ymin,ymax),xlim=c(xmin,xmax),type="l",lwd=2,xlab="",ylab="",main="Moving Window Average and Confidence Interval",bty="n")
         lines(locCenter,res+res4,lty=2)
         lines(locCenter,res-res4,lty=2)
       }
       
     mtext("Center of window",side=1,line=2,cex=0.7)
     mtext("Value",side=2,line=2,cex=0.7)

     plot(locCenter,res3,xlim=c(xmin,xmax),type="l",xlab="",ylab="",main="Variance",lwd=2,bty="n")
     mtext("Center of window",side=1,line=2,cex=0.7)
     mtext("Variance",side=2,line=2,cex=0.7)
     
     plot(locCenter,n,xlim=c(xmin,xmax),type="l",xlab="",ylab="",main="Number of Data Points in Window",bty="n")
     mtext("Center of window",side=1,line=2,cex=0.7)
     mtext("# Points",side=2,line=2,cex=0.7)
   }

    if(output)
     {
       if(CI==0) out = data.frame (cbind (locCenter,res,res2,res3,n) ) 
       if(CI>0) out = data.frame (cbind (locCenter,res,res2,res3,n,res4) ) 
       colnames(out)[1] <- 'Center_win'
       colnames(out)[2] <- 'Average'
       colnames(out)[3] <- 'Median'
       colnames(out)[4] <- 'Variance'
       colnames(out)[5] <- 'Points'
       if(CI>0) colnames(out)[6] <- 'Conf_interval'
       return(out)
     }  

### END function mwStats
}