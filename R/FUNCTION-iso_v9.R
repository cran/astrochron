### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### iso: isolate a portion of the data set for analysis 
###       - (SRM: January 30, 2012; March 29, 2012; Nov. 14, 2012; 
###               May 20, 2013; August 13, 2013; April 3, 2014; April 9-10, 2014
###               May 15, 2014; Nov. 19, 2014; February 4, 2015; May 11, 2016)
###########################################################################

iso <- function (dat,xmin=NULL,xmax=NULL,col=2,logx=F,logy=F,genplot=T,verbose=T)
{
  if(verbose) cat("\n----- ISOLATE STRATIGRAPHIC DATA BY LOCATION -----\n")
  
  dat <- data.frame(dat)
  ncols <- length(dat)
  ipts <- length(dat[,1]) 
  if(verbose) cat(" * Number of data points=", ipts,"\n")
  if(verbose) cat(" * Number of columns=", ncols,"\n")
  mi <- min(dat[1])
  ma <- max(dat[1])
  if(verbose) cat(" * Minimum=",mi,", Maximum=",ma,"\n")


  if(is.null(xmin) && is.null(xmax)) 
      {
        if(verbose) cat("\n * Graphically select portion of data series to isolate\n")
        if(verbose) cat("   Stop by pressing ESC-key (Mac) or STOP button (Windows)\n")
        a=idPts(cbind(dat[,1],dat[,col]),output=1,logx=logx,logy=logy,verbose=F)
        xmin=min(a[1])
        xmax=max(a[1])
      } 
     if(is.null(xmin) && !is.null(xmax)) 
      {
        if(verbose) cat("\n * xmax set to:",xmax,"\n")
        if(verbose) cat(" * Graphically select lowest datum (xmin) to isolate\n")
        if(verbose) cat("   Stop by pressing ESC-key (Mac) or STOP button (Windows)\n")
        a=idPts(cbind(dat[,1],dat[,col]),xmax=xmax,output=1,logx=logx,logy=logy,verbose=F) 
        xmin=min(a[1])
      }  
     if(is.null(xmax) && !is.null(xmin))  
      {
        if(verbose) cat("\n * xmin set to:",xmin,"\n")
        if(verbose) cat(" * Graphically select highest datum (xmax) to isolate\n")
        if(verbose) cat("   Stop by pressing ESC-key (Mac) or STOP button (Windows)\n")
        a=idPts(cbind(dat[,1],dat[,col]),xmin=xmin,output=1,logx=logx,logy=logy,verbose=F) 
        xmax=max(a[1])
      }  
  
  if(verbose) cat(" * Isolating data between",xmin,"and",xmax,"\n")  
  dat <- subset(dat, (dat[1] >= xmin) & (dat[1] <= xmax) ) 

  newpts=length(dat[,1])
  if(verbose) cat(" * Number of data points following culling=",newpts,"\n")

# determine if any missing entries (NA) exist
     numNA=sum(is.na(dat))
     if(numNA > 0) 
       {
         cat("\n  WARNING:", numNA,"empty entries are present in your data frame.\n")
       } 
  
if(length(dat)==2 && genplot==T)
 {

  plotdat = subset(dat, !(dat[2] == "NA"))
### plots
  par(mfrow=c(2,2))
  if(logx && logy) setlog="xy"
  if(logx && !logy) setlog="x"
  if(!logx && logy) setlog="y"
  if(!logx && !logy) setlog=""
  plot(plotdat,cex=0.5,xlab="Location",ylab="Value",main="Stratigraphic Series",log=setlog); lines(plotdat)
### plot the denisty and the histogram together
  hist(plotdat[,2],freq=F,xlab="Value",main="Distribution of Isolated Values"); lines(density(plotdat[,2], bw="nrd"),col="red"); grid()
### boxplot
  boxplot(plotdat[,2],ylab="Value",main="Boxplot for Isolated Values")
### Normal probabilty plot (Normal Q-Q Plot)
  qqnorm(plotdat[,2]); qqline(plotdat[,2], col="red");grid()

  } 
  
  
if(length(dat)>2 && genplot==T)
  {
    ncols1 = length(dat)-1
    nrows = ceiling(sqrt(ncols1-1))
    par(mfrow = c(nrows, ceiling(ncols1/nrows)))

#    for (i in 2:ncols1) 
    for (i in 2:(ncols1+1)) 
      {
        plotdat = subset(dat[,i], !(dat[,i] == "NA"))
        xlab = names(dat)[i]
        hist(plotdat,freq=F,xlab=xlab,main=""); lines(density(plotdat, bw="nrd"),col="red"); grid()
      }
  }  
  
  return(dat)

### END function iso
}
