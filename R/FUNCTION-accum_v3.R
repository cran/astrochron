### This function is a component of astro: An R Package for Astrochronology
### Copyright (C) 2024 Stephen R. Meyers
###
###########################################################################
### function accum - (SRM: January 15, 2014; April 6, 2024; June 10, 2024; Sept. 4 2024)
###
### calculate accumulation rates
###########################################################################


accum <- function (dat,sedrate=NULL,density=NULL,genplot=T,verbose=T)
{

if(verbose) cat("\n----- CALCULATING ACCUMULATION RATES -----\n")

# dat = spatial series (2 columns): depth or time, concentration in wt.% (0-100)
# sedrate = sedimentation rates (2 columns): depth or time, sedrate (cm/ka)
# density = dry bulk density (2 columns): depth or time, density (g/cm3)

# make sure the input series are data frames
   dat=data.frame(dat)
   out=dat
   
   if(!is.null(sedrate))sedrate=data.frame(sedrate)
   if(is.null(sedrate))
    {
      cat("\n**** ERROR: sedimentation rate not provided\n")
      stop("**** TERMINATING NOW!")
    }
    
    if(!is.null(density))density=data.frame(density)
    if(is.null(density))
    {
      cat("\n**** WARNING: dry bulk density not provided. Will use default of 2.65 g/cm3.\n")
      density=data.frame(c(2.65))
    }
    
   npts <- length(dat[,1]) 
   if(verbose) cat(" * Number of geochem data points=", npts,"\n")

   numSed <- length(sedrate[,1])
   if(verbose) cat(" * Number of sedimentation rates=", numSed,"\n")
   
   numDen <- length(density[,1])
   if(verbose) cat(" * Number of density measurements=", numDen,"\n")
   
### sort to ensure increasing depth/height/time
   if(verbose) cat(" * Sorting datasets into ensure increasing order, removing empty entries\n")
   dat <- dat[order(dat[,1],na.last=NA,decreasing=F),]
   if(numSed >1) sedrate <- sedrate[order(sedrate[,1],na.last=NA,decreasing=F),]
   if(numDen >1) density <- density[order(density[,1],na.last=NA,decreasing=F),]
   
### check for duplicate depths/heights in dat
   dx1=dat[2:npts,1]-dat[1:(npts-1),1]
   if(min(dx1) == 0)
     {
       cat("\n**** WARNING: duplicate depth/height datum found in dat\n")
     }  

### check for duplicate depths/heights in sedrate
   if(numSed >1) 
     {
       dx2=sedrate[2:numSed,1]-sedrate[1:(numSed-1),1]
       if(min(dx2) == 0)
        {
          cat("\n**** ERROR: duplicate depth/height datum found in sedimentation rates\n")
          stop("**** TERMINATING NOW!")
        }
      }    

### check for duplicate depths/heights in density
   if (numDen > 1) 
    {
      dx3=density[2:numDen,1]-density[1:(numDen-1),1]
      if(min(dx3) == 0)
       {
         cat("\n**** ERROR: duplicate depth/height datum found in density measurements\n")
         stop("**** TERMINATING NOW!")
       } 
     }   

    if(numSed == 1) cat("\n**** WARNING: one sedimentation rate applied to entire record\n")
    if(numDen == 1) cat("\n**** WARNING: one bulk density applied to entire record\n")

### if only one sedimentation rate and bulk density are supplied
    if(numSed == 1 && numDen == 1)
     {
       bulk=sedrate[1,1]*density[1,1]
       out[2] <- (out[2]/100)*bulk   
# now redefine bulk for plotting
       b1=c(min(out[,1]),max(out[,1]))
       b2=c(bulk,bulk)
       bulk <- data.frame(cbind(b1,b2))
     }     

### use piecewise linear interpolation to estimate sedimentation rate and bulk density
# resample density on dat[,1]
      if(numDen > 1) density<-resample(density,dat[,1],genplot=F,verbose=F)
# resample sedrates on dat[,1]      
      if(numSed > 1) sedrate<-resample(sedrate,dat[,1],genplot=F,verbose=F)

# find min and max depth that are present in all data sets
      minDat=min(dat[,1])
      minSed=min(sedrate[,1])
      minDen=min(density[,1])
      if(numSed>1 && numDen >1) isoMin=max(minDat,minSed,minDen)
      if(numSed>1 && numDen == 1) isoMin=max(minDat,minSed)
      if(numSed==1 && numDen >1) isoMin=max(minDat,minDen)
      maxDat=max(dat[,1])
      maxSed=max(sedrate[,1])
      maxDen=max(density[,1])
      if(numSed>1 && numDen >1) isoMax=min(maxDat,maxSed,maxDen)
      if(numSed>1 && numDen ==1) isoMax=min(maxDat,maxSed)
      if(numSed==1 && numDen >1) isoMax=min(maxDat,maxDen)

# isolate common portion from all records
      if(numSed>1 || numDen >1) out=iso(dat,isoMin,isoMax,genplot=F,verbose=F)
      if(numDen>1) density2=iso(density,isoMin,isoMax,genplot=F,verbose=F)
      if(numSed>1) sedrate2=iso(sedrate,isoMin,isoMax,genplot=F,verbose=F)
      if(numDen ==1) density2=density
      if(numSed ==1) sedrate2=sedrate
        
# calculate accumulation rates
      if(numSed>1 || numDen >1) 
       {
         if(numSed>1 && numDen >1) bulk= sedrate2[2]*density2[2] 
         if(numSed==1 && numDen >1) bulk= sedrate2[1,1]*density2[2] 
         if(numSed>1 && numDen ==1) bulk= sedrate2[2]*density2[1,1]
         out[2] <- (out[2]/100)*bulk[1]      
# now define bulk for plotting
         bulk <- data.frame(cbind(out[,1],bulk[,1]))
       }     
            
     if(genplot)
      {
### plot data series. Note, cex is the factor by which to increase or decrease default symbol size
       par(mfrow=c(3,1))
       plot(dat, xlab="Location",ylab="Concentration (wt.%)",main="Concentration Series",bty="n",lwd=2,cex.axis=1.3,cex.lab=1.3,cex.main=1.4,cex=0.5)
       lines(dat)
       plot(bulk, xlab="Location",ylab="Accumulation Rate (g/cm2/ka)",main="Bulk Accumulation Series",type="l",bty="n",lwd=2,cex.axis=1.3,cex.lab=1.3,cex.main=1.4,cex=0.5)
       plot(out, xlab="Location",ylab="Accumulation Rate (g/cm2/ka)",main="Component Accumulation Series",bty="n",lwd=2,cex.axis=1.3,cex.lab=1.3,cex.main=1.4,cex=0.5)
       lines(out)
      }
     
     return(data.frame(out))

### END function accum
}
