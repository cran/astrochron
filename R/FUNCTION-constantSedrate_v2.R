### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2014 Stephen R. Meyers
###
###########################################################################
### constantSedrate function - (SRM: November 24, 2012; June 18, 2013; Sept. 25, 2013)
###
### Apply a constant sedimentation rate to transform a depth series to elapsed time
###########################################################################


constantSedrate <- function (dat,sedrate,begin=0,timeDir=1,genplot=T,verbose=T)
{

   dat = data.frame(dat)
   npts <- length(dat[,1]) 
   if(verbose)   { cat("\n----- APPLY CONSTANT SEDIMENTATION RATE MODEL-----\n") }   
   if(verbose)   { cat(" * Number of data points=", npts,"\n") }

### sort to ensure increasing depth/height/time
   if(verbose) cat("\n * Sorting data, removing empty entries\n")
   dat <- dat[order(dat[1],na.last=NA,decreasing=F),]

# vector of sampling intervals
   dx <- dat[2:npts,1]-dat[1:(npts-1),1]
   mindx=min(dx)
   maxdx=max(dx)
   if(verbose) cat("\n * Minimum sampling interval (spatial)=", mindx,"\n")
   if(verbose) cat(" * Maximum sampling interval (spatial)=", maxdx,"\n")
   mindt=mindx/sedrate
   maxdt=maxdx/sedrate
   if(verbose) cat("\n * Minimum sampling interval (temporal)=", mindt,"\n")
   if(verbose) cat(" * Maximum sampling interval (temporal)=", maxdt,"\n")   
    
   time <- double(npts)      
   time[1] = begin
   time[2:npts] <- cumsum(dx/sedrate) + begin

### assign to data frame
    out <- as.data.frame(cbind(time,dat[,2]))

# if you want time to decrease with depth/height
    if(timeDir== -1) out = flip(out,begin=begin,genplot=F,verbose=F)
    
    if(genplot)
     {
      par(mfrow=c(2,1))
      plot(dat, cex=.5,xlab="Location",ylab="Value",main="Spatial Series",bty="n")
      lines(dat)
      plot(out, cex=.5,xlab="Time",ylab="Value",main="Tuned Series",bty="n")
      lines(out)
     }
     
     return(data.frame(out))

### END function constantSedrate
}
