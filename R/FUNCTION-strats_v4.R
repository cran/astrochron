### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2014 Stephen R. Meyers
###
###########################################################################
### strats - (SRM: January 25, 2012; March 9, 2012; April 24, 2013; 
###                May 20, 2013; August 6, 2013; August 9, 2013)
### 
### summary statistics for stratigraphic data series
###########################################################################

strats <- function (dat)
{

   cat("\n----- DETERMINING SUMMARY STATISTICS FOR STRATIGRAPHIC SERIES -----\n")
   dat <- data.frame(dat)
   npts <- length(dat[,1]) 
   cat(" * Number of data points=",npts,"\n")

### now evaluate sampling statistics
   t1<-dat[1:(npts-1),1]
   t2<-dat[2:(npts),1]
   dt=t2-t1
   dtMin=min(dt) 
   dtMax=max(dt)
   dtMean=mean(dt)     
   dtMedian=median(dt)

   cat("\n * Mean sampling interval=", dtMean,"\n")
   cat(" * Median sampling interval=",dtMedian,"\n")
   cat(" * Maximum sampling interval=",dtMax,"\n")
   cat(" * Minimum sampling interval=", dtMin,"\n")
   cat(" * Variance of sampling interval=",var(dt),"\n")

   epsm=1e-10
   if( dtMax-dtMin > epsm ) 
    {
### now output summary plots of dt
      par(mfrow=c(2,3))
      plot(t1,dt, cex=.5, xlab="location",ylab="dt",main="dt by location")
      lines(t1,dt)
### plot the denisty and the histogram together
       hist(dt,freq=F, xlab="dt",main="histogram of dt"); lines(density(dt, bw="nrd"),col="red"); grid() 
    }  
    
### do not output summary plots of dt    
   if( dtMax-dtMin < epsm ) par(mfrow=c(2,2))
    
### statistics on data values
    cat("\n * Mean data value=",mean(dat[,2]),"\n")
    cat(" * Median data value=",median(dat[,2]),"\n")
    cat(" * Maximum data value=",max(dat[,2]),"\n") 
    cat(" * Minimum data value=",min(dat[,2]),"\n") 
    cat(" * Variance of data values=",var(dat[,2]),"\n")
    
### now plot y variable of data series. Note, cex is the factor by which to increase or decrease default symbol size
    plot(dat, cex=.5, xlab="location",ylab="data value",main="data")
    lines(dat)
### plot the denisty and the histogram together
    hist(dat[,2],freq=F, ylab="data value",main="distribution of values"); lines(density(dat[,2], bw="nrd"),col="red"); grid()
### boxplot
    boxplot(dat[,2], ylab="data value",main="boxplot of values")
### Normal probabilty plot (Normal Q-Q Plot)
    qqnorm(dat[,2]); qqline(dat[,2], col="red"); grid()
   
### END function strats
}
