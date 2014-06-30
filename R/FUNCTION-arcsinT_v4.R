### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2014 Stephen R. Meyers
###
###########################################################################
### function arcsinT - (SRM: January 25, 2012; April 30, 2012; May 20, 2013;
###                          June 5, 2013)
###
### apply arcsin transformation
###########################################################################

arcsinT <- function (dat,genplot=T,verbose=T) 
{
   if(verbose) cat("\n----- PERFORMING ARCSIN TRANSFORM OF STRATIGRAPHIC SERIES-----\n")
   dat <- data.frame(dat)
   npts <- length(dat[,1]) 
   if(verbose) cat(" * Number of data points=", npts,"\n")

# rescale y from 0-1
# find max and min y
      ymax = max(dat[,2]) 
      ymin = min(dat[,2])
      slope= 1/(ymax-ymin)
      b= 1 - (slope * ymax)
      if(verbose)
        { 
          cat(" * minimum =", ymin,"\n")
          cat(" * maximum =", ymax,"\n")
          cat(" * slope =", slope,"\n")
          cat(" * y-intercept =", b,"\n")
        }
        
      y1 <- double(npts)
      for ( i in 1:npts )
        {
          y1[i] = (slope*dat[i,2] + b)
# correct for possible very small negative numbers due to round off errors
          if (y1[i] < 0) { y1[i] = 0 }
        }

#     arcsin transform
      trans <- double(npts)
      for (i in 1:npts)
       {
        trans[i] <- asin( sqrt(y1[i]) )
       }
 
      out <- cbind(dat[,1],trans)

  if(genplot)
   {
### plot data series. Note, cex is the factor by which to increase or decrease default symbol size
    par(mfrow=c(2,2))
    plot(out, cex=.5,xlab="Location",ylab="Transformed Value", main="New Stratigraphic Series")
    lines(out)
### plot the denisty and the histogram together
    hist(trans,freq=F,main="Histogram of Transformed Data",xlab="Transformed Value") 
    lines(density(trans, bw="nrd"),col="red"); grid()
### boxplot
    boxplot(trans,ylab="Transformed Value",main="Boxplot of Transformed Data")
### Normal probabilty plot (Normal Q-Q Plot)
    qqnorm(trans); qqline(trans, col="red"); grid()
   }
   
   return(data.frame(out))

### END function arcsinT
}
