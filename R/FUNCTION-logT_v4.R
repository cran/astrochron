### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### function logT - (SRM: January 30, 2012; December 12, 2012; April 24, 2013
###                       May 20, 2013; June 16, 2015)
###
### apply log transformation
###########################################################################

logT <- function (dat,c=0,opt=1,genplot=T,verbose=T) 
{

   if(verbose) cat("\n----- PERFORMING log TRANSFORM OF STRATIGRAPHIC SERIES-----\n")
   dat <- data.frame(dat)
# a is a constant to add before log transformation
   npts <- length(dat[,1]) 
   if(verbose) cat(" * Number of data points=", npts,"\n")

   if(opt==1) trans <- log(dat[,2] + c )
   if(opt==2) trans <- log10(dat[,2] + c )
 
   out <- data.frame(cbind(dat[,1],trans))

   if(genplot)
     {
### plot data series. Note, cex is the factor by which to increase or decrease default symbol size
      par(mfrow=c(2,2))
      plot(out, cex=.5)
      lines(out)
### plot the denisty and the histogram together
      hist(trans,freq=F) 
      lines(density(trans, bw="nrd"),col="red"); grid()
### boxplot
      boxplot(trans)
### Normal probabilty plot (Normal Q-Q Plot)
      qqnorm(trans); qqline(trans, col="red");grid()
     }
   
   return(out)

### END function logT
}
