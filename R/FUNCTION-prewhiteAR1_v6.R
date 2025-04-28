### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2025 Stephen R. Meyers
###
###########################################################################
### function prewhiteAR1 - (SRM: January 28, 2012; March 1, 2012; Feb. 27, 2013;
###                          June 27, 2013; Sept. 25, 2013; Sept. 10, 2015; 
###                          July 22, 2016; January 4, 2025)
###
### prewhiten with AR1 filter
###########################################################################

prewhiteAR1 <- function (dat,setrho=NULL,bias=F,genplot=T,verbose=T) 
{

   if(verbose) cat("\n----- PREWHITENING STRATIGRAPHIC SERIES WITH AR1 FILTER -----\n")

   dat <- data.frame(dat)
   npts <- length(dat[,1])
   dt <- dat[2,1]-dat[1,1]

# error checking 
   dtest <- dat[2:npts,1]-dat[1:(npts-1),1] 
   epsm=1e-9
   if( (max(dtest)-min(dtest)) > epsm ) 
     {
       cat("\n**** ERROR: sampling interval is not uniform.\n")
       stop("**** TERMINATING NOW!")
     }

   if(verbose) 
     {
       cat(" * Number of data points=", npts,"\n")
       cat(" * Sampling interval (space or time):",dt,"\n")
     }  

### what is the estimated AR1 coefficient?
   if(is.null(setrho))
    {
# npts is the length of data vector 'dat', ensure mean value is zero
      d0=dat[,2]-mean(dat[,2])
      rho=sum(d0[1:(npts-1)] * d0[2:npts]) / sum(d0^2)
      if(verbose) cat("\n * Estimated AR1 coefficient=",rho,"\n")

      if(bias) 
       { 
# derived from EQ 2.45 of Mulelsee book, page 57
         rho= (rho*(npts-1) + 1 ) / (npts - 4)
         if(verbose) cat(" * Unbiased AR1 =",rho,"\n")
        } 
      setrho=rho
    }
   
### now prewhiten
    prewhite=dat[2:npts,2] - setrho * dat[1:(npts-1),2]
    ipts=npts-1
### what is the new raw estimate for the AR(1) coefficient?
    d0=prewhite-mean(prewhite)
    coeff_est_prewhite=sum(d0[1:(ipts-1)] * d0[2:ipts]) / sum(d0^2)
    if(verbose) cat(" * Prewhitened AR1 coefficient=",coeff_est_prewhite,"\n")
         
    out <- data.frame(cbind(dat[1:ipts,1],prewhite))

### plot data
  if(genplot)
   {    
### plot data series. Note, cex is the factor by which to increase or decrease default symbol size
     par(mfrow=c(2,2))
     plot(dat,cex=.5,main="Original data series",xlab="Position",ylab="Value")
     lines(dat)
     plot(out,cex=.5,main="Whitened data series",xlab="Position",ylab="Whitened value")
     lines(out)
### plot the denisty and the histogram together
     hist(prewhite,freq=F,main="Distribution of whitened data",xlab="Whitened value"); lines(density(prewhite, bw="nrd0"),col="red")
### boxplot
#     boxplot(prewhite,main="Boxplot of whitened data")
### Normal probabilty plot (Normal Q-Q Plot)
     qqnorm(prewhite); qqline(prewhite, col="red")
   }
   
   return(out)

### END function prewhiteAR1
}
