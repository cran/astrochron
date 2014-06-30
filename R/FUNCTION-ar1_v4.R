### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2014 Stephen R. Meyers
###
###########################################################################
### AR1 : make AR1 noise - (SRM: January 24, 2012; April 29, 2012; 
###                         April 27-29, 2013; May 20-25, 2013; July 31, 2013)
###########################################################################

ar1 <- function (npts=1024, dt=1, mean=0, sdev=1, rho=0.9, genplot=T,verbose=T)
{

   if(verbose) cat("\n ----- GENERATING AR1 MODEL -----\n")
### Generate normal deviates, mean = 0
### this is "white noise"
   white <- rnorm(npts,sd=sdev)
### generate AR(1) red noise
   ar1 <- 1:(length(white))
### initialize ar1[1]
   ar1[1] <- white[1]
   for (i in 2:(length(white)))  
### multiply previous value by coeff, then add innovation
      { 
        ar1[i] <- rho*ar1[i-1]+white[i]
      }
### generate time axis
    ta <- 1:(length(white))
### change time axis from unit spacing of 1 to desired value
    ta <- (ta*dt) - dt

### assign to data frame
    noise <- as.data.frame(cbind(ta,ar1))
### plot noise model
    if(genplot)
      {
        par(mfrow=c(2,2))
        plot(noise, cex=.5, xlab="Location", ylab="Noise Value", main="AR1 Noise Series"); lines(noise)
### plot the denisty and the histogram together
        hist(ar1,freq=F,xlab="Noise Value",main="Histogram of Noise Values"); lines(density(ar1, bw="nrd"),col="red"); grid()
### boxplot
        boxplot(ar1,ylab="Noise Value",main="Boxplot of Noise Values")
### Normal probabilty plot (Normal Q-Q Plot)
        qqnorm(ar1,xlab="Noise Value"); qqline(ar1, col="red"); grid()
       }
       
### what is the estimated AR1 coefficient?
    lag0 <- ar1[1:(npts-1)]
    lag1 <- ar1[2:npts]

    rho = cor(lag0,lag1)
    if(verbose) cat(" * Estimated AR1 coefficient =",rho,"\n")
    
# derived from EQ 2.45 of Mulelsee book, page 57
    rho_unbias= (rho*(npts-1) + 1 ) / (npts - 4)
    if(verbose) cat(" * Unbiased AR1 =",rho_unbias,"\n") 
        
    return(noise)
    
### END function ar1
}

