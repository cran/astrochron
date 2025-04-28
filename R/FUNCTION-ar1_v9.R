### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2025 Stephen R. Meyers
###
###########################################################################
### AR1 : make AR1 noise - (SRM: January 24, 2012; April 29, 2012; 
###                         April 27-29, 2013; May 20-25, 2013; July 31, 2013;
###                         November 7, 2014; January 20, 2015; April 9, 2015;
###                         July 22,2016; January 12, 2025)
###########################################################################

ar1 <- function (npts=1024, dt=1, mean=0, sdev=1, rho=0.9, shuffle=F, nsim=1, genplot=T, verbose=T)
{

   if(verbose) cat("\n ----- GENERATING AR1 SURROGATES -----\n")
 
# initialize matrix for simulation results
   red <- double(npts*nsim)
   dim(red) <- c(npts,nsim)
   
# start simulation loop
   for (i in 1:nsim)
     { 
### Generate normal deviates, mean = 0, std dev=1
      if(!shuffle) { white <- rnorm(npts,mean=0,sd=1) }
      if(shuffle)
        {
# generate normal deviates
         r1=rnorm(npts,mean=0,sd=1)
# sort and output index 
         index=sort.int(r1, method=c("shell"),index.return=T)$ix
# generate another set of normal deviates
         r2=rnorm(npts,mean=0,sd=1)
# now reorder r2 using index
         white=r2[index]
        }
### generate AR(1) red noise
      red[1,i] <- white[1]
### multiply previous value by coeff, then add innovation
      for (ii in 2:npts) { red[ii,i] <- rho*red[ii-1,i]+white[ii] }
### standardize to specified mean and variance
###  note, first mean is function call, second instance is desired mean value
       red[, i] = red[, i] - mean(red[, i])
       red[, i] = red[, i] * sdev/sd(red[, i])
       red[, i] = red[, i] + mean
     }
    
   if(nsim==1)
     {     
       if(verbose)
        {
### what is the estimated AR1 coefficient?
# npts is the length of vector 'red', ensure mean value is zero
         d0=red-mean(red)
         rho=sum(d0[1:(npts-1)] * d0[2:npts]) / sum(d0^2)
         cat("\n * Estimated AR1 coefficient=",rho,"\n")
# derived from EQ 2.45 of Mulelsee book, page 57
         rho_unbias= (rho*(npts-1) + 1 ) / (npts - 4)
         if(verbose) cat(" * Unbiased AR1 =",rho_unbias,"\n") 
        }
### generate time axis
       ta <- 1:npts
### change time axis from unit spacing of 1 to desired value
       ta <- (ta*dt) - dt
### assign to data frame
       noise <- as.data.frame(cbind(ta,red))
### plot noise model
      if(genplot)
       {
        par(mfrow=c(2,2))
        plot(noise, cex=.5, xlab="Location", ylab="Noise Value", main="AR1 Noise Series"); lines(noise)
### plot the denisty and the histogram together
        hist(red,freq=F,xlab="Noise Value",main="Histogram of Noise Values"); lines(density(red, bw="nrd0"),col="red"); grid()
### boxplot
        boxplot(red,ylab="Noise Value",main="Boxplot of Noise Values")
### Normal probabilty plot (Normal Q-Q Plot)
        qqnorm(red,xlab="Noise Value"); qqline(red, col="red"); grid()
        }
      return(noise)
     }
     
    if(nsim>1) return(red)
      
    
### END function ar1
}

