### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2021 Stephen R. Meyers
###
###########################################################################
### function periodogram - (SRM: January 26, 2012; April 28, 2012; 
###                         October 10, 2012, November 20, 2012; May 20-24, 2013; 
###                         June 5, 2013; June 13, 2013; July 30-31, 2013;
###                         August 3, 2013; August 7-10, 2013; Nov. 26, 2013;
###                         October 18, 2014; October 22, 2014; January 21, 2015;
###                         September 10, 2015; November 20-29, 2017; January 14, 2021)
###
### simple unwindowed periodogram
###########################################################################

periodogram <- function (dat,padfac=2,demean=T,detrend=F,nrm=1,background=0,output=0,f0=F,fNyq=T,xmin=0,xmax=Nyq,pl=1,genplot=T,verbose=T)
{

if(verbose) cat("\n----- CALCULATING PERIODOGRAM FOR STRATIGRAPHIC SERIES -----\n")

   d <- data.frame(dat)
   npts <- length(d[,1]) 
   dt = d[2,1]-d[1,1]

# error checking 
   if(dt<0)
     { 
       if (verbose) cat("\n * Sorting data into increasing height/depth/time, removing empty entries\n")
       d <- d[order(d[,1], na.last = NA, decreasing = F), ]
       dt <- d[2,1]-d[1,1]
       npts <- length(d[,1])
     }
   dtest <- d[2:npts,1]-d[1:(npts-1),1] 
   epsm=1e-9
   if( (max(dtest)-min(dtest)) > epsm ) 
     {
       cat("\n**** ERROR: sampling interval is not uniform.\n")
       stop("**** TERMINATING NOW!")
     }

  if (verbose) 
    {
      cat(" * Number of data points in stratigraphic series:",npts,"\n")
      cat(" * Stratigraphic series length (space or time):",(npts-1)*dt,"\n")
      cat(" * Sampling interval (space or time):",dt,"\n")
    }

###########################################################################
### remove mean and linear trend
###########################################################################
   if (demean) 
     { 
       dave <- colMeans(d[2])
       d[2] <- d[2] - dave
       if(verbose) cat(" * Mean value removed=",dave,"\n")
     }
     
   if (!demean && verbose) cat(" * Mean value NOT subtracted\n")

### use least-squares fit to remove linear trend
    if (detrend) 
      {
      lm.0 <- lm(d[,2] ~ d[,1])
      d[2] <- d[2] - (lm.0$coeff[2]*d[1] + lm.0$coeff[1])
      if(verbose) cat(" * Linear trend removed. m=",lm.0$coeff[2],"b=",lm.0$coeff[1],"\n")
      }

   if (!detrend && verbose) cat(" * Linear trend NOT subtracted\n") 
 
### Calculate Nyquist
   Nyq <- 1/(2*dt)
#### Calculate Rayleigh Frequency
   Ray <- 1/(npts*dt)

### pad with zeros if like (power of 2 not required!)
### also convert from data frame to numeric
  pad <- as.numeric(d[,2])
  if(padfac>1) pad <- append( pad, rep(0,(npts*padfac-npts)) )
  
### add another zero if we don't have an even number of data points, so Nyquist exists.   
   if((npts*padfac)%%2 != 0) pad <- append(pad,0)

if(verbose)
 {
   cat(" * Nyquist frequency:",Nyq,"\n")
   cat(" * Rayleigh frequency:",Ray,"\n")
   cat(" * Padded to",length(pad),"points\n")
 }

### new resulting frequency grid   
   nf = length(pad)
   df <- 1/(nf*dt)
   freq <- double(nf)
 
### set frequency index vector 
    i <- seq(1,nf,by=1)
     
### take fft
   ft <- fft(pad)
### apply normalization   
   if(nrm == 1) ft=ft/npts
### caculate power 
   pwr <- Mod(ft)^2
### caculate amplitude
   amp <- sqrt(pwr)
### calculate phase
   phase <- atan2(Im(ft),Re(ft))
### now make frequency vector (negative frequencies are not assigned correctly here; we will discard)
   freq <- df*(i-1)
### put all results into a common data frame
   fft.out <- data.frame(cbind(freq,amp,pwr,phase))
### extract results from 0 to positive nyquist (recall that the negative frequencies were not 
###   assigned correctly; they are listed in 'freq' as > Nyquist)
   if(fNyq) fft.out <- subset(fft.out,(fft.out[,1] <= Nyq))
   if(!fNyq) fft.out <- subset(fft.out,(fft.out[,1] < Nyq))
   if(!f0) fft.out <- subset(fft.out,(fft.out[,1] > 0))

   colnames(fft.out) <- c('Frequency','Amplitude','Power','Phase')

# AR(1) noise model spectrum
   if(background==1)
    {
### what is the estimated AR1 coefficient?
      lag0 <- d[1:(npts-1),2]
      lag1 <- d[2:npts,2]
      rho <- cor(lag0,lag1)
### Calculate Raw red noise spectrum
### "So" is the average power. This can be determined from the white noise variance
###  as So = var/(1-rho^2), where rho is the lag-1 coeff
###  We will determine average power directly from measured spectrum
###  Note that we are potentially omitting f0 and/or fNyq, depending upon user selection
      So = mean(fft.out[,3])
      AR = So * (1-(rho^2)) / (  1 - (2*rho*cos(pi*fft.out[,1]/Nyq)) + (rho^2) )
      dofAR = 2
      chiAR <-  (fft.out[,3]/AR) * dofAR
      chiCLAR <- pchisq(chiAR, df=dofAR)
### 90, 95 and 99% CL for power plot
      AR1_90 <- AR*qchisq(0.9, df=dofAR)/dofAR
      AR1_95 <- AR*qchisq(0.95, df=dofAR)/dofAR
      AR1_99 <- AR*qchisq(0.99, df=dofAR)/dofAR
      fft.out <- data.frame(cbind(fft.out[,1],fft.out[,2],fft.out[,3],fft.out[,4],chiCLAR*100,AR,AR1_90,AR1_95,AR1_99))
      colnames(fft.out) <- c('Frequency','Amplitude','Power','Phase','AR1_CL','AR1_Fit','AR1_90_power','AR1_95_power','AR1_99_power')
    }

# power law noise model spectrum
   if(background==2 && f0 && verbose) cat(" * WARNING: Cannot conduct power law fit, because f(0) is included.\n")

   if(background==2 && !f0)
    { 
      specIn=cbind(fft.out[1],fft.out[,3])
      resFit=pwrLawFit(specIn,dof=2,output=1,genplot=F,verbose=F)
      pwrLaw <- resFit$PowerLaw_fit
      pwrLaw_90 <- resFit$CL_90
      pwrLaw_95 <- resFit$CL_95
      pwrLaw_99 <- resFit$CL_99 
      pwrLaw_CL <- resFit$PowerLaw_CL
      fft.out <- data.frame(cbind(fft.out[,1],fft.out[,2],fft.out[,3],fft.out[,4],pwrLaw_CL,pwrLaw,pwrLaw_90,pwrLaw_95,pwrLaw_99))
      colnames(fft.out) <- c('Frequency','Amplitude','Power','Phase','PwrLaw_CL','PwrLaw_Fit','PwrLaw_90_power','PwrLaw_95_power','PwrLaw_99_power')
    }


### output the Fourier coefficients if desired
   if (output==2)
     {
      fc.out <- data.frame(cbind(freq,Re(ft),Im(ft)))
      if(fNyq) fc.out <- subset(fc.out,(fc.out[,1] <= Nyq))
      if(!fNyq) fc.out <- subset(fc.out,(fc.out[,1] < Nyq))
      if(!f0) fc.out <- subset(fc.out,(fc.out[,1] > 0))
      colnames(fc.out)[1] <- 'Frequency'
      colnames(fc.out)[2] <- 'Real Coeff.'
      colnames(fc.out)[3] <- 'Imag. Coeff.'
     }

   if(genplot)
    {
      dev.new(title = "Periodogram results", height = 7, width = 9)
      par(mfrow=c(2,2))
### plot the results
      plot(d,type="l", col="blue",ylab="Value", xlab="Location", main="Stratigraphic Series",bty="n")
      if (pl == 2) 
        {
          plot(fft.out[,1],fft.out[,3], type="l",col="red", ylab="Power", xlim=c(xmin,xmax),xlab="Frequency", main="Periodogram Power",bty="n")
          if(background > 0)
           { 
             lines(fft.out[,1],fft.out[,7],lwd=1,lty=3)
             lines(fft.out[,1],fft.out[,8],lwd=1,lty=3)
             lines(fft.out[,1],fft.out[,9],lwd=1,lty=3)
           }
        }
# do not plot f(0) if using log spectrum
      if (pl == 1)  
        { 
         ii = which(fft.out[,1] > 0)
         plot(fft.out[ii,1],log(fft.out[ii,3]), type="l",col="red", ylab="Log Power", xlim=c(xmin,xmax), xlab="Frequency", main="Log Periodogram Power",bty="n")
         if(background > 0)
           { 
             lines(fft.out[ii,1],log(fft.out[ii,7]),lwd=1,lty=3)
             lines(fft.out[ii,1],log(fft.out[ii,8]),lwd=1,lty=3)
             lines(fft.out[ii,1],log(fft.out[ii,9]),lwd=1,lty=3)
           }
        }
     plot(fft.out[,1],fft.out[,2], type="l",col="red", ylab="Amplitude", xlim=c(xmin,xmax), xlab="Frequency", main="Periodogram Amplitude",bty="n")
     plot(fft.out[,1],fft.out[,4], type="l",col="red", ylab="Phase", xlim=c(xmin,xmax), xlab="Frequency", main="Periodogram Phase",bty="n")
   }
   
if (output==1)  return(fft.out)
if (output==2)  return(fc.out)
   
#### END function periodogram
}