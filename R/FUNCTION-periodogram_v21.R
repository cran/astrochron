### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2025 Stephen R. Meyers
###
###########################################################################
### function periodogram - (SRM: January 26, 2012; April 28, 2012; 
###                         October 10, 2012, November 20, 2012; May 20-24, 2013; 
###                         June 5, 2013; June 13, 2013; July 30-31, 2013;
###                         August 3, 2013; August 7-10, 2013; Nov. 26, 2013;
###                         October 18, 2014; October 22, 2014; January 21, 2015;
###                         September 10, 2015; November 20-29, 2017; 
###                         January 14, 2021; October 4, 2021; May 26, 2022; 
###                         June 18-19, 2022; July 28-31, 2022; August 18, 2022;
###                         November 6, 2024; January 4, 2025)
###
### simple unwindowed periodogram
###########################################################################

periodogram <- function (dat,padfac=2,demean=T,detrend=F,nrm=1,background=0,medsmooth=0.2,bc=F,output=0,f0=F,fNyq=T,xmin=0,xmax=Nyq,pl=1,genplot=T,check=T,verbose=T)
{

if(verbose) cat("\n----- CALCULATING PERIODOGRAM FOR STRATIGRAPHIC SERIES -----\n")

   d <- data.frame(dat)
   npts <- length(d[,1]) 
   dt = d[2,1]-d[1,1]

   if(check)
    {
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

# padfac must be >= 1
   if(padfac <= 0) padfac=1  
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

# Locate real components for zero, nyquist, first neg freq., (zero-df)
   izero = 1
   nyqfreq = 0.5*nf + 1
   negfreq = 0.5*nf + 2
   minusdf = nf
   
### set frequency index vector 
   i <- seq(1,nf,by=1)
### assign positive frequencies out to Nyquist
   freq <- df*(i[1:nyqfreq]-1)
### assign negative frequencies
   f2 <- ( (nf/2) - (1:(minusdf-negfreq+1) ) ) * df * -1
   freq <- append(freq,f2)

### put all results into a common data frame
   fft.out <- data.frame(cbind(freq,amp,pwr,phase))
### extract results from 0 to positive nyquist
   if(fNyq) fft.out <- fft.out[1:nyqfreq,]
   if(!fNyq) fft.out <- fft.out[1:(nyqfreq-1),]
   if(!f0) fft.out <- fft.out[-1,]

   colnames(fft.out) <- c('Frequency','Amplitude','Power','Phase')


### nominal number of independent frequencies, ignoring f(0) and f(Nyq)
   nf1=Nyq/Ray 
   if(!bc) probs=c(0.9,0.95,0.99)
   if(bc) probs=c(1.0-((1-0.9)/nf1),1.0-((1-0.95)/nf1),1.0-((1-0.99)/nf1))
   dof=2

# AR(1) noise model spectrum
   if(background==1)
    {
### what is the estimated conventional AR1 coefficient?
# npts is the length of data vector 'd', which may or may not have been
# detrended/demeaned, so ensure mean value is zero
      d0=d[,2]-mean(d[,2])
      rho=sum(d0[1:(npts-1)] * d0[2:npts]) / sum(d0^2)
      if(verbose) cat(" * Estimated AR1 coefficient =",rho,"\n")
### Calculate Raw red noise spectrum
### "So" is the average power. This can be determined from the white noise variance
###  as So = var/(1-rho^2), where rho is the lag-1 coeff
###  We will determine average power directly from measured spectrum
###  Note that we are potentially omitting f0 and/or fNyq, depending upon user selection
      So = mean(fft.out[,3])
      AR = So * (1-(rho^2)) / (  1 - (2*rho*cos(pi*fft.out[,1]/Nyq)) + (rho^2) )
      chiAR <-  (fft.out[,3]/AR) * dof
      chiCLAR <- pchisq(chiAR, df=dof)
      if(bc) chiCLAR <- 1-pmin(1,(1-chiCLAR)*nf1)
### 90, 95 and 99% CL for power plot, 2 degrees of freedom
      AR1_90 <- AR*qchisq(probs[1], df=dof)/dof
      AR1_95 <- AR*qchisq(probs[2], df=dof)/dof
      AR1_99 <- AR*qchisq(probs[3], df=dof)/dof
      fft.out <- data.frame(cbind(fft.out[,1],fft.out[,2],fft.out[,3],fft.out[,4],chiCLAR*100,AR,AR1_90,AR1_95,AR1_99))
      colnames(fft.out) <- c('Frequency','Amplitude','Power','Phase','AR1_CL','AR1_Fit','AR1_90_power','AR1_95_power','AR1_99_power')
    }


# power law noise model spectrum
   if(background==2 && f0 && verbose) cat(" * WARNING: Cannot conduct power law fit, because f(0) is included.\n")

   if(background==2 && !f0)
    { 
     specfit=cbind(fft.out[,1],fft.out[,3])
#     if(fNyq) specfit=cbind(fft.out[-length(fft.out[,1]),1],fft.out[-length(fft.out[,1]),3])
      
# fit line to log10(power) and log10(freq), which is the equivalent
#  of an exponential continuum model
#   predict power given frequency
     lm.1 <- lm(log10(specfit[,2]) ~ log10(specfit[,1]))

     if (verbose) 
      {
       cat(" * Slope from log(power) vs. log(frequency) fit (m):",lm.1$coefficients[2],"\n")
       cat(" * Y-intercept from log(power) vs. log(frequency) fit (b):",lm.1$coefficients[1],"\n")
      }
# save linear fit, which is in log power
     lm.1.fit= as.vector(lm.1$fit)

# see Vaughan (2005, pg. 3) and Abramowitz & Stegun (1964): Here we want the 
#  expectation value of the log of the spectrum, which is not the expectation 
#  of log(spectrum). We must apply a bias correction for y = mx + b, 
#  and the power law relationship P=Nf^-beta
#  beta = - m 
     beta =-1*lm.1$coefficients[2]
# log(N) = b - bias
     bias = (digamma(dof/2)-log(dof/2))/log(10)
     logN = lm.1$coefficients[1] - bias

     if (verbose) 
      {
       cat(" * beta =",beta,"\n")
       cat(" * log(N) =",logN,"\n")
       cat(" * estimated bias =",bias,"\n")
      }
# this is the bias corrected power law fit, as log10(power)
     pwrLaw = logN - beta*log10(specfit[,1])
# convert PLfit from log10 power to power
     pwrLaw=10^pwrLaw
     chiPL <- (specfit[,2]/pwrLaw) * dof
     pwrLaw_CL <- pchisq(chiPL, df=dof)
     if(bc) pwrLaw_CL <- 1-pmin(1,(1-pwrLaw_CL)*nf1)
### 90, 95 and 99% confidence levels
     pwrLaw_90 <- pwrLaw*qchisq(probs[1], df=dof)/dof
     pwrLaw_95 <- pwrLaw*qchisq(probs[2], df=dof)/dof
     pwrLaw_99 <- pwrLaw*qchisq(probs[3], df=dof)/dof
     fft.out <- data.frame(cbind(fft.out[,1],fft.out[,2],fft.out[,3],fft.out[,4],pwrLaw_CL*100,pwrLaw,pwrLaw_90,pwrLaw_95,pwrLaw_99))
     colnames(fft.out) <- c('Frequency','Amplitude','Power','Phase','PwrLaw_CL','PwrLaw_Fit','PwrLaw_90_power','PwrLaw_95_power','PwrLaw_99_power')
    }

# Mann and Lees Robust AR(1) noise model spectrum (as modified in Patterson et al., 2014)
   if(background==3)
    {
      freqSmooth = Nyq * medsmooth
# number of frequencies per each smoothing window
      nptsSmooth = as.integer( round(freqSmooth/df, digits=1) )
# check to see if nptsSmooth is even
      if(nptsSmooth %% 2 == 0) 
       {
          nptsSmooth=nptsSmooth+1
          if(verbose) cat(" * Median smoothing window increased by 1 point\n")
       }
      if(verbose) cat(" * Number of smoothing points=",nptsSmooth,"\n")
      pwrMedian <- runmed(fft.out[,3],nptsSmooth,endrule="median")
# determine average power (So)
      So <- mean(pwrMedian)
# check to ensure logarthim exists
      if(min(pwrMedian) < 0) 
       {
         cat("**** ERROR IN ANALYTIC FIT TO MEDIAN SMOOTHER: Spectrum has negative power\n")
         stop("**** TERMINATING NOW!")
       }
# define function for Brent's method
       rednoise1 <- function (rho)
         {
           rhospec=So * (1-rho^2)/(1 - (2*rho*cos(pi*fft.out[,1]/Nyq)) + rho^2)
           sum((log(pwrMedian)-log(rhospec))^2)
         }
# npts is the length of data vector 'd', which may or may not have been
# detrended/demeaned, so ensure mean value is zero
       d0=d[,2]-mean(d[,2])
       rho_raw=sum(d0[1:(npts-1)] * d0[2:npts]) / sum(d0^2)
       if(verbose) cat(" * Calculating analytic fit of AR1 to median smoothed spectrum using Brent's method\n")
       rho=optim(par=rho_raw,rednoise1,method="Brent",lower=0,upper=1)$par    
       if(verbose) cat(" * Estimated Robust AR1 coefficient =",rho,"\n")

       ML96 = So * (1-(rho^2)) / (  1 - (2*rho*cos(pi*fft.out[,1]/Nyq)) + (rho^2) )
       chiML96 <-  (fft.out[,3]/ML96) * dof
       chiCL_ML96 <- pchisq(chiML96, df=dof)
       if(bc) chiCL_ML96 <- 1-pmin(1,(1-chiCL_ML96)*nf1)
# 90, 95 and 99% CL for power
       ML96_90 <- ML96*qchisq(probs[1], df=dof)/dof
       ML96_95 <- ML96*qchisq(probs[2], df=dof)/dof
       ML96_99 <- ML96*qchisq(probs[3], df=dof)/dof
       fft.out <- data.frame(cbind(fft.out[,1],fft.out[,2],fft.out[,3],fft.out[,4],chiCL_ML96*100,ML96,ML96_90,ML96_95,ML96_99))
       colnames(fft.out) <- c('Frequency','Amplitude','Power','Phase','ML96_CL','ML96_Fit','ML96_90_power','ML96_95_power','ML96_99_power')
     }
     

### output the Fourier coefficients if desired
   if (output==2)
     {
# note that all results are frequencies are output
      fc.out <- data.frame(cbind(freq,Re(ft),Im(ft)))
      colnames(fc.out)[1] <- 'Frequency'
      colnames(fc.out)[2] <- 'Real Coeff.'
      colnames(fc.out)[3] <- 'Imag. Coeff.'
     }

### output the noise coefficients if desired
# conventional AR1 and robust AR1
if (background == 1 || background == 3) 
      {
        ncoeffs.out=data.frame(cbind(rho,So))
        rownames(ncoeffs.out)[1] <- c("AR1_fit")
      }  
# power law fit
    if(background == 2) 
      {
        ncoeffs.out=data.frame(cbind(beta,logN,bias))
        rownames(ncoeffs.out)[1] <- c("PL_fit")   
      }    


   if(genplot && output != 2)
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
             lines(fft.out[,1],fft.out[,6],lwd=2)
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
             lines(fft.out[ii,1],log(fft.out[ii,6]),lwd=2)
             lines(fft.out[ii,1],log(fft.out[ii,7]),lwd=1,lty=3)
             lines(fft.out[ii,1],log(fft.out[ii,8]),lwd=1,lty=3)
             lines(fft.out[ii,1],log(fft.out[ii,9]),lwd=1,lty=3)
           }
        }
     plot(fft.out[,1],fft.out[,2], type="l",col="red", ylab="Amplitude", xlim=c(xmin,xmax), xlab="Frequency", main="Periodogram Amplitude",bty="n")
     plot(fft.out[,1],fft.out[,4], type="l",col="red", ylab="Phase", xlim=c(xmin,xmax), xlab="Frequency", main="Periodogram Phase",bty="n")
   }

   if(genplot && output == 2)
    {
      dev.new(title = "Fourier Coefficients", height = 7, width = 9)
      par(mfrow=c(2,1))
# sort frequencies for plotting (preserve original FFT order for output)
      fc.pl <- fc.out[order(fc.out[,1],na.last=NA,decreasing=FALSE),]
### plot the results
      plot(fc.pl[,1],fc.pl[,2],type="l", col="blue",ylab="Real Coefficient (a)", xlab="Frequency", main="Real Fourier Coefficients",bty="n")
      plot(fc.pl[,1],fc.pl[,3],type="l", col="blue",ylab="Imagiary Coefficient (b)", xlab="Frequency", main="Imaginary Fourier Coefficients",bty="n")
    }
   
if (output==1)  return(fft.out)
if (output==2)  return(fc.out)
if (output==3 && background != 0) return(ncoeffs.out)
if (output==4 && background != 0) return( list(fft.out,ncoeffs.out) )

 
#### END function periodogram
}