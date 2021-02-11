### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2021 Stephen R. Meyers
###
###########################################################################
### MTM function - (SRM: February 28, 2012; March 29, 2012; 
###                      September 18, 2012; Oct 8, 2012; Oct. 12, 2012; Nov. 23, 2012
###                      May 20-21, 2013; May 23, 2013; May 27, 2013; June 5, 2013; 
###                      June 13, 2013; August 5, 2013; August 12, 2013; Nov. 26, 2013;
###                      July 31, 2014; January 31, 2015; February 1-3, 2015; 
###                      February 26, 2015; March 6, 2015; June 30, 2015; Sept. 10, 2015;
###                      December 14, 2015; May 20, 2016; August 22, 2016; Oct. 4, 2016;
###                      March 20, 2017; November 20, 2017; August 17, 2018; 
###                      January 14, 2021)
###
### uses multitaper library and built in functions from R
###########################################################################

mtm <- function (dat,tbw=3,ntap=NULL,padfac=5,demean=T,detrend=F,siglevel=0.9,ar1=T,output=0,CLpwr=T,xmin=0,xmax=Nyq,pl=1,sigID=T,genplot=T,verbose=T)
{

if(verbose) cat("\n----- PERFORMING Multitaper Spectral Analysis -----\n")

dat <- data.frame(dat)
npts <- length(dat[,1])
dt <- dat[2,1]-dat[1,1]

# error checking 
   if(dt<0)
     { 
       if (verbose) cat("\n * Sorting data into increasing height/depth/time, removing empty entries\n")
       dat <- dat[order(dat[,1], na.last = NA, decreasing = F), ]
       dt <- dat[2,1]-dat[1,1]
       npts <- length(dat[,1])
     }
   dtest <- dat[2:npts,1]-dat[1:(npts-1),1] 
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

numtap=trunc((2*tbw)-1)
if(is.null(ntap)) 
 {
  ntap=numtap
  if (verbose) cat(" * Will use default setting of",ntap,"DPSS tapers\n")
 }

if(ntap > numtap)
 {
  ntap=numtap
  if (verbose) cat("**** WARNING: The number of DPSS tapers specified is too large. ntap reset to default value of",ntap,"\n")
 } 

if(ntap < 2)
 {
  ntap=numtap
  if (verbose) cat("**** WARNING: The number of DPSS tapers specified is too small. ntap reset to default value of",ntap,"\n")
 } 

  
###########################################################################
### MTM Power spectrum using 'multitaper' library
###########################################################################
### this version computes the adaptive multitaper spectrum

# remove mean and linear trend if requested
if (demean) 
  { 
    dave <- colMeans(dat[2])
    dat[2] <- dat[2] - dave
    if(verbose) cat(" * Mean value subtracted=",dave,"\n")
  }

if (!demean && verbose) cat(" * Mean value NOT subtracted\n")

### use least-squares fit to remove linear trend
if (detrend) 
  {
    lm.0 <- lm(dat[,2] ~ dat[,1])
    dat[2] <- dat[2] - (lm.0$coeff[2]*dat[1] + lm.0$coeff[1])
    if(verbose) cat(" * Linear trend subtracted. m=",lm.0$coeff[2],"b=",lm.0$coeff[1],"\n")
  }

    if (!detrend && verbose) cat(" * Linear trend NOT subtracted\n")

### calculate Nyquist freq
Nyq <- 1/(2*dt)
### calculate Rayleigh frequency
Ray <- 1/(dt*npts)
npad=as.integer(npts*padfac)
### add another zero if we don't have an even number of data points, so Nyquist exists.   
if((npad*padfac)%%2 != 0) npad = npad + 1
# padded frequency grid
df = 1/(npad*dt)

if(verbose)
 {
   cat(" * Nyquist frequency:",Nyq,"\n")
   cat(" * Rayleigh frequency:",Ray,"\n")
   cat(" * MTM Power spectrum bandwidth resolution (halfwidth):",tbw/(npts*dt),"\n")
   cat(" * Padded to",npad,"points\n")
 }

# make dat a time series object, here with unit sampling interval
datTS<-as.ts(dat[,2])
spec <- spec.mtm(datTS,nw=tbw,k=ntap,Ftest=T,nFFT=npad,taper=c("dpss"),centre=c("none"),jackknife=F,returnZeroFreq=F,plot=F)

# assign correct frequencies to spec$freq (note: this variable is returned if output = 4)
spec$freq <- spec$freq/dt
# note: no zero frequency present, also remove Nyquist now
nfreq = length(spec$freq) - 1
freq <- spec$freq[1:nfreq]

# normalize power (divided by npts in spec.mtm)
pwrRaw <- spec$spec[1:nfreq]/npts
FtestRaw <- spec$mtm$Ftest[1:nfreq]

# AR(1) noise model spectrum
if(ar1)
 {
### what is the estimated AR1 coefficient?
    lag0 <- dat[1:(npts-1),2]
    lag1 <- dat[2:npts,2]
    rho_raw <- cor(lag0,lag1)
### Calculate Raw red noise spectrum
### "So" is the average power. This can be determined from the white noise variance
###  as So = var/(1-rho^2), where rho is the lag-1 coeff
###  We will determine average power directly from measured spectrum
    So = mean(pwrRaw)
    RawAR = So * (1-(rho_raw^2)) / (  1 - (2*rho_raw*cos(pi*freq/Nyq)) + (rho_raw^2) )
    dofAR = (2*ntap)
    chiRawAR <-  (pwrRaw/RawAR) * dofAR
    chiCLRawAR <- pchisq(chiRawAR, df=dofAR)
### 90, 95 and 99% CL for power plot
    if(CLpwr)
     {
       AR1_90 <- RawAR*qchisq(0.9, df=dofAR)/dofAR
       AR1_95 <- RawAR*qchisq(0.95, df=dofAR)/dofAR
       AR1_99 <- RawAR*qchisq(0.99, df=dofAR)/dofAR
    }
 }

### f-test CL
dof=2*ntap
prob <- pf(FtestRaw,2,dof-2)


###########################################################################
### Identify "significant" frequencies
###########################################################################

if(verbose) 
  {
    cat("\n * Searching for significant spectral peaks that satisfy",siglevel*100,"% CL\n")
    cat("     requirements outlined in Meyers (2012):\n") 
  }

### identify the harmonic f-test peaks
res <- peak(cbind(freq,prob),level=siglevel,genplot=F,verbose=F)

numpeak = length(res[,1])
freqloc = res[,1]
probmax = res[,3]
# Frequency and Harmonic_CL will be written over later if ar1=T
Frequency <- freq[freqloc]
Harmonic_CL <- prob[freqloc]

if(!ar1 && verbose) 
  { 
    cat(" * Number of significant F-test peaks identified =",numpeak,"\n")
    cat("ID  / Frequency / Period / Harmonic_CL \n")
    for(i in 1:numpeak) cat(i," ", Frequency[i]," ",1/Frequency[i]," ",Harmonic_CL[i]*100,"\n")
  }  


# if ar1 selected, also evaluate red noise
if(ar1)
 {

# FORTRAN wrapper
peakfilter <- function (numpeak,nfreq,tbwRay,siglevel,freqloc,probmax,freq,background,pwr,cl) 
 { 
    F_dat = .Fortran('peakfilter_r',
    
    numpeak=as.integer(numpeak),nfreq=as.integer(nfreq),tbwRay=as.double(tbwRay),
    siglevel=as.double(siglevel),freqloc=as.integer(freqloc),probmax=as.double(probmax),
    freq=as.double(freq),background=as.double(background),pwr=as.double(pwr),
    cl=as.double(cl),
    
    loc=integer(as.integer(numpeak)), nout=integer(1)
    )    
# return the results
    return(F_dat)
 }

# identify maxima of peaks
tbwRay=tbw*Ray
res2 <- peakfilter(numpeak,nfreq,tbwRay,siglevel,freqloc,probmax,freq,RawAR,pwrRaw,chiCLRawAR)
numpeak2=res2$nout
loc=res2$loc[1:numpeak2]
Frequency <- freq[loc]
Harmonic_CL <- prob[loc]
Red_Noise_CL <- chiCLRawAR[loc] 

if(verbose) 
  {
    cat(" * Number of significant F-test peaks identified =",numpeak2,"\n")
    cat("ID  / Frequency / Period / Harmonic_CL / Rednoise_CL\n")
    for(i in 1:numpeak2) cat(i," ", Frequency[i]," ",1/Frequency[i]," ",Harmonic_CL[i]*100," ",Red_Noise_CL[i]*100,"\n")
  }  
# reassign numpeak  
numpeak=numpeak2
# end ar1=T section
 }



### generate plots
if(genplot)
 {

# first plot power spectrum, with red noise model and confidence levels if requested
   if(ar1) 
    {
      par(mfrow=c(3,1))
      if(!CLpwr) mtitle=c("MTM Power (black) & AR1 fit (red)")
      if(CLpwr) mtitle=c("MTM Power (black); AR1 fit (red); 90%CL, 95%CL, 99%CL (dotted)")
    }   
   if(!ar1) 
    {
     par(mfrow=c(2,1))
     mtitle=c("MTM Power")
    }
   if(pl == 1) logxy="y"
   if(pl == 2) logxy=""
   if(pl == 3) logxy="xy"
   if(pl == 4) logxy="x"
   if(pl == 3 || pl == 4) xmin=freq[1]
   plot(freq,pwrRaw,type="l", col="black", xlim=c(xmin,xmax), xlab="Frequency",ylab="Power",main=mtitle,cex.axis=1.1,cex.lab=1.1,lwd=2,bty="n",log=logxy)
      if(ar1) 
        {
          lines(freq,RawAR,col="red",lwd=2)
          if(CLpwr) 
            {
              lines(freq,AR1_90,col="red",lwd=1,lty=3)
              lines(freq,AR1_95,col="red",lwd=1,lty=3)
              lines(freq,AR1_99,col="red",lwd=1,lty=3)
             }
         }     

### plot "significant" frequencies on power spectrum
    if(sigID && (numpeak) > 0)
     {
      plfreq=double(numpeak)
      pltext=double(numpeak)
      for (k in 1:(numpeak))
        {
           plfreq[k]=Frequency[k]
           pltext[k]=k
        }
      abline(v=plfreq,col="gray",lty=3)
      mtext(pltext[seq(from=1,to=(numpeak),by=2)], side=3,line=0.25,at=plfreq[seq(from=1,to=(numpeak),by=2)],cex=0.5,font=4)
      if((numpeak) > 1) mtext(pltext[seq(from=2,to=(numpeak),by=2)], side=3,line=-0.25,at=plfreq[seq(from=2,to=(numpeak),by=2)],cex=0.5,font=4)
    }

# next plot harmonic F-test results or AR1 confidence levels
   if(pl == 1) logxy=""
   if(pl == 2) logxy=""
   if(pl == 3) logxy="x"
   if(pl == 4) logxy="x"   

  if(!ar1) 
   {
     plot(freq,prob*100,type="l",col="red",xlim=c(xmin,xmax),ylim=c(80,100),cex.axis=1.1,cex.lab=1.1,xlab="Frequency",ylab="Confidence Level",main="MTM Harmonic F-Test Confidence Level Estimates",bty="n",lwd=2,log=logxy)
     abline(h=c(90,95,99),col="black",lty=3)
     if(sigID && numpeak > 0)
      {
        abline(v=plfreq,col="gray",lty=3)
        mtext(pltext[seq(from=1,to=numpeak,by=2)], side=3,line=0.25,at=plfreq[seq(from=1,to=numpeak,by=2)],cex=0.5,font=4)
        if(numpeak > 1) mtext(pltext[seq(from=2,to=numpeak,by=2)], side=3,line=-0.25,at=plfreq[seq(from=2,to=numpeak,by=2)],cex=0.5,font=4)
      }
   }  
  if(ar1) 
   { 
     plot(freq,chiCLRawAR*100,type="l",col="red",xlim=c(xmin,xmax),ylim=c(0,100),cex.axis=1.1,cex.lab=1.1,lwd=2,xlab="Frequency",ylab="Confidence Level",main="AR1 Confidence Level Estimates",bty="n",log=logxy)
     abline(h=c(90,95,99),col="black",lty=3)
     if(sigID && numpeak > 0)
      {
        abline(v=plfreq,col="gray",lty=3)
        mtext(pltext[seq(from=1,to=numpeak,by=2)], side=3,line=0.25,at=plfreq[seq(from=1,to=numpeak,by=2)],cex=0.5,font=4)
        if(numpeak > 1) mtext(pltext[seq(from=2,to=numpeak,by=2)], side=3,line=-0.25,at=plfreq[seq(from=2,to=numpeak,by=2)],cex=0.5,font=4)
      }
     plot(freq,prob*100,type="l",col="red",xlim=c(xmin,xmax),ylim=c(80,100),cex.axis=1.1,cex.lab=1.1,xlab="Frequency",ylab="Confidence Level",main="Harmonic F-Test Confidence Level Estimates",bty="n",lwd=2,log=logxy)
     abline(h=c(90,95,99),col="black",lty=3)
     if(sigID && numpeak > 0)
      {
        abline(v=plfreq,col="gray",lty=3)
        mtext(pltext[seq(from=1,to=numpeak,by=2)], side=3,line=0.25,at=plfreq[seq(from=1,to=numpeak,by=2)],cex=0.5,font=4)
        if(numpeak > 1) mtext(pltext[seq(from=2,to=numpeak,by=2)], side=3,line=-0.25,at=plfreq[seq(from=2,to=numpeak,by=2)],cex=0.5,font=4)
      }
   } 
# end genplot section
 }
  

if (output==1) 
 {
   if(!ar1) 
     { 
       spectrum <- data.frame(cbind(freq,pwrRaw,prob*100))
       colnames(spectrum)[1] <- 'Frequency'
       colnames(spectrum)[2] <- 'Power'
       colnames(spectrum)[3] <- 'Harmonic_CL'
      }
   if(ar1)
     {
       spectrum <- data.frame(cbind(freq,pwrRaw,prob*100,chiCLRawAR*100,RawAR,AR1_90,AR1_95,AR1_99))
       colnames(spectrum)[1] <- 'Frequency'
       colnames(spectrum)[2] <- 'Power'
       colnames(spectrum)[3] <- 'Harmonic_CL'
       colnames(spectrum)[4] <- 'AR1_CL'
       colnames(spectrum)[5] <- 'AR1_fit'
       colnames(spectrum)[6] <- 'AR1_90_power'
       colnames(spectrum)[7] <- 'AR1_95_power'
       colnames(spectrum)[8] <- 'AR1_99_power'
     }  
   return(spectrum)
 }

if (output==2) 
 {
    sigfreq <- data.frame(Frequency[1:numpeak])
    colnames(sigfreq) <- 'Frequency'
    return(sigfreq)
 }

if (output==3) 
 {
    sigfreq <- data.frame(Frequency[1:numpeak],Harmonic_CL[1:numpeak])
    colnames(sigfreq)[1] <- 'Frequency'
    colnames(sigfreq)[2] <- 'Harmonic_CL'
    return(sigfreq)
 }

if (output==4) 
 {
   return(spec)
 }

#### END function mtm
}
