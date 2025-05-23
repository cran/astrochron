### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2025 Stephen R. Meyers
###
###########################################################################
### lowspec function - (SRM: September 28, 2012; October 8,11,15 2012; May 20-21, 2013; 
###                           May 23-25, 2013; May 27, 2013; June 5, 2013; June 13, 2013;
###                           June 21, 2013; Nov. 26. 2013; January 31, 2015; 
###                           February 1-3, 2015; February 26, 2015; March 5-6, 2015;
###                           September 10, 2015; August 21-22, 2016; March 20, 2017;
###                           June 12, 2017; November 19-20, 2017; August 21, 2018;
###                           January 14, 2021; January 4, 2025)
###
### "Seeing red" algorithm as an R function
###########################################################################

lowspec <- function (dat,decimate=NULL,tbw=3,padfac=5,detrend=F,siglevel=0.9,setrho=rho_raw_sig,lowspan=1,b_tun=-1,output=0,CLpwr=T,xmin=0,xmax=Nyq,pl=1,sigID=T,genplot=T,verbose=T)
{

if(verbose) cat("\n----- PERFORMING Robust Locally-Weighted Regression Spectral Background Estimation -----\n")

if(tbw !=3 && tbw !=2) {stop("**** ERROR: time-bandwidth product must be 2 or 3")}
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

if(!is.null(decimate) && decimate<dt) {stop("**** ERROR: Cannot add data, only remove")}
if(!is.null(decimate) && decimate==dt) {if(verbose) cat(" * Selected value for decimation is identical to data sampling interval. Will not decimate.\n")}
if(!is.null(decimate) && decimate>dt)
 {
   if(verbose) cat(" * Decimating stratigraphic series from sampling interval of",dt,"to",decimate,"\n")
   dat <- linterp(dat,dt=decimate,genplot=F,verbose=F)
   dt = decimate
   npts = length(dat[,1])
 }

if((npts*lowspan) <= 100) {stop("**** ERROR: The series has too few data points.  REQUIREMENT: (# points)*(lowspan) > 100")}


###########################################################################
### Prewhiten with Raw AR(1) filter
###########################################################################
### what is the estimated AR1 coefficient?
# npts is the length of data vector 'dat', which may or may not have been 
# detrended/demeaned, so ensure mean value is zero
    d0=dat[,2]-mean(dat[,2])
    rho_raw_sig=sum(d0[1:(npts-1)] * d0[2:npts]) / sum(d0^2)
    if(verbose) cat("\n * Estimated AR1 coefficient=",rho_raw_sig,"\n")

if (setrho != 0)
  {
    if(verbose) cat(" * Prewhitening with AR1 coefficient of",setrho,"\n")
### now prewhiten
    prewhite=dat[2:npts,2] - setrho * dat[1:(npts-1),2]
    nnpts=npts-1
### what is the new raw estimate for the AR(1) coefficient?
    d0=prewhite-mean(prewhite)
    coeff_est_prewhite=sum(d0[1:(nnpts-1)] * d0[2:nnpts]) / sum(d0^2)
    if(verbose) cat(" * Prewhitened AR1 coefficient=",coeff_est_prewhite,"\n")
   }
   
if (setrho == 0) 
   {
     prewhite <- dat[,2]
     nnpts=npts
   } 

### use least-squares fit to remove linear trend
if (detrend) 
  {
    lm.0 <- lm(prewhite ~ dat[1:nnpts,1])
    prewhite <- prewhite - (lm.0$coeff[2]*dat[1:nnpts,1] + lm.0$coeff[1])
    if(verbose) cat("\n * Linear trend removed. m=",lm.0$coeff[2],"b=",lm.0$coeff[1],"\n")
  }


###########################################################################
### MTM Power spectrum of prewhitend signal using 'multitaper' library
###########################################################################

### calculate Nyquist freq
Nyq <- 1/(2*dt)
### calculate Rayleigh frequency
Ray <- 1/(dt*nnpts)
### set number of DPSS sequences to use
numtap <- (2*tbw)-1
### number of points in padded series
npad=nnpts*padfac
### add another zero if we don't have an even number of data points, so Nyquist exists.   
if((npad*padfac)%%2 != 0) npad = npad + 1
# padded frequency grid
df = 1/(npad*dt)

if(verbose)
 {
   cat("\n * Nyquist frequency:",Nyq,"\n")
   cat(" * Rayleigh frequency:",Ray,"\n")
   cat(" * MTM Power spectrum bandwidth resolution (halfwidth):",tbw/(nnpts*dt),"\n")
   cat(" * Padded to",npad,"points\n")
 }

# make prewhite a time series object, here with unit sampling interval
prewhiteTS<-as.ts(prewhite)
specwhite <- spec.mtm(prewhiteTS,nw=tbw,k=numtap,Ftest=T,nFFT=npad,taper=c("dpss"),centre=c("Slepian"),jackknife=F,returnZeroFreq=F,plot=F)

# note: no zero frequency present, also remove Nyquist now
nfreq = length(specwhite$freq) - 1
freq <- specwhite$freq[1:nfreq]/dt
# normalize power (divided by npts in spec.mtm)
pwr <- specwhite$spec[1:nfreq]/nnpts
Ftest <- specwhite$mtm$Ftest[1:nfreq]


###########################################################################
### Estimate "smooth" continum using lowess with large span (=1) to
###   prewhitened signal
###########################################################################

### set tuning parameter for LOWESS fit
if(b_tun < 0)
 {
  if(tbw == 2) { tun = sqrt(nnpts)*0.35268 - 0.31158 }      
  if(tbw == 3) { tun = sqrt(nnpts)*0.245136 + 0.023780 }  
  if(verbose) cat(" * LOWSPEC tuning parameter =",tun,"\n")
 } 


if(b_tun >= 0) {tun = b_tun}
  
### fit a LOWESS smoother to power spectrum using function rfbaseline in package IDPmisc
   baseout <- rfbaseline(freq,pwr,span=lowspan,NoXP=NULL,maxit=c(2,2),b=tun,weight=NULL,Scale = function(r) median(abs(r))/0.6745,delta = NULL,SORT=FALSE,DOT = FALSE, init = NULL)
   y.predict <- as.vector(baseout$fit)

     
###########################################################################
### Determine confidence levels of measured spectrum versus smoothed
###   estimate
###########################################################################

### degrees of freedom
dof = (2*numtap)

### create plot of CL versus frequency 
### chi square value estimate = (pwr/background) * DOFs
chi <-  (pwr/y.predict) * dof

### calculate CL
chiCL <- pchisq(chi, df=dof)

### 90, 95 and 99% CL 
CL90 <- y.predict*qchisq(0.9, df=dof)/dof
CL95 <- y.predict*qchisq(0.95, df=dof)/dof
CL99 <- y.predict*qchisq(0.99, df=dof)/dof

### calculate f-test CL for prewhitened spectrum
prob <- pf(Ftest,2,dof-2)



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
res2 <- peakfilter(numpeak,nfreq,tbwRay,siglevel,freqloc,probmax,freq,y.predict,pwr,chiCL)
numpeak2=res2$nout
loc=res2$loc[1:numpeak2]
Frequency <- freq[loc]
Harmonic_CL <- prob[loc]
Red_Noise_CL <- chiCL[loc] 


if(verbose) 
  {
    cat(" * Number of significant F-test peaks identified =",numpeak2,"\n")
    cat("ID  / Frequency / Period / Harmonic_CL / LOWSPEC_CL\n")
    for(i in 1:numpeak2) cat(i," ", Frequency[i]," ",1/Frequency[i]," ",Harmonic_CL[i]*100," ",Red_Noise_CL[i]*100,"\n")
  }  


if(genplot)
  {
    if(pl == 1) logxy="y"
    if(pl == 2) logxy=""
    if(pl == 3) logxy="xy"
    if(pl == 4) logxy="x"
    if(pl == 3 || pl == 4) xmin=freq[1]
    par(mfrow=c(3,1))
    plot(freq,pwr, type="l",col="black",xlim=c(xmin,xmax),cex.axis=1.1,cex.lab=1.1,lwd=2,xlab="Frequency",ylab="Power",main="MTM Power and LOWSPEC Continuum Estimates",bty="n",log=logxy)
    lines(freq,y.predict, col="red", lwd=2)
    if(CLpwr) 
        {
              lines(freq,CL90,col="red",lwd=1,lty=3)
              lines(freq,CL95,col="red",lwd=1,lty=3)
              lines(freq,CL99,col="red",lwd=1,lty=3)
        }

### plot "significant" frequencies on power spectrum
    if(sigID && (numpeak2) > 0)
     {
      plfreq=double(numpeak2)
      pltext=double(numpeak2)
      for (k in 1:(numpeak2))
        {
           plfreq[k]=Frequency[k]
           pltext[k]=k
        }
      abline(v=plfreq,col="gray",lty=3)
      mtext(pltext[seq(from=1,to=(numpeak2),by=2)], side=3,line=0.25,at=plfreq[seq(from=1,to=(numpeak2),by=2)],cex=0.5,font=4)
      if((numpeak2) > 1) mtext(pltext[seq(from=2,to=(numpeak2),by=2)], side=3,line=-0.25,at=plfreq[seq(from=2,to=(numpeak2),by=2)],cex=0.5,font=4)
     }
     
    if(pl == 1) logxy=""
    if(pl == 2) logxy=""
    if(pl == 3) logxy="x"
    if(pl == 4) logxy="x"  
# plot LOWSPEC confidence levels
    plot(freq,chiCL*100,type="l",col="red",xlim=c(xmin,xmax),lwd=2,cex.axis=1.1,cex.lab=1.1,xlab="Frequency",ylab="Confidence Level",main="LOWSPEC Confidence Level Estimates",bty="n",log=logxy)
    abline(h=c(90,95,99),col="black",lty=3)

### plot "significant" frequencies on LOWSPEC confidence level spectrum
    if(sigID && (numpeak2) > 0)
     {
      abline(v=plfreq,col="gray",lty=3)
      mtext(pltext[seq(from=1,to=(numpeak2),by=2)], side=3,line=0.25,at=plfreq[seq(from=1,to=(numpeak2),by=2)],cex=0.5,font=4)
      if((numpeak2) > 1) mtext(pltext[seq(from=2,to=(numpeak2),by=2)], side=3,line=-0.25,at=plfreq[seq(from=2,to=(numpeak2),by=2)],cex=0.5,font=4)
     }

# plot f-test confidence levels
    plot(freq,prob*100,type="l",col="red",xlim=c(xmin,xmax),ylim=c(80,100),cex.axis=1.1,cex.lab=1.1,xlab="Frequency",ylab="Confidence Level",main="MTM Harmonic F-test Confidence Level Estimates",bty="n",lwd=2,log=logxy)
    abline(h=c(90,95,99),col="black",lty=3)

### plot "significant" frequencies on f-test condfidence levels
    if(sigID && (numpeak2) > 0)
     {
      abline(v=plfreq,col="gray",lty=3)
      mtext(pltext[seq(from=1,to=(numpeak2),by=2)], side=3,line=0.25,at=plfreq[seq(from=1,to=(numpeak2),by=2)],cex=0.5,font=4)
      if((numpeak2) > 1) mtext(pltext[seq(from=2,to=(numpeak2),by=2)], side=3,line=-0.25,at=plfreq[seq(from=2,to=(numpeak2),by=2)],cex=0.5,font=4)
     }  
  }


if (output==1) 
{
    spectrum <- data.frame(cbind(freq,pwr,prob*100,chiCL*100,y.predict,CL90,CL95,CL99))
    colnames(spectrum)[1] <- 'Frequency'
    colnames(spectrum)[2] <- 'Prewhite_power'
    colnames(spectrum)[3] <- 'Harmonic_CL'
    colnames(spectrum)[4] <- 'LOWSPEC_CL'
    colnames(spectrum)[5] <- 'LOWSPEC_back'
    colnames(spectrum)[6] <- 'LOWSPEC_90_power'
    colnames(spectrum)[7] <- 'LOWSPEC_95_power'
    colnames(spectrum)[8] <- 'LOWSPEC_99_power'
    return(spectrum)
}

if (output==2) 
{
    sigfreq <- data.frame(Frequency[1:numpeak2])
    colnames(sigfreq) <- 'Frequency'
    return(sigfreq)
}

#### END function lowspec
}
