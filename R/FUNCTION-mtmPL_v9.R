### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2021 Stephen R. Meyers
###
###########################################################################
### MTM with power law fit and confidence levels - (SRM: Ocotober 29, 2015;
###                                    November 1, 2015; November 15-29, 2017;
###                                    December 4, 2017; April 13, 2018;
###                                    July 8, 2018; August 21, 2018; January 14, 2021)
###
### uses multitaper library and built in functions from R
### this follows the recipe of Vaughan et al. (2005):
###  (1) fit log10 power/ log10 frequency by least squares
###  (2) estimate alpha and N from fit
###  (3) bias correction of N
###  (4) estimate confidence levels with Chi-sq distribution
###
### we also need to avoid the lowest frequences in the fit, which are biased for MTM.
###########################################################################


mtmPL <- function (dat,tbw=3,ntap=NULL,padfac=5,demean=T,detrend=F,siglevel=0.9,flow=NULL,fhigh=NULL,output=0,CLpwr=T,xmin=0,xmax=Nyq,pl=1,sigID=F,genplot=T,verbose=T)
{

if(verbose) cat("\n----- PERFORMING Multitaper Spectral Analysis with 1/f noise test -----\n")

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
   cat(" * MTM Power spectrum bandwidth resolution:",tbw/(npts*dt),"\n")
   cat(" * Padded to",npad,"points\n")
 }

# make dat a time series object, here with unit sampling interval
datTS<-as.ts(dat[,2])
spec <- spec.mtm(datTS,nw=tbw,k=ntap,Ftest=T,nFFT=npad,taper=c("dpss"),centre=c("none"),jackknife=F,returnZeroFreq=F,plot=F)

# note: no zero frequency present, also remove Nyquist now
nfreq = length(spec$freq) - 1
freq <- spec$freq[1:nfreq]/dt

# by default, avoid lowermost frequencies <= mtm-halfwidth, as they are 
#   biased (see McCoy et al., 1998 and Huybers & Curry, 2006).
#if(is.null(flow)) flow=freq[which(freq==(tbw/(npts*dt)))+1]
# modified to allow for round-off error
if(is.null(flow)) flow=freq[as.integer(tbw/(npts*dt)/df) + 1]

if(is.null(fhigh)) fhigh=freq[nfreq]

# remove values that are not in range from flow to fhigh
ii=which((freq >= flow) & (freq <= fhigh))
nfreq = length(spec$freq[ii])
freq <- spec$freq[ii]/dt
# normalize power (divided by npts in spec.mtm)
pwrRaw <- spec$spec[ii]/npts
FtestRaw <- spec$mtm$Ftest[ii]

# power law fit and confidence level estimation
dof = (2*ntap)
specIn=data.frame(cbind(freq,pwrRaw))
resFit=pwrLawFit(specIn,dof=dof,flow=flow,fhigh=fhigh,output=1,genplot=F,verbose=verbose)

resFreq=resFit[,1]
unbiasedFit=resFit[,4]
CL_90=resFit[,5]
CL_95=resFit[,6]
CL_99=resFit[,7]
chiCL <- resFit[,3]/100

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
# Frequency will be rewritten over below
Frequency <- freq[freqloc]
Harmonic_CL <- prob[freqloc]

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
res2 <- peakfilter(numpeak,nfreq,tbwRay,siglevel,freqloc,probmax,freq,unbiasedFit,pwrRaw,chiCL)
numpeak2=res2$nout
loc=res2$loc[1:numpeak2]
Frequency <- freq[loc]
Harmonic_CL <- prob[loc]
Red_Noise_CL <- chiCL[loc] 

if(verbose) 
  {
    cat(" * Number of significant F-test peaks identified =",numpeak2,"\n")
    cat("ID  / Frequency / Period / Harmonic_CL / PowerLaw_CL\n")
    for(i in 1:numpeak2) cat(i," ", Frequency[i]," ",1/Frequency[i]," ",Harmonic_CL[i]*100," ",Red_Noise_CL[i]*100,"\n")
  }  
# reassign numpeak  
numpeak=numpeak2


### generate plots
if(genplot)
 {
   par(mfrow=c(3,1))
   if(!CLpwr) mtitle=c("1/f fit (blue), unbiased 1/f fit (red)")
   if(CLpwr) mtitle=c("1/f fit (blue), unbiased 1/f fit (red), 90%CL, 95%CL, 99%CL (dotted)")
   if(pl == 1) logxy="y"
   if(pl == 2) logxy=""
   if(pl == 3) logxy="xy"
   if(pl == 4) logxy="x"
   if(pl == 3 || pl == 4) xmin=freq[1]

# first plot power spectrum, with red noise model and confidence levels
   plot(freq,pwrRaw,type="l", col="black", xlim=c(xmin,xmax), xlab="Frequency",ylab="Power",main=mtitle,cex.axis=1.1,cex.lab=1.1,lwd=2,bty="n",log=logxy)
   lines(resFreq,unbiasedFit,col="red",lwd=2)
   if(CLpwr) 
        {
              lines(resFreq,CL_90,col="red",lwd=1,lty=3)
              lines(resFreq,CL_95,col="red",lwd=1,lty=3)
              lines(resFreq,CL_99,col="red",lwd=1,lty=3)
        }
 
### plot "significant" frequencies on power spectrum
   if(sigID && (numpeak) > 0)
    {
### plot "significant" F-test frequencies (on power plot first)
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

### plot powerLaw confidence levels
   if(pl == 1) logxy=""
   if(pl == 2) logxy=""
   if(pl == 3) logxy="x"
   if(pl == 4) logxy="x"   

   plot(resFreq,chiCL*100,type="l",col="red",xlim=c(xmin,xmax),ylim=c(0,100),cex.axis=1.1,cex.lab=1.1,lwd=2,xlab="Frequency",ylab="Confidence Level",main="1/f Confidence Level Estimates",bty="n",log=logxy)
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
# end genplot section  
   }


if (output==1) 
 {
   probfit <- data.frame(cbind(freq,prob))
   probfit <- subset(probfit, (probfit[1] >= flow) & (probfit[1] <= fhigh) )
   specfit <- subset(specIn, (specIn[1] >= flow) & (specIn[1] <= fhigh) )
   spectrum <- data.frame(cbind(resFreq,specfit[,2],probfit[,2]*100,chiCL*100,unbiasedFit,CL_90,CL_95,CL_99))
   colnames(spectrum)[1] <- 'Frequency'
   colnames(spectrum)[2] <- 'Power'
   colnames(spectrum)[3] <- 'Harmonic_CL'
   colnames(spectrum)[4] <- 'PowerLaw_CL'
   colnames(spectrum)[5] <- 'PowerLaw_fit'
   colnames(spectrum)[6] <- 'PowerLaw_90_power'
   colnames(spectrum)[7] <- 'PowerLaw_95_power'
   colnames(spectrum)[8] <- 'PowerLaw_99_power'
   return(spectrum)
 }


if (output==2) 
 {
    sigfreq <- data.frame(Frequency[1:numpeak2])
    colnames(sigfreq) <- 'Frequency'
    return(sigfreq)
 }

if (output==3) 
 {
    sigfreq <- data.frame(Frequency[1:numpeak2],Harmonic_CL[1:numpeak2])
    colnames(sigfreq)[1] <- 'Frequency'
    colnames(sigfreq)[2] <- 'Harmonic_CL'
    return(sigfreq)
 }


if (output==4) 
 {
   return(spec)
 }

#### END function mtmPL
}
