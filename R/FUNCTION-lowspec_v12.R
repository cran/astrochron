### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2014 Stephen R. Meyers
###
###########################################################################
### lowspec function - (SRM: September 28, 2012; October 8,11,15 2012; May 20-21, 2013; 
###                           May 23-25, 2013; May 27, 2013; June 5, 2013; June 13, 2013;
###                           June 21, 2013; Nov. 26. 2013)
###
### "Seeing red" algorithm as an R function
###########################################################################

lowspec <- function (dat,decimate=NULL,tbw=3,padfac=5,detrend=F,siglevel=0.9,xmin=0,xmax=Nyq,setrho=rho_raw_sig,lowspan=1,b_tun=-1,output=0,sigID=T,pl=1,genplot=T,verbose=T)
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
       dat <- dat[order(dat[1], na.last = NA, decreasing = F), ]
       dt <- dat[2,1]-dat[1,1]
       npts <- length(dat[,1])
     }
   dtest <- dat[2:npts,1]-dat[1:(npts-1),1] 
   epsm=1e-10
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


### if generating PDF plots
if(output == 3) 
  {
     pdf(file="LOWSPEC-results.pdf")
     par(mfrow=c(3,1))
  }
if(output != 3 && genplot) par(mfrow=c(3,1))


###########################################################################
### Prewhiten with Raw AR(1) filter
###########################################################################
### what is the estimated AR1 coefficient?
    lag0 <- dat[1:(npts-1),2]
    lag1 <- dat[2:npts,2]
    rho_raw_sig <- cor(lag0,lag1)
    if(verbose) cat("\n * Raw AR1 =",rho_raw_sig,"\n")

if (setrho != 0)
  {
    if(verbose) cat(" * Prewhitening with AR1 coefficient of",setrho,"\n")
    prewhite <- 1:(npts-1)
### now prewhiten
    j=1
    for (i in 2:npts )
     { 
       prewhite[j]=dat[i,2]-setrho*dat[i-1,2] 
       j=j+1
     } 

### what is the new raw estimate for the AR(1) coefficient?
    lag0 <- prewhite[1:(npts-2)]
    lag1 <- prewhite[2:(npts-1)]
    coeff_est_prewhite <- cor(lag0,lag1)
    if(verbose) cat(" * Prewhitened AR1 =",coeff_est_prewhite,"\n")
    nnpts=npts-1
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
   cat(" * MTM Power spectrum bandwidth resolution:",tbw/(nnpts*dt),"\n")
 }

# make prewhite a time series object, here with unit sampling interval
prewhiteTS<-as.ts(prewhite)
specwhite <- spec.mtm(prewhiteTS,nw=tbw,k=numtap,Ftest=T,nFFT=npad,taper=c("dpss"),centre=c("Slepian"),jackknife=F,returnZeroFreq=F,plot=F)

### save frequency and power to freq and pwr
freq2 <- specwhite$freq/dt
pwr2 <- specwhite$spec
Ftest2 <- specwhite$mtm$Ftest


###########################################################################
### Estimate "smooth" continum using lowess with large span (=1) to
###   prewhitened signal
###########################################################################
### remove Nyquist for fitting, f(0) already removed
nfreq=length(freq2) - 1

tpwr <- as.matrix(specwhite$spec[1:nfreq])
### convert pwr to row matrix (necessary for baseline)
tpwr <- t(tpwr)

### set tuning parameter for LOWESS fit
if(b_tun < 0)
 {
  if(tbw == 2) { tun = sqrt(nnpts)*0.35268 - 0.31158 }      
  if(tbw == 3) { tun = sqrt(nnpts)*0.245136 + 0.023780 }  
  if(verbose) cat(" * LOWSPEC tuning parameter =",tun,"\n")
 } 


if(b_tun >= 0) {tun = b_tun}
  
### fit a LOWESS smoother to power spectrum using library baseline
   baseout <- baseline(tpwr, method = 'rfbaseline', span=lowspan, b=tun)
   y.predict <- as.vector(baseout@baseline)

     
###########################################################################
### Determine confidence levels of measured spectrum versus smoothed
###   estimate
###########################################################################

### degrees of freedom
dof = (2*numtap)

### create plot of CL versus frequency 
### chi square value estimate = (pwr/background) * DOFs
chi <-  (pwr2[1:nfreq]/y.predict) * dof

### calculate CL
chiCL <- pchisq(chi, df=dof)

### calculate f-test CL for prewhitened spectrum
prob <- pf(Ftest2[1:nfreq],2,dof-2)



###########################################################################
### Identify "significant" frequencies
###########################################################################
if(verbose) 
  {
    cat("\n * Searching for significant spectral peaks that satisfy",siglevel*100,"% CL\n")
    cat("     requirements outlined in Meyers (2012):\n") 
    cat("ID  / Frequency / Period / Harmonic_CL / Rednoise_CL\n") 
  }

### identify the harmonic f-test peaks
### initialize variables and arrays
probmax <- 1:nfreq
freqloc <- 1:nfreq
plateau <- 0
ij <- 1 
for ( j in 1:(nfreq-2) )
  {
# Check for plateaus
    if (prob[j] == prob[j+1] || prob[j+1] == prob [j+2] ) 
      {
         plateau=plateau + 1
         if(verbose) cat("**** WARNING: plateau detected at",freq2[j+1],"Probability=",prob[j+1],"\n")
      }
    if ( prob[j] < prob[j+1] && prob[j+1] > prob[j+2] )
      {
# save peak
         probmax[ij] <- prob[j+1]
         freqloc[ij] <- j+1
         ij=ij+1
       }  
  }
numpeak <- ij-1


Frequency<- double(numpeak)
ii = 1
for (j in 1:numpeak)
     {
      ifreq = freqloc[j]
      test = 0
      if(probmax[j] >= siglevel )
        {
###  examine power spectrum +/- tbw*Ray to see if there is power at required Red Noise CL;   
###  also require that F-test peak is on a power peak relative the local background as
###  estimated with the lowess smoother.
         for  (k in 1:nfreq)
           {
              lowband = freq2[ifreq] - (tbw*Ray)
              highband = freq2[ifreq] + (tbw*Ray)
              if(freq2[k] >= lowband && freq2[k] <= highband )
                {
                      if(chiCL[k] >= siglevel && pwr2[ifreq] >= y.predict[ifreq]) {test = 1}                
                }
           }
         if(test == 1)
           {      
                Frequency[ii] <- freq2[ifreq]
                Harmonic_CL <- prob[ifreq]
                Red_Noise_CL <- chiCL[ifreq] 
                if(verbose) cat(ii," ", Frequency[ii]," ",1/Frequency[ii]," ",Harmonic_CL," ",Red_Noise_CL,"\n")
                ii=ii+1
           }
        }  
     }


if (plateau > 0 && verbose) 
  {
   cat(" * Number of plateaus detected=",plateau,"\n")
   }


if(output ==3 || genplot)
  {
    if(pl == 1) 
       {
         plot(freq2[1:nfreq],log(pwr2[1:nfreq]), type="l",col="black",xlim=c(xmin,xmax),cex.axis=1.1,cex.lab=1.1,lwd=2,xlab="Frequency",ylab="Log Power",main="MTM Power and LOWSPEC Continuum Estimates",bty="n")
         lines(freq2[1:nfreq],log(y.predict), col="red", lwd=2)
        } 
    if(pl == 2) 
       {
         plot(freq2[1:nfreq],pwr2[1:nfreq], type="l",col="black",xlim=c(xmin,xmax),cex.axis=1.1,cex.lab=1.1,lwd=2,xlab="Frequency",ylab="Linear Power",main="MTM Power and LOWSPEC Continuum Estimates",bty="n")
         lines(freq2[1:nfreq],y.predict, col="red", lwd=2)
       }

### plot "significant" frequencies on power spectrum
    if(sigID)
     {
      setline=0.25
      for  (k in 1:(ii-1))
        {
           abline(v=Frequency[k],col="green",lty=3)
           mtext(k, side=3,line=setline,at=Frequency[k],cex=0.5,font=4)
           setline=setline*-1
        }
     }
     
# plot LOWSPEC confidence levels
    plot(freq2[1:nfreq],chiCL*100,type="l",col="red",xlim=c(xmin,xmax),lwd=2,cex.axis=1.1,cex.lab=1.1,xlab="Frequency",ylab="Confidence Level",main="LOWSPEC Confidence Level Estimates",bty="n")
    abline(h=90,col="black",lty=3)
    abline(h=95,col="black",lty=3)
    abline(h=99,col="black",lty=3)

### plot "significant" frequencies on LOWSPEC confidence level spectrum
    if(sigID)
     {
      setline=0.25
      for  (k in 1:(ii-1))
        {
           abline(v=Frequency[k],col="green",lty=3)
           mtext(k, side=3,line=setline,at=Frequency[k],cex=0.5,font=4)
           setline=setline*-1
        }
     }

# plot f-test confidence levels
    plot(freq2[1:nfreq],prob*100,type="l",col="red",xlim=c(xmin,xmax),ylim=c(80,100),cex.axis=1.1,cex.lab=1.1,xlab="Frequency",ylab="Confidence Level",main="MTM Harmonic Test Confidence Level Estimates",bty="n",lwd=2)
    abline(h=90,col="black",lty=3)
    abline(h=95,col="black",lty=3)
    abline(h=99,col="black",lty=3)

### plot "significant" frequencies on f-test condfidence levels
    if(sigID)
     {
      setline=0.25
      for  (k in 1:(ii-1))
        {
           abline(v=Frequency[k],col="green",lty=3)
           mtext(k, side=3,line=setline,at=Frequency[k],cex=0.5,font=4)
           setline=setline*-1
        }
     }  
  }


if (output==1) 
{

    spectrum <- data.frame(cbind(freq2[1:nfreq],pwr2[1:nfreq],y.predict[1:nfreq],chiCL[1:nfreq]*100,prob[1:nfreq]*100))
    colnames(spectrum)[1] <- 'Frequency'
    colnames(spectrum)[2] <- 'Prewhite_power'
    colnames(spectrum)[3] <- 'LOWSPEC_back'
    colnames(spectrum)[4] <- 'LOWSPEC_CL'
    colnames(spectrum)[5] <- 'F-test_CL'

  return(spectrum)
}

if (output==2) 
{
    sigfreq <- data.frame(Frequency[1:(ii-1)])
    colnames(sigfreq) <- 'Frequency'
    return(sigfreq)
}

if(output == 3) 
  {
    if(verbose) cat(" * LOWSPEC-results.pdf file generated and saved to working directory\n")
    dev.off()
   }

#### END function lowspec
}
