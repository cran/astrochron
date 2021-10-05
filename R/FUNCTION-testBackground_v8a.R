### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2018 Stephen R. Meyers
###
##############################################################################
### testBackground function - (SRM: November 19-29, 2017; December 4-5, 2017
###                                 October 4, 2021)
###
### Evaluate false positive rate and frequency distribution for
###  a range of different noise types and background estimation
###  approaches.
###########################################################################


testBackground <- function (npts=1001,dt=5,noiseType="ar1",coeff=NULL,method="periodogramAR1",opt=NULL,demean=T,detrend=F,low=0,tbw=3,multi=F,iter=2000,output=F,genplot=F,verbose=T)
{

###########################################################################
### 1. initialize all variables/vectors
###########################################################################
# padding
padfac=1

### calculate Nyquist freq
Nyq <- 1/(2*dt)
### calculate rayleigh frequency
Ray <- 1/(dt*npts)

### set the random number seed to allow duplication of simulations
#set.seed(403)
#set.seed(400)

# frac90, frac95 and frac99 record the fraction of the frequencies in 
#  each spectrum (iteration) that achieve the specified confidence level
frac90 <- double(iter)
frac95 <- double(iter)
frac99 <- double(iter)

# sum90, sum95 and sum99 record the total number of false positives at each 
#  frequency, across all iterations. it is dimensioned larger than needed.
sum90 <- double(npts)
sum95 <- double(npts)
sum99 <- double(npts)

# tick90all, tick95all and tick99all record the total number of false positives
#  across all frequencies and all iterations
tick90all = 0
tick95all = 0
tick99all = 0

# these record the total number of false positives across all 
#  frequencies and all iterations
# FDR-BH
fdrBH90all = 0
fdrBH95all = 0
fdrBH99all = 0
# FDR-BY
fdrBY90all = 0
fdrBY95all = 0
fdrBY99all = 0
# Hommel
hommel90all = 0
hommel95all = 0
hommel99all = 0
# Hochberg
hochberg90all = 0
hochberg95all = 0
hochberg99all = 0
# Holm
holm90all = 0
holm95all = 0
holm99all = 0
# Bonferroni
bonferroni90all = 0
bonferroni95all = 0
bonferroni99all = 0

# these record the number of spectra (iterations) that have false detections 
#  (regardless of number per spectrum).
# FDR-BH
fdrBH90any = 0
fdrBH95any = 0
fdrBH99any = 0
# FDR-BY
fdrBY90any = 0
fdrBY95any = 0
fdrBY99any = 0
# Hommel
hommel90any = 0
hommel95any = 0
hommel99any = 0
# Hochberg
hochberg90any = 0
hochberg95any = 0
hochberg99any = 0
# Holm
holm90any = 0
holm95any = 0
holm99any = 0
# Bonferroni
bonferroni90any = 0
bonferroni95any = 0
bonferroni99any = 0

if(noiseType=="ar1" && is.null(coeff)) coeff=0.9
if(noiseType=="pwrLaw" && is.null(coeff)) coeff=2

if (method=="mtmML96" && is.null(opt)) opt=0.2
if (method=="lowspec" && is.null(opt)) opt=1
if (method=="periodogramPL" && is.null(opt)) opt=0.25
if (method=="periodogramAR1" && is.null(opt)) opt=0.25

if(verbose) 
  {  
    cat("\n SURROGATES:\n")
    if(noiseType=="ar1") cat(" Using AR1 surrogates (function ar1) with rho of", coeff,"\n\n")
    if(noiseType=="pwrLaw") cat(" Using Power Law surrogates (function pwrLaw) with beta of",coeff,"\n\n")
    
    cat(" SPECTRAL APPROACH:\n")
    if(method=="mtmAR1") cat(" Using MTM with Conventional AR1 background fit (function mtm)\n")
    if(method=="mtmML96") cat(" Using MTM with ML96 AR1 background fit (function mtmML96): medsmooth =",opt,"\n")
    if(method=="lowspec") cat(" Using MTM LOWSPEC background fit (function lowspec), with lowspan =",opt, "\n")
    if(method=="mtmPL") cat(" Using MTM with Power Law background fit (function mtmPL)\n")
    if(method=="periodogramPL" && opt>0) cat(" Using Periodogram with Power Law background fit (function periodogram):",100*opt,"% cosine taper \n")
    if(method=="periodogramPL" && opt==0) cat(" Using Periodogram with Power Law background fit (function periodogram) \n")
    if(method=="periodogramAR1" && opt>0) cat(" Using Periodogram with Conventional AR1 background fit (function periodogram):",100*opt,"% cosine taper \n")
    if(method=="periodogramAR1" && opt==0) cat(" Using Periodogram with Conventional AR1 background fit (function periodogram) \n")


    cat("\n * PLEASE WAIT: Performing Simulations\n")
    cat("\n0%       25%       50%       75%       100%\n")
# create a progress bar
    progress = utils::txtProgressBar(min = 0, max = iter, style = 1, width=43)
  }

#### start loop
for (ii in 1:iter)
{
 if(verbose) utils::setTxtProgressBar(progress, ii)

###########################################################################
### 2. Make noise 
###########################################################################
# choose ar1 or power law noise.
shuffle=F
if(noiseType=="ar1") noise=ar1(npts=npts,dt=dt,rho=coeff,shuffle=shuffle,genplot=F,verbose=F)
if(noiseType=="pwrLaw") noise=pwrLaw(npts=npts,dt=dt,beta=coeff,genplot=F,verbose=F)
if(low > 0) noise=noLow(noise,smooth=low,genplot=F,verbose=F)

#noise=prewhiteAR1(noise,genplot=F,verbose=F)

###########################################################################
### 3. estimate spectrum
###########################################################################
####  (1) MTM-AR1
####  (2) MTM-ML96
####  (3) MTM-lowspec
####  (4) MTM-PL
####  (5) Periodogram-PL
####  (6) Periodogram-AR1
if (method=="mtmAR1") spec <- mtm(noise,tbw=tbw,padfac=padfac,demean=demean,detrend=detrend,output=1,verbose=F,genplot=F)
if (method=="mtmML96") spec <- mtmML96(noise,tbw=tbw,medsmooth=opt,padfac=padfac,demean=demean,detrend=detrend,output=1,verbose=F,genplot=F)
# no demean option for lowspec  
if (method=="lowspec") spec <- lowspec(noise,tbw=tbw,lowspan=opt,padfac=padfac,detrend=detrend,output=1,verbose=F,genplot=F)
if (method=="mtmPL") spec <- mtmPL(noise,tbw=tbw,padfac=padfac,demean=demean,detrend=detrend,output=1,verbose=F,genplot=F)

if (method == "mtmAR1" || method == "mtmML96" || method == "lowspec" || method == "mtmPL")
 {
   freq <- spec[,1]
   pwr <- spec[,2]
   ifreq <- length(freq)
   smoothCL <- spec[,4]
   smooth <- spec[,5]
   smooth90 <- spec[,6]
   smooth95 <- spec[,7]
   smooth99 <- spec[,8]
 }

### periodogram power law fit (25% cos taper)
if (method=="periodogramPL")
 {
   if(opt > 0) 
     { 
       noise=cosTaper(noise,demean=demean,detrend=detrend,verbose=F,genplot=F)
       spec <- periodogram(noise,padfac=padfac,demean=F,detrend=F,f0=F,fNyq=F,background=2,output=1,verbose=F,genplot=F)
     }  
   if(opt == 0) spec <- periodogram(noise,padfac=padfac,demean=demean,detrend=detrend,f0=F,fNyq=F,background=2,output=1,verbose=F,genplot=F)
 
   freq <- spec[,1]
   pwr <- spec[,3]
   ifreq <- length(freq)
   smooth <- spec[,6]
   smoothCL <- spec[,5]
   smooth90 <- spec[,7]
   smooth95 <- spec[,8]
   smooth99 <- spec[,9]
 }

if (method=="periodogramAR1")
 {
   if(opt > 0) 
     { 
       noise=cosTaper(noise,demean=demean,detrend=detrend,verbose=F,genplot=F)
       spec <- periodogram(noise,padfac=padfac,demean=F,detrend=F,f0=F,fNyq=F,background=1,output=1,verbose=F,genplot=F)
     }  
   if(opt == 0) spec <- periodogram(noise,padfac=padfac,demean=demean,detrend=detrend,f0=F,fNyq=F,background=1,output=1,verbose=F,genplot=F)

   freq <- spec[,1]
   pwr <- spec[,3]
   ifreq <- length(freq)
   smooth <- spec[,6]
   smoothCL <- spec[,5]
   smooth90 <- spec[,7]
   smooth95 <- spec[,8]
   smooth99 <- spec[,9]
 }
 
###########################################################################
### 4. estimate false positive rates and multiple correction procedures
###########################################################################
### initialize tickers
 tick90 = 0
 tick95 = 0
 tick99 = 0
  
# speed this up by vectorization! 
### determine what fraction of spectrum exceeds 90%, 95% and 99% Confidence Levels
 for ( j in 1:ifreq )
  {
   if (pwr[j] > smooth90[j]) { tick90 = tick90 + 1 ; sum90[j] = sum90[j] + 1}
   if (pwr[j] > smooth95[j]) { tick95 = tick95 + 1 ; sum95[j] = sum95[j] + 1}
   if (pwr[j] > smooth99[j]) { tick99 = tick99 + 1 ; sum99[j] = sum99[j] + 1}
  }
   frac90[ii]=100*tick90/ifreq
   frac95[ii]=100*tick95/ifreq
   frac99[ii]=100*tick99/ifreq

   tick90all = tick90all + tick90
   tick95all = tick95all + tick95
   tick99all = tick99all + tick99

if(multi)
 {
  multiRes=multiTest(cbind(freq,smoothCL),genplot=F,verbose=F)
# FDR-BH
  fdrBH90all= fdrBH90all + sum(multiRes[,4]<.1)
  fdrBH95all= fdrBH95all + sum(multiRes[,4]<.05)
  fdrBH99all= fdrBH99all + sum(multiRes[,4]<.01)
  fdrBH90any= fdrBH90any + sum(any(multiRes[,4]<.1))
  fdrBH95any= fdrBH95any + sum(any(multiRes[,4]<.05))
  fdrBH99any= fdrBH99any + sum(any(multiRes[,4]<.01))
# FDR-BY
  fdrBY90all= fdrBY90all + sum(multiRes[,5]<.1)
  fdrBY95all= fdrBY95all + sum(multiRes[,5]<.05)
  fdrBY99all= fdrBY99all + sum(multiRes[,5]<.01)
  fdrBY90any= fdrBY90any + sum(any(multiRes[,5]<.1))
  fdrBY95any= fdrBY95any + sum(any(multiRes[,5]<.05))
  fdrBY99any= fdrBY99any + sum(any(multiRes[,5]<.01))
# Hommel
  hommel90all= hommel90all + sum(multiRes[,6]<.1)
  hommel95all= hommel95all + sum(multiRes[,6]<.05)
  hommel99all= hommel99all + sum(multiRes[,6]<.01)
  hommel90any= hommel90any + sum(any(multiRes[,6]<.1))
  hommel95any= hommel95any + sum(any(multiRes[,6]<.05))
  hommel99any= hommel99any + sum(any(multiRes[,6]<.01))
# Hochberg
  hochberg90all= hochberg90all + sum(multiRes[,7]<.1)
  hochberg95all= hochberg95all + sum(multiRes[,7]<.05)
  hochberg99all= hochberg99all + sum(multiRes[,7]<.01)
  hochberg90any= hochberg90any + sum(any(multiRes[,7]<.1))
  hochberg95any= hochberg95any + sum(any(multiRes[,7]<.05))
  hochberg99any= hochberg99any + sum(any(multiRes[,7]<.01))
# Holm
  holm90all= holm90all + sum(multiRes[,8]<.1)
  holm95all= holm95all + sum(multiRes[,8]<.05)
  holm99all= holm99all + sum(multiRes[,8]<.01)
  holm90any= holm90any + sum(any(multiRes[,8]<.1))
  holm95any= holm95any + sum(any(multiRes[,8]<.05))
  holm99any= holm99any + sum(any(multiRes[,8]<.01))
# Bonferroni
  bonferroni90all= bonferroni90all + sum(multiRes[,9]<.1)
  bonferroni95all= bonferroni95all + sum(multiRes[,9]<.05)
  bonferroni99all= bonferroni99all + sum(multiRes[,9]<.01)
  bonferroni90any= bonferroni90any + sum(any(multiRes[,9]<.1))
  bonferroni95any= bonferroni95any + sum(any(multiRes[,9]<.05))
  bonferroni99any= bonferroni99any + sum(any(multiRes[,9]<.01))
 }
   
if(genplot)
 {
  par(mfrow=c(1,1))
  plot(freq,log(pwr),col="red",type="l",main="Monte Carlo Simulation",ylab="Log(power)",xlab="Frequency(cycles/ka)")
  lines(freq,log(smooth))
 }
 
### end iteration loop
}

# close progress bar
if(verbose) close(progress) 

###########################################################################
### 5. Save results and plot
###########################################################################

tick90all=100*tick90all/(ifreq*iter)
tick95all=100*tick95all/(ifreq*iter)
tick99all=100*tick99all/(ifreq*iter)

med99=median(frac99)
med95=median(frac95)
med90=median(frac90)

if(verbose)
 {
   cat("\n STANDARD CONFIDENCE LEVEL SIMULATION RESULTS:\n")
   cat("   Median estimate 90% CL =", med90,"\n")
   cat("   Median estimate 95% CL =", med95,"\n")
   cat("   Median estimate 99% CL =", med99,"\n\n")
   cat("   Total percent > 90% CL =", tick90all,"\n")
   cat("   Total percent > 95% CL =", tick95all,"\n")
   cat("   Total percent > 99% CL =", tick99all,"\n")
 }

if(multi)
 {
# FDR-BH
   fdrBH90all = 100*fdrBH90all/iter
   fdrBH95all = 100*fdrBH95all/iter
   fdrBH99all = 100*fdrBH99all/iter
   fdrBH90any = 100*fdrBH90any/iter
   fdrBH95any = 100*fdrBH95any/iter
   fdrBH99any = 100*fdrBH99any/iter
# FDR-BY
   fdrBY90all = 100*fdrBY90all/iter
   fdrBY95all = 100*fdrBY95all/iter
   fdrBY99all = 100*fdrBY99all/iter
   fdrBY90any = 100*fdrBY90any/iter
   fdrBY95any = 100*fdrBY95any/iter
   fdrBY99any = 100*fdrBY99any/iter
# Hommel
   hommel90all = 100*hommel90all/iter
   hommel95all = 100*hommel95all/iter
   hommel99all = 100*hommel99all/iter
   hommel90any = 100*hommel90any/iter
   hommel95any = 100*hommel95any/iter
   hommel99any = 100*hommel99any/iter
# Hochberg
   hochberg90all = 100*hochberg90all/iter
   hochberg95all = 100*hochberg95all/iter
   hochberg99all = 100*hochberg99all/iter
   hochberg90any = 100*hochberg90any/iter
   hochberg95any = 100*hochberg95any/iter
   hochberg99any = 100*hochberg99any/iter
# Holm
   holm90all = 100*holm90all/iter
   holm95all = 100*holm95all/iter
   holm99all = 100*holm99all/iter
   holm90any = 100*holm90any/iter
   holm95any = 100*holm95any/iter
   holm99any = 100*holm99any/iter
# Bonferroni
   bonferroni90all = 100*bonferroni90all/iter
   bonferroni95all = 100*bonferroni95all/iter
   bonferroni99all = 100*bonferroni99all/iter
   bonferroni90any = 100*bonferroni90any/iter
   bonferroni95any = 100*bonferroni95any/iter
   bonferroni99any = 100*bonferroni99any/iter

   cat("\n MULTIPLE-TEST CORRECTION PROCEDURES: \n")
   cat("\n FDR-BH (Benjamini & Hochberg, 1995): \n")
   cat("   Total percent > 90% CL =", fdrBH90all,"   % spectra =",fdrBH90any,"\n")
   cat("   Total percent > 95% CL =", fdrBH95all,"   % spectra =",fdrBH95any,"\n")
   cat("   Total percent > 99% CL =", fdrBH99all,"   % spectra =",fdrBH99any,"\n\n")   
   cat(" FDR-BY (Benjamini & Yekutieli, 2001): \n")
   cat("   Total percent > 90% CL =", fdrBY90all,"   % spectra =",fdrBY90any,"\n")
   cat("   Total percent > 95% CL =", fdrBY95all,"   % spectra =",fdrBY95any,"\n")
   cat("   Total percent > 99% CL =", fdrBY99all,"   % spectra =",fdrBY99any,"\n\n")  
   cat(" Hommel (1988): \n")
   cat("   Total percent > 90% CL =", hommel90all,"   % spectra =",hommel90any,"\n")
   cat("   Total percent > 95% CL =", hommel95all,"   % spectra =",hommel95any,"\n")
   cat("   Total percent > 99% CL =", hommel99all,"   % spectra =",hommel99any,"\n\n")
   cat(" Hochberg (1988): \n")
   cat("   Total percent > 90% CL =", hochberg90all,"   % spectra =",hochberg90any,"\n")
   cat("   Total percent > 95% CL =", hochberg95all,"   % spectra =",hochberg95any,"\n")
   cat("   Total percent > 99% CL =", hochberg99all,"   % spectra =",hochberg99any,"\n\n")  
   cat(" Holm (1979): \n")
   cat("   Total percent > 90% CL =", holm90all,"   % spectra =",holm90any,"\n")
   cat("   Total percent > 95% CL =", holm95all,"   % spectra =",holm95any,"\n")
   cat("   Total percent > 99% CL =", holm99all,"   % spectra =",holm99any,"\n\n")
   cat(" Bonferroni: \n")
   cat("   Total percent > 90% CL =", bonferroni90all,"   % spectra =",bonferroni90any,"\n")
   cat("   Total percent > 95% CL =", bonferroni95all,"   % spectra =",bonferroni95any,"\n")
   cat("   Total percent > 99% CL =", bonferroni99all,"   % spectra =",bonferroni99any,"\n\n")
 }


### plot results
den90=density(frac90)
den95=density(frac95)
den99=density(frac99)
par(mfrow=c(2,1))
plot(den90,lwd=1.5,type="l", col="red",main="Percentage of spectrum exceeding 99%, 95% and 90% CL", xlab="Percent",xlim=c(0,100), ylim=c(0,max(den90$y,den95$y,den99$y)))
lines(den95,lwd=1.5,col="purple")
lines(den99,lwd=1.5,col="green")
### draw vertical line at true value
abline(v=10,lwd=1.5,col="red",lty=3)
abline(v=5,lwd=1.5,col="purple",lty=3)
abline(v=1,lwd=1.5,col="green",lty=3)
mtext(paste("Median Estimates=",round(med99,digits=2),"%,",round(med95,digits=2),"%,", round(med90,digits=2),"%      Total=",round(tick99all,digits=2),"%,",round(tick95all,digits=2),"%,", round(tick90all,digits=2),"%"),side=3,line=0,at=50)

p90=100*sum90[1:ifreq]/iter
p95=100*sum95[1:ifreq]/iter
p99=100*sum99[1:ifreq]/iter
plot(freq,p90,type="l", col="red",main="Distribution of false positives across spectrum",ylab="% False Positives", xlab="Frequency",xlim=c(0,Nyq),ylim=c(0,max(p90,p95,p99)))
lines(freq,p95, col="purple")
lines(freq,p99, col="green")
mtext(paste("Total number of simulations=",iter),side=3,line=0,at=mean(freq))
abline(h=10,col="red",lty=3)
abline(h=5,col="purple",lty=3)
abline(h=1,col="green",lty=3)

if(output)
 {
  frameout1 <- cbind(frac90,frac95,frac99)
  print(frameout1)
  frameout2 <- data.frame(cbind(freq,sum90,sum95,sum99))
  return(frameout2)
  }

# end function testBackground
}