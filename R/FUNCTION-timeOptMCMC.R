##################################################################################
### function timeOptMCMC - (SRM: November 1-29, 2016; December 5, 2016; 
###                              April 1, 2017; April 10, 2017; 
###                              April 14, 2017; April 16, 2017; April 17, 2017;
###                              April 20-21, 2017; May 11, 2018)
##################################################################################

# prior name: FUNCTION-timeOptMCMCgaus_v27.R

timeOptMCMC <- function (dat,iopt=1,sedmin=0.5,sedmax=5,sedstart=NULL,gAve=NULL,gSd=NULL,gstart=NULL,kAve=NULL,kSd=NULL,kstart=NULL,rhomin=0,rhomax=0.9999,rhostart=NULL,sigmamin=NULL,sigmamax=NULL,sigmastart=NULL,ran=F,fit=1,ftol=0.01,roll=10^3,nsamples=1000,epsilon=NULL,test=F,burnin=-1,detrend=T,output=1,savefile=F,genplot=1,verbose=T)
{

#######################################################################################
#### (1) data preparation
#######################################################################################
if(verbose) 
 {
   cat("\n----- TimeOptMCMC: Assessment of Amplitude Modulation & Bundling -----\n")
   cat("--    Optimization with uncertainties using Metropolis-Hastings     --\n")
   cat("-----            Markov Chain-Monte Carlo Algorithm              -----\n\n")
 }

if(test) cat("\n**** WARNING: epsilon testing mode activated!\n\n\n")

cormethod=1

if(fit !=1) stop("**** this option not yet functional!")

# prepare data array
   dat = data.frame(dat)      
   npts <- length(dat[,1]) 
   dx <- dat[2,1]-dat[1,1]

# initial error checking 
   if(dx<0)
     { 
       if (verbose) cat("\n * Sorting data into increasing height/depth/time, removing empty entries\n")
       dat <- dat[order(dat[1], na.last = NA, decreasing = F), ]
       dx <- dat[2,1]-dat[1,1]
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
   cat(" * Stratigraphic series length (meters):",(npts-1)*dx,"\n")
   cat(" * Sampling interval (meters):",dx,"\n\n")
 }

# detrend
if (detrend) 
  {
    lm.1 <- lm(dat[,2] ~ dat[,1])
    dat[2] <- dat[2] - (lm.1$coeff[2]*dat[1] + lm.1$coeff[1])
    if(verbose) cat(" * Linear trend subtracted. m=",lm.1$coeff[2],"b=",lm.1$coeff[1],"\n")
  }

# standardize data series
   dat[2]=dat[2]-colMeans(dat[2])
   dat[2]=dat[2]/sapply(dat[2],sd)


#######################################################################################
#### (2) definition of key FUNCTIONS: calcPeriods, randnrw, rwrand, genCycles, 
####      fitIt, calcLogLH
#######################################################################################
# ---------- FUNCTION calcPeriods ----------
# function to calculate periods in ka, given g and k    
calcPeriods <- function(g,k)
 {
# eccentricity periods in ka
# in order of relative amplitude
# 405.6 ka: g2-g5
  e1 = 1296/(g[2]-g[5])
# 94.9 ka: g4-g5
  e2 = 1296/(g[4]-g[5])
# 123.9 ka: g4-g2
  e3= 1296/(g[4]-g[2])
# 98.9 ka: g3-g5
  e4= 1296/(g[3]-g[5])
# 130.7 ka: g3-g2
  e5= 1296/(g[3]-g[2])
# precession periods in ka
# in order of relative amplitude
# 23.7 ka: k+g5
  p1= 1296/(k+g[5])
# 22.4 ka: k+g2
  p2= 1296/(k+g[2])
# 18.9 ka: k+g4
  p3= 1296/(k+g[4])
# 19.1 ka: k+g[3]
  p4= 1296/(k+g[3]) 
# 23.1 ka: k+g1
  p5= 1296/(k+g[1])
# output from longest to shortest period  
  return(c(e1,e5,e3,e4,e2,p1,p5,p2,p4,p3)) 
 } 

# ---------- FUNCTION randnrw ----------
# randnrw is an R translation of Alberto Malinverno's randnrw MATLAB function.
#  it returns a new value sampled in a random walk from the current value
#  xcurr.  The random walk samples are distributed as in a normal pdf.  
#  must call separately for each innovation. 
#  on first call, set xcurr to (xmax+xmin)/2
# epsilon controls how big the jump is between each iteration, and 
#  is the width of the uniform proposal distribution as a fraction 
#  of twice the standard deviation.
randnrw <- function(xcurr,xmean,xsdev,epsilon)
 {
# value of x if xcand is not accepted
   x=xcurr
   width=epsilon*2*xsdev
   xcand=xcurr+(width*runif(1))-width/2
   deltaxcurr=xcurr-xmean
   deltaxcand=xcand-xmean
   exparg=(deltaxcand^2-deltaxcurr^2)/(2*xsdev^2)
   probacc=exp(-exparg)
   if (probacc>runif(1)) x=xcand
   return(x)
  }

# ---------- FUNCTION rwrand ----------
# rwrand is a translation of Alberto Malinverno's MATLAB function.
#  it will perform a markov chain random walk that samples a 
#  uniform distribution. you must call it separately for each innovation. 
#  on first call, set xcurr to (xmax+xmin)/2
# epsilon controls how big the jump is between each iteration.
#  a value of 0.2 will yield jumps up to +/- 10% of the total
#  range (xmax-xmin).

rwrand <- function(xcurr,xmin,xmax,epsilon)
{
# value of x if xcand is not accepted
  x=xcurr
# sample only if the range is > 0  
  if (xmax>xmin)
   {
 	 halfwidth=0.5*epsilon*(xmax-xmin)
	 xcandhigh=xcurr+halfwidth
	 xcandlow=xcurr-halfwidth
	 xcand=xcandlow+runif(1)*(xcandhigh-xcandlow)
     if (xcand<=xmax && xcand>=xmin) x=xcand
     return(x)
   }
# if xmax=xmin, we simply return input value. this added for case when
#  we want to hold certain parameters constant.   
  if(xmin==xmax) return(x)
}  

# ---------- FUNCTION genCycles ----------
# function to generate cos (real) and sin (imaginary) terms for each target period, 
#   and convert to spatial cycles, given a particular sed rate in m/ka
genCycles <- function(sedrate1, targetIn, n) 
  {
    x <- matrix(0, n, 2*length(targetIn))
    for (i in 1:length(targetIn)) 
      {
        x[,2*i-1] <- cos( (2*pi)/(targetIn[i]) * (dx/sedrate1) * (1:n))
        x[,2*i] <- sin( (2*pi)/(targetIn[i]) * (dx/sedrate1) * (1:n))
      }
    return(x)
  }

# ---------- FUNCTION fitIt ----------
# function to perform fitting and calculate likelihood
#  dx, npts, flow, fhigh, roll passed into function transparently
#  g is a vector of five fundamental frequencies in the following order: g1,g2,g3,g4,g5
#  k contains the precession frequency
fitIt <- function(dat,sedrate1,rhoSpec,sigmaE,rhoEnv,sigmaEenv,g,k) 
  {
# CALIBRATE DEPTH SERIES (m) TO TIME (ka)
    timeSeries = dat
# create new time vector
# it is the index vector for time
    it <- seq(1,npts,by=1)
    time = (dx/sedrate1) * (it-1)
    timeSeries[1] = time
# use g and k to determine eccentricity and precession periods
#  these are output in order from longest period to shortest
    targetIn=calcPeriods(g=g,k=k)
    targetE=targetIn[1:5]
    targetP=targetIn[6:10]
    iflag=F
    if(min(targetIn)<0) iflag=T
#### spectral power fit
    xm <- genCycles(sedrate1, targetIn, npts)
    lm.0 <- lm(timeSeries[,2] ~ xm)  
    logLHspec = calcLogLH(lm.0$residuals,rhoSpec,sigmaE)
# here we calculate rho based on lm.0$residuals, simply for comparison
    rho_spec_m_NOW= cor(lm.0$residuals[1:(npts-1)],lm.0$residuals[2:npts])
    if(iopt == 1)
     {
#### envelope fit
# bandpass precession or short eccentricity band
      bp = taner(timeSeries,padfac=2,flow=(1/targetP[1])-ftol,fhigh=(1/targetP[5])+ftol,roll=roll,demean=T,detrend=F,addmean=F,genplot=F,verbose=F)
# hilbert transform for instantaneous amplitude
      hil = hilbert(bp,padfac=2,demean=T,detrend=F,addmean=F,genplot=F,verbose=F)
# standardize hil to unit variance
      hil[2]=hil[2]-colMeans(hil[2])
      hil[2]=hil[2]/sd(hil[,2])    
# the first 10 columns of xm contain the predicted eccentricity cos and sine terms
      xm1 <- xm[,1:10]
      lm.1 <- lm(hil[,2] ~ xm1)
      logLHenv=calcLogLH(lm.1$residuals,rhoEnv,sigmaEenv)
# here we calculate rho based on lm.1$residuals, simply for comparison
      rho_env_m_NOW= cor(lm.1$residuals[1:(npts-1)],lm.1$residuals[2:npts]) 
      logLH=logLHspec+logLHenv
      return(cbind(logLH,sd(lm.0$residuals),rho_spec_m_NOW,sd(lm.1$residuals),rho_env_m_NOW,iflag))
    }
    if(iopt==2) 
     {
       logLH=logLHspec
       return(cbind(logLH,sd(lm.0$residuals),rho_spec_m_NOW,iflag))
     }
   } 

# ---------- FUNCTION calcLogLH ----------
# function to perform log-likelihood calculation, including
#  assessment of correlated residuals, assuming AR1 model
#  npts passed into function transparently
calcLogLH <- function(residuals,rho,sigma)
 {
# sigma coming into function is log(sigma)
  sigma2=exp(sigma)
# calculate Re^-1 (ReInv), as in EQ. A-7 of Malinverno & Briggs (2004)
# set up array, with zeros
   ReInv <- double(npts*npts)
   dim(ReInv) <-c(npts,npts)
# put 1+rho^2 on diagonal
   ReInv[row(ReInv)==col(ReInv)] = 1+rho^2
# except at 1,1 and 100,100, which have a value of 1
   ReInv[1,1] = 1
   ReInv[npts,npts] = 1
# put -rho on subdiagonal  
   ReInv[(row(ReInv)-1)==col(ReInv)] = -1*rho
# put -rho on superdiagonal  
   ReInv[(row(ReInv)+1)==col(ReInv)] = -1*rho
# now multiple matrix by 1/(1-rho^2)
   ReInv = ReInv/(1-rho^2)
# calculate log-likelihood 
   logLH2 = ( npts*log(2*pi) ) + ( 2*npts*log(sigma2) ) + ( (npts-1)*log(1-rho^2) )
   logLH2 = logLH2 + (1/(sigma2^2)) * t(residuals) %*% ReInv %*% residuals 
   logLH2 = -0.5*logLH2
   return(logLH2)
 }


#######################################################################################
#### (3) ensure resolution requirements are satisfied
#######################################################################################
# convert sedmin and sedmax from cm/ka to m/ka for processing
   sedmin=sedmin/100
   sedmax=sedmax/100

# check minimum and maximum sedimentation rates. sedmin is now in m/ka, dx is in meters.
   NyqFreq=sedmin/(2*dx)
# freqHigh identifies the frequency of the shortest period in the target, given the 
#  stated uncertainties in g and k. we add ftol for the maximum possible upper half power 
#  point in the taner filtering.  we will take 2 standard deviations of mean as limit.
   freqHigh=1/min(calcPeriods((gAve+2*gSd),(kAve+2*kSd)))+ftol
   if(freqHigh>NyqFreq)
    {
     sedmin = 2*dx*freqHigh
     if(verbose) cat("\n**** WARNING: minimum sedimentation rate is too low.\n")
     if(verbose) cat("              sedmin reset to",100*sedmin,"cm/ka\n")
    }

# check maximum sedimentation rate. sedmax is in m/ka. dx is in meters.
    RayFreq = sedmax/(npts*dx)
# freqLow identifies the frequency of the longest period in the target, given the 
#  stated uncertainties in g. we will not worry about ftol here. we will take 2 standard 
#  deviations of mean as limit.
   freqLow=1/max(calcPeriods((gAve-2*gSd),(kAve-2*kSd)))
   if(RayFreq>freqLow)
     {
       sedmax = npts*dx*freqLow
       if(verbose) cat("\n**** WARNING: maximum sedimentation rate is too high.\n")
       if(verbose) cat("              sedmax reset to",100*sedmax,"cm/ka\n")
     }

   if(sedmin>sedmax)
     {
       cat("\n**** ERROR: sedmin > sedmax\n")
       stop("**** TERMINATING NOW!")
     }


#######################################################################################
#### (4) implement Metropolis-Hastings algorithm for Markov-Chain Monte Carlo. 
#### This is translated from Alberto Malinverno's MATLAB code.
#######################################################################################

### SET UP INITIAL PARAMETERS
# set default epsilon values to 0.2
if(is.null(epsilon)) epsilon=rep(0.2,11)

# intial values to kick-off simulations. note that sedrate is in m/ka,
#  and frequencies are in arcsec/yr
# recall that sedstart is in cm/ka, and we want m/ka
if(!is.null(sedstart) && sedstart>0) mstart=sedstart/100
# sedmin and sedmax have already been coverted from cm/ka to m/ka
if(!is.null(sedstart) && sedstart<=0) mstart=runif(1,min=sedmin,max=sedmax)
# set default sedstart
if(is.null(sedstart)) mstart=(sedmax+sedmin)/2

m=mstart
mcand=mstart

if(!is.null(kstart) && kstart <=0) kstart=rnorm(1,mean=kAve,sd=kSd)
if(is.null(kstart)) kstart=kAve
k=kstart
kcand=kstart

if(!is.null(gstart))
 {
   if(gstart[1] <=0) g1start=rnorm(1,mean=gAve[1],sd=gSd[1])
   if(gstart[1] >0) g1start=gstart[1]
   if(gstart[2] <=0) g2start=rnorm(1,mean=gAve[2],sd=gSd[2])
   if(gstart[2] >0) g2start=gstart[2]
   if(gstart[3] <=0) g3start=rnorm(1,mean=gAve[3],sd=gSd[3])
   if(gstart[3] >0) g3start=gstart[3]
   if(gstart[4] <=0) g4start=rnorm(1,mean=gAve[4],sd=gSd[4])
   if(gstart[4] >0) g4start=gstart[4]
   if(gstart[5] <=0) g5start=rnorm(1,mean=gAve[5],sd=gSd[5])
   if(gstart[5] >0) g5start=gstart[5]
 }  
 
if(is.null(gstart))
 {
   g1start=gAve[1]
   g2start=gAve[2]
   g3start=gAve[3]
   g4start=gAve[4]
   g5start=gAve[5]
 }  

g1=g1start
g1cand=g1start

g2=g2start
g2cand=g2start

g3=g3start
g3cand=g3start

g4=g4start
g4cand=g4start

g5=g5start
g5cand=g5start

# this will be used to initialize both spectral and envelope searches
if(!is.null(rhostart)) 
 {
   if(rhostart <=0) mstartRho=runif(1,min=rhomin,max=rhomax)
   if(rhostart >0) mstartRho=rhostart
 }

if(is.null(rhostart)) mstartRho=0.5
mRho=mstartRho
mcandRho=mstartRho
mRhoEnv=mstartRho
mcandRhoEnv=mstartRho

# this will be used to initialize both spectral and envelope searches
# we will use the raw data variance to estimate a plausible maximum for residual sigma
#  (if the record is entirely noise, this would be the value). 
# note that dat is no longer standardized to unit variance
# we will sample from log(sigma), to ensure uniform distribution
sdDat=sd(dat[,2])
if(is.null(sigmamax)) sigmaEmax=log(sdDat)
# here we assume that the signal/noise ratio is always < 10000
if(is.null(sigmamin)) sigmaEmin=log(sdDat*1/10000)
if(!is.null(sigmamax)) sigmaEmax=log(sigmamax)
if(!is.null(sigmamin)) sigmaEmin=log(sigmamin)
 
if(!is.null(sigmastart)) 
 {
   if(sigmastart <=0) sigmaEstart=runif(1,min=sigmaEmin,max=sigmaEmax)
   if(sigmastart >0) sigmaEstart=log(sigmastart)
 }
# default guess: signal/noise ratio of unity (half of sigma is noise)
if(is.null(sigmastart)) sigmaEstart=log(sdDat/2)
msigmaE=sigmaEstart
mcandsigmaE=sigmaEstart
msigmaEenv=sigmaEstart
mcandsigmaEenv=sigmaEstart

if(verbose) 
  {  
    cat("\n * Parameter bounds for prior distributions:\n")
    cat("     Sedimentation rate (cm/ka)=",sedmin*100,"to",sedmax*100,"\n")
    cat("     k=",kAve,"+/-",kSd,"\n")
    cat("     g1=",gAve[1],"+/-",gSd[1],"\n")
    cat("     g2=",gAve[2],"+/-",gSd[2],"\n")
    cat("     g3=",gAve[3],"+/-",gSd[3],"\n")
    cat("     g4=",gAve[4],"+/-",gSd[4],"\n")
    cat("     g5=",gAve[5],"+/-",gSd[5],"\n")
    
    cat("\n * Hyperparameter bounds for prior distributions:\n")
    cat("     Spectral Power Fit Residual rho=",rhomin,"to",rhomax,"\n")
    cat("     Spectral Power Fit Residual sigma=",exp(sigmaEmin),"to",exp(sigmaEmax),"\n")
    cat("     Envelope Fit Residual rho=",rhomin,"to",rhomax,"\n")
    cat("     Envelope Fit Residual sigma=",exp(sigmaEmin),"to",exp(sigmaEmax),"\n")
    
    cat("\n * Initial values for each parameter/hyperparameter:\n")
    cat("     Sedimentation rate (cm/ka)=",mstart*100,"\n")
    cat("     k=",kstart,"\n")
    cat("     g1=",g1start,"\n")
    cat("     g2=",g2start,"\n")
    cat("     g3=",g3start,"\n")
    cat("     g4=",g4start,"\n")
    cat("     g5=",g5start,"\n")
    cat("     Spectral Power Fit:\n") 
    cat("        Residual rho=",mstartRho,"\n")
    cat("        Residual sigma=",exp(sigmaEstart),"\n")
    if(iopt==1) 
     {
       cat("     Envelope Fit:\n") 
       cat("        Residual rho=",mstartRho,"\n")
       cat("        Residual sigma=",exp(sigmaEstart),"\n")    
     }
  }

# this is the initial value for logpdf
if(!test) 
 {
   fitItRes=fitIt(dat,mstart,mstartRho,sigmaEstart,mstartRho,sigmaEstart,c(g1start,g2start,g3start,g4start,g5start),kstart) 
   logpdf=fitItRes[1]
   sigma_spec_m = fitItRes[2]
   rho_spec_m = fitItRes[3]
   if(iopt==1)
     {
       sigma_env_m = fitItRes[4]
       rho_env_m = fitItRes[5]
     }  
# we should set these to something here, so the code doesn't crash. let's choose -1 as a 'dummy' place holder.    
   if(iopt==2)
     {
       sigma_env_m = -1
       rho_env_m = -1
     }  
# stop now if your initial values of g and k are physcially impossible
   if(iopt==1 && fitItRes[6]) 
      {
        cat("\n**** ERROR: some of the periods are negative!\n")
        stop("**** TERMINATING NOW!")
      }
    if(iopt==2 && fitItRes[4])
      {
        cat("\n**** ERROR: some of the periods are negative!\n")
        stop("**** TERMINATING NOW!")
      }  
 }
  
# for epsilon testing, set logliklihood to 1   
if(test) 
 {
   logpdf=1
# we aren't interested in the terms below, but will set them so the code doesn't crash. let's choose -1 as a 'dummy' place holder.  
   sigma_env_m = -1
   rho_env_m = -1
   sigma_env_m = -1
   rho_env_m = -1   
 }  

# prepare for loop
if(verbose) 
  {  
    cat("\n * PLEASE WAIT: Performing Optimization\n")
    cat("\n0%       25%       50%       75%       100%\n")
# create a progress bar
    progress = utils::txtProgressBar(min = 0, max = nsamples, style = 1, width=43)
  }

# if genplot == 2, also display a progress plot. this is useful for testing parameters,
#  particularly epsilon values and their impact on stabilty and efficiency
if(genplot==2) 
{
#   dev.new(title="MCMC Progress")
# for faster plot updates, use x11. note that this is not suitable for all platforms (Windows).
#   x11()
   dev.new()
   par(mfrow = c(2, 1),mar = c(5.1, 4.1, 4.1, 2.1))
}   


# initialize nacc, which will keep track of the number of accepted candidates
nacc=0

# now for each parameter (this only useful when ran=T)
nacc_Sed=0
nacc_k=0
nacc_g1=0
nacc_g2=0
nacc_g3=0
nacc_g4=0
nacc_g5=0
nacc_rho=0
nacc_sigma=0
nacc_rhoEnv=0
nacc_sigmaEnv=0

# keep track of each time the parameters was selected (this is only useful when ran=T)
n_Sed=0
n_k=0
n_g1=0
n_g2=0
n_g3=0
n_g4=0
n_g5=0
n_rho=0
n_sigma=0
n_rhoEnv=0
n_sigmaEnv=0

# initialize matrix to store simulations
samplemat=double(nsamples*17)
dim(samplemat) = c(nsamples,17)
# the first column will contain loglikelihood
# the second columns will give sedimentation rate
# third column is k
# fourth-eigth columns are g terms
# ninth column is rho for spectral fit
# tenth column is sigmaE for spectral fit
# eleventh column is rho for envelope fit
# twelfth column is sigmaE for envelope fit
# thirteenth column contains number for each successful iteration (can monitor for efficiency)
# fourteenth column is measured rho for spectral fit
# fifteenth column is measured sigmaE for spectral fit
# sixteenth column is measured rho for envelope fit
# seventeenth column is measured sigmaE for envelope fit

# assign column titles
colnames(samplemat) <- c("log-likelihood","sedrate","k","g1","g2","g3","g4","g5","rho_spec","sigma_spec","rho_env","sigma_env","iteration","rho_spec_m","sigma_spec_m","rho_env_m","sigma_env_m")


#############################################################
####  Metropolis-Hastings sampling loop #####################
for(iter in 1:nsamples)
 {
# report progress
    if(verbose) utils::setTxtProgressBar(progress, iter)

#   the vector epsilon has the following order: sedrate,k,g1,g2,g3,g4,g5,spec_rho,spec_sigma,env_rho,env_sigma
# if desired, randomly select variable to update, ranging from 1 to 11
#   if(ran && iopt == 1) select = as.integer(floor(runif(1,min=1,max=12)))
#   if(ran && iopt == 2) select = as.integer(floor(runif(1,min=1,max=10)))


#### APRIL 20, 2017: sample parameters with varying frequency, based on sensitiviy tests
   if(ran) pick = runif(1,min=0,max=1)
# 25% prob of sedrate   
   if(pick >= 0 && pick <0.25) select = 1
# 25% prob of k  
   if(pick >= 0.25 && pick <0.5) select = 2
# 30% prob of g's  
   if(pick >= 0.5 && pick <0.8) select = 3
# 10% prob of sigma_spec and rho_spec  
   if(pick >= 0.8 && pick <0.9) select = 4       
# 10% prob of sigma_env and rho_env  
   if(pick >= 0.9 && pick <1) select = 5 


# to hardwire search for sedrate and k only, use next line
#   if(ran && iopt == 2) select = as.integer(floor(runif(1,min=1,max=3)))  
   if(!ran) select = 0

# generate a new candidate sedimentation rate
   if(select == 0 || select == 1) 
    {
      mcand=rwrand(m,xmin=sedmin,xmax=sedmax,epsilon=epsilon[1])
      n_Sed=n_Sed+1
    }  
# generate a new candidate k
   if(select == 0 || select == 2) 
    {
      kcand=randnrw(k,xmean=kAve,xsdev=kSd,epsilon=epsilon[2])
      n_k=n_k+1
    } 
# generate new g terms
   if(select == 0 || select == 3) 
    {
      g1cand=randnrw(g1,xmean=gAve[1],xsdev=gSd[1],epsilon=epsilon[3])
      n_g1=n_g1+1
    }
   if(select == 0 || select == 3) 
    {
      g2cand=randnrw(g2,xmean=gAve[2],xsdev=gSd[2],epsilon=epsilon[4])
      n_g2=n_g2+1
    }  
   if(select == 0 || select == 3) 
    {
      g3cand=randnrw(g3,xmean=gAve[3],xsdev=gSd[3],epsilon=epsilon[5])
      n_g3=n_g3+1
    }
   if(select == 0 || select == 3) 
    {
      g4cand=randnrw(g4,xmean=gAve[4],xsdev=gSd[4],epsilon=epsilon[6])
      n_g4=n_g4+1
    }  
   if(select == 0 || select == 3) 
    {
      g5cand=randnrw(g5,xmean=gAve[5],xsdev=gSd[5],epsilon=epsilon[7])
      n_g5=n_g5+1
    }  
# generate a new candidate rho
   if(select == 0 || select == 4) 
    {
      mcandRho=rwrand(mRho,xmin=rhomin,xmax=rhomax,epsilon=epsilon[8])
      n_rho=n_rho+1
     } 
   if(select == 0 || select == 4) 
    {
      mcandsigmaE=rwrand(msigmaE,xmin=sigmaEmin,xmax=sigmaEmax,epsilon=epsilon[9])
      n_sigma=n_sigma+1
    }  
   if(iopt==1)
    {
      if(select == 0 || select == 5) 
       {
         mcandRhoEnv=rwrand(mRhoEnv,xmin=0,xmax=0.9999,epsilon=epsilon[10])
         n_rhoEnv=n_rhoEnv+1
        } 
      if(select == 0 || select == 5) 
       {
         mcandsigmaEenv=rwrand(msigmaEenv,xmin=sigmaEmin,xmax=sigmaEmax,epsilon=epsilon[11])
         n_sigmaEnv=n_sigmaEnv+1
       }  
    }
   
# evaluate candidates with fitIt
   if(!test) 
    {
      fitItRes= fitIt(dat,mcand,mcandRho,mcandsigmaE,mcandRhoEnv,mcandsigmaEenv,c(g1cand,g2cand,g3cand,g4cand,g5cand),kcand)
      logpdfcand= fitItRes[1] 

# ensure that the astronomical periods are positive. if not, skip this iteration, and use prior candidate
      if(iopt==1 && fitItRes[6] || iopt==2 && fitItRes[4]) 
       {
         cat("\n**** WARNING: Some orbital periods are negative. Skipping this candiate! (iter =)",iter,"\n")
# on error, return to previous values         
         logpdfcand=logpdf
         mcand=m
         mcandRho=mRho
         mcandsigmaE=msigmaE
         mcandRhoEnv=mRhoEnv
         mcandsigmaEenv=msigmaEenv
         kcand=k
         g1cand=g1
         g2cand=g2
         g3cand=g3
         g4cand=g4
         g5cand=g5
         logpdfcand=logpdf
         nacc=nacc-1    
         fitItRes[2]= sigma_spec_m
         fitItRes[3]=rho_spec_m
         if(iopt == 1)
           {
             fitItRes[4]=sigma_env_m
             fitItRes[5]=rho_env_m 
           }           
       }
    }  

# for epsilon testing, set logliklihood to 1 
   if(test) logpdfcand=1

# determine probability of acceptance. this is equivalent to the ratio of the two likelihoods
    probacc=exp(logpdfcand-logpdf)
#    cat(iter,"sedrate=",mcand,"; old=",logpdf,"; new=",logpdfcand,"; prob. accept.=",probacc,"\n")
     
# generate a random number between 0 and 1, compare your observed prob of acceptance to it.
#  note that runif does not generate the extreme values (0 and 1).   
    if (probacc>runif(1))
       {
          m=mcand
          mRho=mcandRho
          msigmaE=mcandsigmaE
          mRhoEnv=mcandRhoEnv
          msigmaEenv=mcandsigmaEenv
          k=kcand
          g1=g1cand
          g2=g2cand
          g3=g3cand
          g4=g4cand
          g5=g5cand
          logpdf=logpdfcand
          nacc=nacc+1
# fitIt returns: logLH,lm.0$residuals,rho_spec_m_NOW,lm.1$residuals,rho_env_m_NOW
          if(!test)
           {
             sigma_spec_m = fitItRes[2]
             rho_spec_m = fitItRes[3]
             if(iopt == 1)
               {
                 sigma_env_m = fitItRes[4]
                 rho_env_m = fitItRes[5]
               }
             if(iopt == 2)
               {
# assign some placeholder values here
                 sigma_env_m = -1
                 rho_env_m = -1
               }   
            }     
          if(test)
           {
# assign some placeholder values here
             sigma_spec_m = -1
             rho_spec_m = -1
             sigma_env_m = -1
             rho_env_m = -1
            }     
            
# update acceptance rates for each parameter
#  NOTE: the acceptance rates are not correct if any of the periods above were <0.
#  could modify this portion of the code to allow for such instances.            

          if(select == 1) nacc_Sed=nacc_Sed+1
          if(select == 2) nacc_k=nacc_k+1
 
          if(select == 3) nacc_g1=nacc_g1+1
          if(select == 3) nacc_g2=nacc_g2+1
          if(select == 3) nacc_g3=nacc_g3+1
          if(select == 3) nacc_g4=nacc_g4+1
          if(select == 3) nacc_g5=nacc_g5+1

          if(select == 4) nacc_rho=nacc_rho+1
          if(select == 4) nacc_sigma=nacc_sigma+1
          
          if(select == 5) nacc_rhoEnv=nacc_rhoEnv+1
          if(select == 5) nacc_sigmaEnv=nacc_sigmaEnv+1   
          
#          if(select == 4) nacc_g2=nacc_g2+1
#          if(select == 5) nacc_g3=nacc_g3+1
#          if(select == 6) nacc_g4=nacc_g4+1
#          if(select == 7) nacc_g5=nacc_g5+1
#          if(select == 8) nacc_rho=nacc_rho+1
#          if(select == 9) nacc_sigma=nacc_sigma+1
#          if(select == 10) nacc_rhoEnv=nacc_rhoEnv+1
#          if(select == 11) nacc_sigmaEnv=nacc_sigmaEnv+1        

# plot is placed in this loop to avoid generating too many graphs       
# updated after every other accepted candidate  
          if (genplot == 2 && nacc > 1 && nacc %% 2 == 0)
            {
# note that this is always one candidate behind, but that's fine for our progress assessment            
              plot(samplemat[1:(iter-1),1],xlab="MCMC sample number",ylab="Log-likelihood",main="Stability")
              plot(1:(iter-1),100*samplemat[1:(iter-1),13]/(1:(iter-1)),xlab="MCMC sample number",ylab="% of Candidates Accepted",main="Efficiency")
            }         
        } 
# save results to matrix samplemat
    samplemat[iter,1]=logpdf
# convert sedimentation rate back to cm/ka
    samplemat[iter,2]=m*100
    samplemat[iter,3]=k
    samplemat[iter,4]=g1
    samplemat[iter,5]=g2   
    samplemat[iter,6]=g3    
    samplemat[iter,7]=g4    
    samplemat[iter,8]=g5 
    samplemat[iter,9]=mRho  
# convert log(sigma) to sigma
    samplemat[iter,10]=exp(msigmaE)   
    samplemat[iter,11]=mRhoEnv
# convert log(sigma) to sigma      
    samplemat[iter,12]=exp(msigmaEenv)  
    samplemat[iter,13]=nacc  
# these are the actual noise (residual) parameters estimated from the fit (not the candidate values)
    samplemat[iter,14]=rho_spec_m
    samplemat[iter,15]=sigma_spec_m 
    samplemat[iter,16]=rho_env_m
    samplemat[iter,17]=sigma_env_m

# if savefile=T, output results after every 1000 candidates have been evaluated
if(savefile && iter == 1000) write.table(file ="MCMCsamples.csv", samplemat[1:1000,], sep = ",", row.names = FALSE, append=FALSE)
if(savefile && iter > 1000 && iter %% 1000 == 0) write.table(file ="MCMCsamples.csv", samplemat[(iter-999):iter,], sep = ",", row.names = FALSE, append=TRUE, col.names=FALSE)
# note if you choose iter that is not an even thousand, the last iterations will not be output

# end Metropolis sampling loop
}
if(verbose) close(progress)


#######################################################################################
#### (5) post-process MCMC results, generate summaries
#######################################################################################

# Find MAP (maximium a posteriori point) and detect nburnin
# extract column 1 of samplemat
# imap is the index
imap=which.max(samplemat[,1])
logpdfmap=samplemat[imap,1] 
 
logpdfstart=samplemat[1,1]

if(!test)
 {
# burn in detection, using fraction threshold criterion
  if(burnin>=0)
   {
     burnintol=burnin  
     logpdfburnin=logpdfstart+burnintol*(logpdfmap-logpdfstart)
# this loop could be vectorized
     for(ii in 1:nsamples)
      {
        if (samplemat[ii,1]>logpdfburnin) 
         {
           nburnin=ii 
           break
         } 
      }   
   }  

# burn in detection, using median value from second half of simulated likelihoods
 if(burnin<0) 
  {
    logpdfburnin=median(samplemat[(nsamples/2):nsamples,1])
# identify the first value to above it, use that as the start of the burn-in  
    nburnin=which(samplemat[,1]>logpdfburnin)[1]
  }
 }
 
if(test) 
 {
# really this should be set to 0, but will likely cause some plots below to fail! 
   nburnin = 1
   logpdfburnin = 1
   cat("\n EPSILON TESTING ACTIVATED: no burnin applied \n")
 }  


cat("\n", nacc, "candidates accepted out of", nsamples,"samples =",100*nacc/nsamples,"%\n")
if(ran)
 {
    cat("\n   INDIVIDUAL PARAMETER ACCEPTANCE RATES:\n")
    cat("\n     sed rate ( out of",n_Sed,")=",100*nacc_Sed/n_Sed,"%\n")
    cat("\n     k ( out of",n_k,")=",100*nacc_k/n_k,"%\n")
    cat("\n     g1 ( out of",n_g1,")=",100*nacc_g1/n_g1,"%\n")
    cat("\n     g2 ( out of",n_g2,")=",100*nacc_g2/n_g2,"%\n")
    cat("\n     g3 ( out of",n_g3,")=",100*nacc_g3/n_g3,"%\n")    
    cat("\n     g4 ( out of",n_g4,")=",100*nacc_g4/n_g4,"%\n")    
    cat("\n     g5 ( out of",n_g5,")=",100*nacc_g5/n_g5,"%\n")
    cat("\n     rho_spec ( out of",n_rho,")=",100*nacc_rho/n_rho,"%\n")
    cat("\n     sigma_spec ( out of",n_sigma,")=",100*nacc_sigma/n_sigma,"%\n")
    if(iopt==1)
     {
       cat("\n     rho_env ( out of",n_rhoEnv,")=",100*nacc_rhoEnv/n_rhoEnv,"%\n")
       cat("\n     sigma_env ( out of",n_sigmaEnv,")=",100*nacc_sigmaEnv/n_sigmaEnv,"%\n")
     }
 }

cat(" Metropolis: initial value of log-likelihood =",logpdfstart,"\n")
cat(" Metropolis: assumed value of log-likelihood at burnin =",logpdfburnin,"\n")
cat(" Number of viable MCMC candidates post-burnin=",nsamples-nburnin,"\n\n")

cat(" Metropolis: maximum log-likelihood is",logpdfmap, "at following parameter values:\n")
cat("  Sedimentation rate (m/ka)=",samplemat[imap,2],"\n")
cat("  Precession constant k (arcsec/yr)=",samplemat[imap,3],"\n")
cat("  g1 frequency (arcsec/yr)=",samplemat[imap,4],"\n")
cat("  g2 frequency (arcsec/yr)=",samplemat[imap,5],"\n")
cat("  g3 frequency (arcsec/yr)=",samplemat[imap,6],"\n")
cat("  g4 frequency (arcsec/yr)=",samplemat[imap,7],"\n")
cat("  g5 frequency (arcsec/yr)=",samplemat[imap,8],"\n")
cat("  Spectral power residual rho=",samplemat[imap,9],"\n")
cat("  Spectral power residual sigma=",samplemat[imap,10],"\n")
if(iopt==1) cat("  Envelope residual rho=",samplemat[imap,11],"\n")
if(iopt==1) cat("  Envelope residual sigma=",samplemat[imap,12],"\n\n")

# calculate optimal periods
MAPperiods=calcPeriods(g=c(samplemat[imap,4],samplemat[imap,5],samplemat[imap,6],samplemat[imap,7],samplemat[imap,8]),k=samplemat[imap,3])
cat("  Reconstructed Eccentricity Periods:\n")
cat("  ",MAPperiods[1],",",MAPperiods[2],",",MAPperiods[3],",",MAPperiods[4],",",MAPperiods[5])
cat("\n  Reconstructed Precession Periods:\n")
cat("  ",MAPperiods[6],",",MAPperiods[7],",",MAPperiods[8],",",MAPperiods[9],",",MAPperiods[10])


#######################################################################################
##### summary plots
#######################################################################################

if(genplot == 1 || genplot ==2 )cexset= 0.85
if(genplot == 3 )
  {
    pngwidth= 800
    pngheight= 800
    pngres= 130
    cexset= 0.75
  }
   
if(genplot > 0)
 {   
   if(genplot == 1 || genplot == 2) dev.new(title="MCMC Stability and Efficiency")
   if(genplot == 3) png(filename = "stability.png", width = pngwidth, height = pngheight, units = "px", res=pngres)

   par(mfrow = c(2, 1),mar = c(5.1, 4.1, 4.1, 2.1))
# zoom-in for plotting
#   nburninplot=min(nsamples,nburnin*10)
# let's see all results
   nburninplot=nsamples
# log-likelihood versus iteration
   plot(samplemat[1:nburninplot,1],xlab="Sample number",ylab="Log-likelihood",main="Stability",cex=cexset)
   points(samplemat[1:nburnin,1],col="red",cex=cexset)
# monitor efficiency
#   plot(1:nburninplot,samplemat[1:nburninplot,13],xlab="Sample number",ylab="Iteration number",main="Efficiency")
#   points(1:nburnin,samplemat[1:nburnin,13],col="red")
   plot(1:nburninplot,100*samplemat[1:nburninplot,13]/(1:nsamples),xlab="Sample number",ylab="% of Candidates Accepted",main="Efficiency",cex=cexset)
   points(1:nburnin,100*samplemat[1:nburnin,13]/(1:nburnin),col="red",cex=cexset)

   if(genplot==3) dev.off()

# sampled values (including burn-in)
  if(genplot == 1 || genplot == 2) dev.new(title="Evolution of Candidate Values throughout MCMC Simulation")
  if(genplot == 3) png(filename = "evolution_with_burnin.png", width = pngwidth, height = pngheight*1.1, units = "px",res=pngres)
  par(mfrow = c(3, 3),mar = c(5.1, 4.1, 4.1, 2.1))
  plot(samplemat[,2],samplemat[,1],xlab="cm/ka",ylab="log-likelihood",main="Sedimentation rate",col="blue",cex=cexset)
  points(samplemat[1:nburnin,2],samplemat[1:nburnin,1],col="red",cex=cexset)
  points(samplemat[imap,2],samplemat[imap,1],cex=2,col="white",pch=21)
  points(samplemat[imap,2],samplemat[imap,1],cex=2,col="mediumblue",pch=16)
  abline(v=c(sedmin,sedmax)*100,lty=4,col="slategray",lwd=1.5)

  plot(samplemat[,3],samplemat[,1],xlab="arcsec/yr",ylab="log-likelihood",main="k frequency",col="blue",cex=cexset)
  points(samplemat[1:nburnin,3],samplemat[1:nburnin,1],col="red",cex=cexset)
  points(samplemat[imap,3],samplemat[imap,1],cex=2,col="white",pch=21)
  points(samplemat[imap,3],samplemat[imap,1],cex=2,col="mediumblue",pch=16)
  abline(v=c(kAve-kSd,kAve+kSd),lty=4,col="slategray",lwd=1.5)

  plot(samplemat[,4],samplemat[,1],xlab="arcsec/yr",ylab="log-likelihood",main="g1 frequency",col="blue",cex=cexset)
  points(samplemat[1:nburnin,4],samplemat[1:nburnin,1],col="red",cex=cexset)
  points(samplemat[imap,4],samplemat[imap,1],cex=2,col="white",pch=21)
  points(samplemat[imap,4],samplemat[imap,1],cex=2,col="mediumblue",pch=16)
  abline(v=c(gAve[1]-gSd[1],gAve[1]+gSd[1]),lty=4,col="slategray",lwd=1.5)

  plot(samplemat[,5],samplemat[,1],xlab="arcsec/yr",ylab="log-likelihood",main="g2 frequency",col="blue",cex=cexset)
  points(samplemat[1:nburnin,5],samplemat[1:nburnin,1],col="red",cex=cexset)
  points(samplemat[imap,5],samplemat[imap,1],cex=2,col="white",pch=21)
  points(samplemat[imap,5],samplemat[imap,1],cex=2,col="mediumblue",pch=16)
  abline(v=c(gAve[2]-gSd[2],gAve[2]+gSd[2]),lty=4,col="slategray",lwd=1.5)

  plot(samplemat[,6],samplemat[,1],xlab="arcsec/yr",ylab="log-likelihood",main="g3 frequency",col="blue",cex=cexset)
  points(samplemat[1:nburnin,6],samplemat[1:nburnin,1],col="red",cex=cexset)
  points(samplemat[imap,6],samplemat[imap,1],cex=2,col="white",pch=21)
  points(samplemat[imap,6],samplemat[imap,1],cex=2,col="mediumblue",pch=16)
  abline(v=c(gAve[3]-gSd[3],gAve[3]+gSd[3]),lty=4,col="slategray",lwd=1.5)

  plot(samplemat[,7],samplemat[,1],xlab="arcsec/yr",ylab="log-likelihood",main="g4 frequency",col="blue",cex=cexset)
  points(samplemat[1:nburnin,7],samplemat[1:nburnin,1],col="red",cex=cexset)
  points(samplemat[imap,7],samplemat[imap,1],cex=2,col="white",pch=21)
  points(samplemat[imap,7],samplemat[imap,1],cex=2,col="mediumblue",pch=16)
  abline(v=c(gAve[4]-gSd[4],gAve[4]+gSd[4]),lty=4,col="slategray",lwd=1.5)

  plot(samplemat[,8],samplemat[,1],xlab="arcsec/yr",ylab="log-likelihood",main="g5 frequency",col="blue",cex=cexset)
  points(samplemat[1:nburnin,8],samplemat[1:nburnin,1],col="red",cex=cexset)
  points(samplemat[imap,8],samplemat[imap,1],cex=2,col="white",pch=21)
  points(samplemat[imap,8],samplemat[imap,1],cex=2,col="mediumblue",pch=16)
  abline(v=c(gAve[5]-gSd[5],gAve[5]+gSd[5]),lty=4,col="slategray",lwd=1.5)
  
  xlim2=c(min(samplemat[,9],samplemat[,11]),max(samplemat[,9],samplemat[,11]))
  plot(samplemat[,9],samplemat[,1],xlab=expression(rho),ylab="log-likelihood",main=expression(paste("Residual ", rho)),xlim=xlim2,col="blue",cex=cexset)
  mtext(c("(blue=power; purple=envelope)"),side=3,cex=0.75)
  points(samplemat[1:nburnin,9],samplemat[1:nburnin,1],col="red",cex=cexset)
  points(samplemat[imap,9],samplemat[imap,1],cex=2,col="white",pch=21)
  points(samplemat[imap,9],samplemat[imap,1],cex=2,col="mediumblue",pch=16)
  if(iopt==1)
   {
     points(samplemat[,11],samplemat[,1],col="purple",cex=cexset)
     points(samplemat[1:nburnin,11],samplemat[1:nburnin,1],col="red",cex=cexset)
     points(samplemat[imap,11],samplemat[imap,1],cex=2,col="white",pch=21)
     points(samplemat[imap,11],samplemat[imap,1],cex=2,col="purple4",pch=16)
   }
  abline(v=c(rhomin,rhomax),lty=4,col="slategray",lwd=1.5)

  xlim2=c(min(samplemat[,10],samplemat[,12]),max(samplemat[,10],samplemat[,12]))  
  plot(samplemat[,10],samplemat[,1],xlab=expression(sigma),ylab="log-likelihood",main=expression(paste("Residual ", sigma)),xlim=xlim2,col="blue",cex=cexset)
  mtext(c("(blue=power; purple=envelope)"),side=3,cex=0.75)
  points(samplemat[1:nburnin,10],samplemat[1:nburnin,1],col="red",cex=cexset)
  points(samplemat[imap,10],samplemat[imap,1],cex=2,col="white",pch=21)
  points(samplemat[imap,10],samplemat[imap,1],cex=2,col="mediumblue",pch=16)
  if(iopt==1)
   {
     points(samplemat[,12],samplemat[,1],col="purple",cex=cexset)
     points(samplemat[1:nburnin,12],samplemat[1:nburnin,1],col="red",cex=cexset)
     points(samplemat[imap,12],samplemat[imap,1],cex=2,col="white",pch=21)
     points(samplemat[imap,12],samplemat[imap,1],cex=2,col="purple4",pch=16)
   }
  abline(v=c(exp(sigmaEmin),exp(sigmaEmax)),lty=4,col="slategray",lwd=1.5)

  if(genplot == 3) dev.off()
  
# sampled values (post burn-in)
  if(genplot == 1 || genplot == 2) dev.new(title="Evolution of Candidate Values Post Burn-in")
  if(genplot == 3) png(filename = "evolution_no_burnin.png", width = pngwidth, height = pngheight*1.1, units = "px", res=pngres)
  par(mfrow = c(3, 3),mar = c(5.1, 4.1, 4.1, 2.1))
  plot(samplemat[(nburnin+1):nsamples,2],samplemat[(nburnin+1):nsamples,1],xlab="cm/ka",ylab="log-likelihood",main="Sedimentation rate",col="blue",cex=cexset)
  points(samplemat[imap,2],samplemat[imap,1],cex=2,col="white",pch=21)
  points(samplemat[imap,2],samplemat[imap,1],cex=2,col="mediumblue",pch=16)
  abline(v=c(sedmin,sedmax)*100,lty=4,col="slategray",lwd=1.5)

  plot(samplemat[(nburnin+1):nsamples,3],samplemat[(nburnin+1):nsamples,1],xlab="arcsec/yr",ylab="log-likelihood",main="k frequency",col="blue",cex=cexset)
  points(samplemat[imap,3],samplemat[imap,1],cex=2,col="white",pch=21)
  points(samplemat[imap,3],samplemat[imap,1],cex=2,col="mediumblue",pch=16)
  abline(v=c(kAve-kSd,kAve+kSd),lty=4,col="slategray",lwd=1.5)  

  plot(samplemat[(nburnin+1):nsamples,4],samplemat[(nburnin+1):nsamples,1],xlab="arcsec/yr",ylab="log-likelihood",main="g1 frequency",col="blue",cex=cexset)
  points(samplemat[imap,4],samplemat[imap,1],cex=2,col="white",pch=21)
  points(samplemat[imap,4],samplemat[imap,1],cex=2,col="mediumblue",pch=16)
  abline(v=c(gAve[1]-gSd[1],gAve[1]+gSd[1]),lty=4,col="slategray",lwd=1.5)

  plot(samplemat[(nburnin+1):nsamples,5],samplemat[(nburnin+1):nsamples,1],xlab="arcsec/yr",ylab="log-likelihood",main="g2 frequency",col="blue",cex=cexset)
  points(samplemat[imap,5],samplemat[imap,1],cex=2,col="white",pch=21)
  points(samplemat[imap,5],samplemat[imap,1],cex=2,col="mediumblue",pch=16)
  abline(v=c(gAve[2]-gSd[2],gAve[2]+gSd[2]),lty=4,col="slategray",lwd=1.5)

  plot(samplemat[(nburnin+1):nsamples,6],samplemat[(nburnin+1):nsamples,1],xlab="arcsec/yr",ylab="log-likelihood",main="g3 frequency",col="blue",cex=cexset)
  points(samplemat[imap,6],samplemat[imap,1],cex=2,col="white",pch=21)
  points(samplemat[imap,6],samplemat[imap,1],cex=2,col="mediumblue",pch=16)
  abline(v=c(gAve[3]-gSd[3],gAve[3]+gSd[3]),lty=4,col="slategray",lwd=1.5)

  plot(samplemat[(nburnin+1):nsamples,7],samplemat[(nburnin+1):nsamples,1],xlab="arcsec/yr",ylab="log-likelihood",main="g4 frequency",col="blue",cex=cexset)
  points(samplemat[imap,7],samplemat[imap,1],cex=2,col="white",pch=21)
  points(samplemat[imap,7],samplemat[imap,1],cex=2,col="mediumblue",pch=16)
  abline(v=c(gAve[4]-gSd[4],gAve[4]+gSd[4]),lty=4,col="slategray",lwd=1.5)

  plot(samplemat[(nburnin+1):nsamples,8],samplemat[(nburnin+1):nsamples,1],xlab="arcsec/yr",ylab="log-likelihood",main="g5 frequency",col="blue",cex=cexset)
  points(samplemat[imap,8],samplemat[imap,1],cex=2,col="white",pch=21)
  points(samplemat[imap,8],samplemat[imap,1],cex=2,col="mediumblue",pch=16)
  abline(v=c(gAve[5]-gSd[5],gAve[5]+gSd[5]),lty=4,col="slategray",lwd=1.5)
  
  xlim2=c(min(samplemat[(nburnin+1):nsamples,9],samplemat[(nburnin+1):nsamples,11]),max(samplemat[(nburnin+1):iter,9],samplemat[(nburnin+1):iter,11]))
  plot(samplemat[(nburnin+1):nsamples,9],samplemat[(nburnin+1):nsamples,1],xlab=expression(rho),ylab="log-likelihood",main=expression(paste("Residual ", rho)),xlim=xlim2,col="blue",cex=cexset)
  mtext(c("(blue=power; purple=envelope)"),side=3,cex=0.75)
  points(samplemat[imap,9],samplemat[imap,1],cex=2,col="white",pch=21)
  points(samplemat[imap,9],samplemat[imap,1],cex=2,col="mediumblue",pch=16)
  if(iopt==1)
   {
     points(samplemat[(nburnin+1):nsamples,11],samplemat[(nburnin+1):nsamples,1],col="purple",cex=cexset)
     points(samplemat[imap,11],samplemat[imap,1],cex=2,col="white",pch=21)
     points(samplemat[imap,11],samplemat[imap,1],cex=2,col="purple4",pch=16)
   }
  abline(v=c(rhomin,rhomax),lty=4,col="slategray",lwd=1.5)

  xlim2=c(min(samplemat[(nburnin+1):nsamples,10],samplemat[(nburnin+1):nsamples,12]),max(samplemat[(nburnin+1):nsamples,10],samplemat[(nburnin+1):nsamples,12]))  
  plot(samplemat[(nburnin+1):nsamples,10],samplemat[(nburnin+1):nsamples,1],xlab=expression(sigma),ylab="log-likelihood",main=expression(paste("Residual ", sigma)),xlim=xlim2,col="blue",cex=cexset)
  mtext(c("(blue=power; purple=envelope)"),side=3,cex=0.75)
  points(samplemat[imap,10],samplemat[imap,1],cex=2,col="white",pch=21)
  points(samplemat[imap,10],samplemat[imap,1],cex=2,col="mediumblue",pch=16)
  if(iopt==1)
   {
     points(samplemat[(nburnin+1):nsamples,12],samplemat[(nburnin+1):nsamples,1],col="purple",cex=cexset)
     points(samplemat[imap,12],samplemat[imap,1],cex=2,col="white",pch=21)
     points(samplemat[imap,12],samplemat[imap,1],cex=2,col="purple4",pch=16)
   }
  abline(v=c(exp(sigmaEmin),exp(sigmaEmax)),lty=4,col="slategray",lwd=1.5)

  if(genplot==3) dev.off()

  if(genplot == 1 || genplot == 2) dev.new(title="Kernel Densities for MCMC Parameter Estimates")
  if(genplot == 3) png(filename = "kernel.png", width = pngwidth, height = pngheight*1.05, units = "px", res=pngres)
  par(mfrow = c(3, 3),mar = c(5.1, 4.1, 4.1, 2.1))
  plot(density(samplemat[(nburnin+1):nsamples,2]),xlab="cm/ka",main="Sedimentation rate",lwd=2,col="blue",xlim=c(sedmin,sedmax)*100)
  polygon(density(samplemat[(nburnin+1):nsamples,2]),col="#0000FF64",border=NA)
  abline(v=samplemat[imap,2],lty=3,lwd=2,col="red")
  abline(v=c(sedmin,sedmax)*100,lty=4,col="slategray",lwd=1.5)
  plot(density(samplemat[(nburnin+1):nsamples,3]),xlab="arcsec/yr",main="k frequency",lwd=2,col="blue",xlim=c(kAve-3*kSd,kAve+3*kSd))
  polygon(density(samplemat[(nburnin+1):nsamples,3]),col="#0000FF64",border=NA)
  abline(v=samplemat[imap,3],lty=3,lwd=2,col="red")
  abline(v=c(kAve-kSd,kAve+kSd),lty=4,col="slategray",lwd=1.5)
  plot(density(samplemat[(nburnin+1):nsamples,4]),xlab="arcsec/yr",main="g1 frequency",lwd=2,col="blue",xlim=c(gAve[1]-3*gSd[1],gAve[1]+3*gSd[1]))
  polygon(density(samplemat[(nburnin+1):nsamples,4]),col="#0000FF64",border=NA)
  abline(v=samplemat[imap,4],lty=3,lwd=2,col="red")
  abline(v=c(gAve[1]-gSd[1],gAve[1]+gSd[1]),lty=4,col="slategray",lwd=1.5)
  plot(density(samplemat[(nburnin+1):nsamples,5]),xlab="arcsec/yr",main="g2 frequency",lwd=2,col="blue",xlim=c(gAve[2]-3*gSd[2],gAve[2]+3*gSd[3]))
  polygon(density(samplemat[(nburnin+1):nsamples,5]),col="#0000FF64",border=NA)
  abline(v=samplemat[imap,5],lty=3,lwd=2,col="red")
  abline(v=c(gAve[2]-gSd[2],gAve[2]+gSd[2]),lty=4,col="slategray",lwd=1.5)
  plot(density(samplemat[(nburnin+1):nsamples,6]),xlab="arcsec/yr",main="g3 frequency",lwd=2,col="blue",xlim=c(gAve[3]-3*gSd[3],gAve[3]+3*gSd[3]))
  polygon(density(samplemat[(nburnin+1):nsamples,6]),col="#0000FF64",border=NA)
  abline(v=samplemat[imap,6],lty=3,lwd=2,col="red")
  abline(v=c(gAve[3]-gSd[3],gAve[3]+gSd[3]),lty=4,col="slategray",lwd=1.5)
  plot(density(samplemat[(nburnin+1):nsamples,7]),xlab="arcsec/yr",main="g4 frequency",lwd=2,col="blue",xlim=c(gAve[4]-3*gSd[4],gAve[4]+3*gSd[4]))
  polygon(density(samplemat[(nburnin+1):nsamples,7]),col="#0000FF64",border=NA)
  abline(v=samplemat[imap,7],lty=3,lwd=2,col="red")
  abline(v=c(gAve[4]-gSd[4],gAve[4]+gSd[4]),lty=4,col="slategray",lwd=1.5)
  plot(density(samplemat[(nburnin+1):nsamples,8]),xlab="arcsec/yr",main="g5 frequency",lwd=2,col="blue",xlim=c(gAve[5]-3*gSd[5],gAve[5]+3*gSd[5]))
  polygon(density(samplemat[(nburnin+1):nsamples,8]),col="#0000FF64",border=NA)
  abline(v=samplemat[imap,8],lty=3,lwd=2,col="red")
  abline(v=c(gAve[5]-gSd[5],gAve[5]+gSd[5]),lty=4,col="slategray",lwd=1.5)
  den1=density(samplemat[(nburnin+1):nsamples,9])
  if(iopt==1) den2=density(samplemat[(nburnin+1):nsamples,11])
  if(iopt==2) den2=den1
#  xlim2=c(min(den1$x,den2$x),max(den1$x,den2$x))
  xlim2=c(rhomin,rhomax)
  ylim2=c(min(den1$y,den2$y),max(den1$y,den2$y))
  plot(density(samplemat[(nburnin+1):nsamples,9]),xlab=expression(rho),main=expression(paste("Residual ", rho)),xlim=xlim2,ylim=ylim2,lwd=2,col="blue")
  mtext(c("(blue=power; purple=envelope)"),side=3,cex=0.75)
  polygon(density(samplemat[(nburnin+1):nsamples,9]),col="#0000FF64",border=NA)
  if(iopt==1)
   {
    lines(density(samplemat[(nburnin+1):nsamples,11]),col="purple",lwd=2)
    polygon(density(samplemat[(nburnin+1):nsamples,11]),col="#A020F064",border=NA)
    abline(v=samplemat[imap,11],lty=3,lwd=2,col="red")
    }
  abline(v=samplemat[imap,9],lty=3,lwd=2,col="red")  
  abline(v=c(rhomin,rhomax),lty=4,col="slategray",lwd=1.5)
  den1=density(samplemat[(nburnin+1):nsamples,10])
  if(iopt==1) den2=density(samplemat[(nburnin+1):nsamples,12])
  if(iopt==2) den2=den1
#  xlim2=c(min(den1$x,den2$x),max(den1$x,den2$x))
  xlim2=c(exp(sigmaEmin),exp(sigmaEmax))
  ylim2=c(min(den1$y,den2$y),max(den1$y,den2$y))
  plot(density(samplemat[(nburnin+1):nsamples,10]),xlab=expression(sigma),main=expression(paste("Residual ", sigma)),xlim=xlim2,ylim=ylim2,lwd=2,col="blue")
  polygon(density(samplemat[(nburnin+1):nsamples,10]),col="#0000FF64",border=NA)
  mtext(c("(blue=power; purple=envelope)"),side=3,cex=0.75)
  if(iopt==1)
   {
    lines(density(samplemat[(nburnin+1):nsamples,12]),col="purple",lwd=2)
    polygon(density(samplemat[(nburnin+1):nsamples,12]),col="#A020F064",border=NA)
    abline(v=samplemat[imap,12],lty=3,lwd=2,col="red")
   }
  abline(v=samplemat[imap,10],lty=3,lwd=2,col="red") 
  abline(v=c(exp(sigmaEmin),exp(sigmaEmax)),lty=4,col="slategray",lwd=1.5)

  if(genplot==3) dev.off()

  if(genplot == 1 || genplot == 2) dev.new(title="Histograms for MCMC Parameter Estimates")
  if(genplot == 3) png(filename = "histogram.png", width = pngwidth, height = pngheight*1.2, units = "px", res=pngres)
  par(mfrow = c(4, 3),mar = c(5.1, 4.1, 4.1, 2.1))
  hist(samplemat[(nburnin+1):nsamples,2],xlab="cm/ka",main="Sedimentation rate",lwd=2,col="#0000FF64",breaks="FD",xlim=c(sedmin,sedmax)*100)
  abline(v=samplemat[imap,2],lty=3,lwd=2,col="red")
  abline(v=c(sedmin,sedmax)*100,lty=4,col="slategray",lwd=1.5)
  hist(samplemat[(nburnin+1):nsamples,3],xlab="arcsec/yr",main="k frequency",lwd=2,col="#0000FF64",breaks="FD",xlim=c(kAve-3*kSd,kAve+3*kSd))
  abline(v=samplemat[imap,3],lty=3,lwd=2,col="red")
  abline(v=c(kAve-kSd,kAve+kSd),lty=4,col="slategray",lwd=1.5)
  hist(samplemat[(nburnin+1):nsamples,4],xlab="arcsec/yr",main="g1 frequency",lwd=2,col="#0000FF64",breaks="FD",xlim=c(gAve[1]-3*gSd[1],gAve[1]+3*gSd[1]))
  abline(v=samplemat[imap,4],lty=3,lwd=2,col="red")
  abline(v=c(gAve[1]-gSd[1],gAve[1]+gSd[1]),lty=4,col="slategray",lwd=1.5)
  hist(samplemat[(nburnin+1):nsamples,5],xlab="arcsec/yr",main="g2 frequency",lwd=2,col="#0000FF64",breaks="FD",xlim=c(gAve[2]-3*gSd[2],gAve[2]+3*gSd[2]))
  abline(v=samplemat[imap,5],lty=3,lwd=2,col="red")
  abline(v=c(gAve[2]-gSd[2],gAve[2]+gSd[2]),lty=4,col="slategray",lwd=1.5)
  hist(samplemat[(nburnin+1):nsamples,6],xlab="arcsec/yr",main="g3 frequency",lwd=2,col="#0000FF64",breaks="FD",xlim=c(gAve[3]-3*gSd[3],gAve[3]+3*gSd[3]))
  abline(v=samplemat[imap,6],lty=3,lwd=2,col="red")
  abline(v=c(gAve[3]-gSd[3],gAve[3]+gSd[3]),lty=4,col="slategray",lwd=1.5)
  hist(samplemat[(nburnin+1):nsamples,7],xlab="arcsec/yr",main="g4 frequency",lwd=2,col="#0000FF64",breaks="FD",xlim=c(gAve[4]-3*gSd[4],gAve[4]+3*gSd[4]))
  abline(v=samplemat[imap,7],lty=3,lwd=2,col="red")
  abline(v=c(gAve[4]-gSd[4],gAve[4]+gSd[4]),lty=4,col="slategray",lwd=1.5)
  hist(samplemat[(nburnin+1):nsamples,8],xlab="arcsec/yr",main="g5 frequency",lwd=2,col="#0000FF64",breaks="FD",xlim=c(gAve[5]-3*gSd[5],gAve[5]+3*gSd[5]))
  abline(v=samplemat[imap,8],lty=3,lwd=2,col="red")
  abline(v=c(gAve[5]-gSd[5],gAve[5]+gSd[5]),lty=4,col="slategray",lwd=1.5)
  hist(samplemat[(nburnin+1):nsamples,9],xlab=expression(rho),main=expression(paste("Spectral Power Residual ", rho)),lwd=2,col="#0000FF64",breaks="FD",xlim=c(rhomin,rhomax))
  abline(v=samplemat[imap,9],lty=3,lwd=2,col="red")
  abline(v=c(rhomin,rhomax),lty=4,col="slategray",lwd=1.5)
  if(iopt==1)
   {
     hist(samplemat[(nburnin+1):nsamples,11],xlab=expression(rho),main=expression(paste("Envelope Residual ", rho)),lwd=2,col="#A020F064",breaks="FD",xlim=c(rhomin,rhomax))
     abline(v=samplemat[imap,11],lty=3,lwd=2,col="red")
     abline(v=c(rhomin,rhomax),lty=4,col="slategray",lwd=1.5)
   }
  hist(samplemat[(nburnin+1):nsamples,10],xlab=expression(sigma),main=expression(paste("Spectral Power Residual ", sigma)),lwd=2,col="#0000FF64",breaks="FD",xlim=c(exp(sigmaEmin),exp(sigmaEmax)))
  abline(v=samplemat[imap,10],lty=3,lwd=2,col="red")
  abline(v=c(exp(sigmaEmin),exp(sigmaEmax)),lty=4,col="slategray",lwd=1.5)
  if(iopt==1)
   {
     hist(samplemat[(nburnin+1):nsamples,12],xlab=expression(sigma),main=expression(paste("Envelope Residual ", sigma)),lwd=2,col="#A020F064",breaks="FD",xlim=c(exp(sigmaEmin),exp(sigmaEmax)))
     abline(v=samplemat[imap,12],lty=3,lwd=2,col="red")
     abline(v=c(exp(sigmaEmin),exp(sigmaEmax)),lty=4,col="slategray",lwd=1.5)
   }
  if(genplot==3) dev.off()
  if(genplot == 3) cat("\n * Summary plots saved in directory", getwd(), "\n")

# end plot section
}


if(output == 1) return(samplemat)
 
### END function timeOptMCMC
}