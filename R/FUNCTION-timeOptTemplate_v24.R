### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2021 Stephen R. Meyers
###
###########################################################################
### function timeOptTemplate - (SRM: May 28, 2012; Oct. 14, 2014; Oct. 17, 2014; 
###                          Oct. 21, 2014; Jan. 8, 2015; March 9, 2015; Sept. 20, 2015
###                          Sept. 22, 2015; February 16-18, 2016; March 28, 2016; 
###                          May 4, 2016; May 10, 2016; May 17, 2016; 
###                          December 14-18, 2017; December 21, 2017; 
###                          November 22, 2018; November 24, 2018; December 2, 2018
###                          January 7, 2019; January 14, 2021)
###########################################################################

timeOptTemplate <- function (dat,template=NULL,sedmin=0.5,sedmax=5,difmin=NULL,difmax=NULL,fac=NULL,numsed=50,linLog=1,limit=T,fit=1,fitModPwr=T,iopt=3,flow=NULL,fhigh=NULL,roll=NULL,targetE=NULL,targetP=NULL,cormethod=1,detrend=T,detrendTemplate=F,flipTemplate=F,ncores=1,output=0,genplot=1,check=T,verbose=1)
{

#######################################################################################
# (1) prepare arrays, perform checks, initialize parameters
#######################################################################################
  if(verbose>0) 
   {
     cat("\n----- TimeOptTemplate: Assessment of Amplitude Modulation & Bundling----- \n\n")
     cat("        Evaluate variable sedimentation template, including: \n")
     cat("          - differential accumulation rate across bedding couplets \n")
     cat("          - linear sedimentation rate increase/decrease through interval \n")
     cat("          - step changes in sedimentation rate \n")
     cat("          - hiatus of unknown duration \n\n")
   }  

# Prepare data array
  dat = data.frame(dat)      
  npts = length(dat[,1]) 
  dx = dat[2,1]-dat[1,1]
  space = (npts-1)*dx

  if(check)
   {
# initial error checking 
     if(dx<0)
       { 
         if (verbose>0) cat("\n * Sorting data into increasing height/depth/time, removing empty entries\n")
         dat <- dat[order(dat[,1], na.last = NA, decreasing = F), ]
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
   }

  if(!is.null(template)) template = data.frame(template) 
  if(is.null(template)) 
   {
     if(verbose>0) cat(" * Using stratigraphic series (dat) as default template for differential accumulation optimization")
     template <- dat
   }  
# make sure each depth of template matches a depth in dat   
  if(sum(template[,1]-dat[,1]) != 0) stop("***** ERROR: each location in dat must have a matching location in template")

  if(flipTemplate)
    {
      template[2]=template[2]*-1
      if(verbose>0) cat("\n * Sedimentation template flipped.\n")
    }

  if (verbose>0) 
   {
     cat("\n * Number of data points in stratigraphic series:",npts,"\n")
     cat(" * Stratigraphic series length (space or time):",space,"\n")
     cat(" * Sampling interval (space or time):",dx,"\n\n")
   }

# detrend data series
  if(detrend) 
    {
      lm.1 <- lm(dat[,2] ~ dat[,1])
      dat[2] <- dat[2] - (lm.1$coeff[2]*dat[1] + lm.1$coeff[1])
      if(verbose>0) cat(" * Linear trend subtracted from data series. m=",lm.1$coeff[2],"b=",lm.1$coeff[1],"\n")
    }

# detrend template
  if(detrendTemplate) 
    {
      lm.2 <- lm(template[,2] ~ template[,1])
      template[2] <- template[2] - (lm.2$coeff[2]*template[1] + lm.2$coeff[1])
      if(verbose>0) cat(" * Linear trend subtracted from sedimentation template. m=",lm.2$coeff[2],"b=",lm.2$coeff[1],"\n")
    }

# standardize data series
  dat[2]=dat[2]-colMeans(dat[2])
  dat[2]=dat[2]/sapply(dat[2],sd)

# set error when (difmin,difmax) and fac specified
  if(all((!is.null(difmin)||!is.null(difmax)),!is.null(fac)))
   {
     cat("\n**** ERROR: You must specify either difmin & difmax, or fac\n")
     stop("**** TERMINATING NOW!")
   }  

  bound=1
  if(!is.null(difmax)) plotMax=difmax 
  if(!is.null(fac)) 
   {
     plotMax=sedmax*fac
     bound=2
   }  

# set default difmax
  if(is.null(difmax) && is.null(fac))
   {
     difmax=sedmax+(sedmax*0.1)
     plotMax=difmax
     if(verbose) cat("\n * difmax set to",difmax,"cm/kyr\n")
   }

# set default difmin
  if(is.null(difmin) && is.null(fac))
   {
     difmin=sedmin-(sedmin*0.1)
     if(verbose) cat(" * difmin set to",difmin,"cm/kyr\n")
   }
        
# convert difmin, difmax, sedmin and sedmax from cm/ka to m/ka for processing
  if(!is.null(difmin)) difmin=difmin/100
  if(!is.null(difmax)) difmax=difmax/100 
  sedmin=sedmin/100
  sedmax=sedmax/100

# set up default bandpass frequencies and targets
#  first for precession modulations
  if(fit == 1)
   {
     if(is.null(flow)) 
      {
        flow = 0.035
        if(verbose>0) cat(" * Using default flow =",flow,"\n")
      }  
     if(is.null(fhigh)) 
      {
        fhigh = 0.065
        if(verbose>0) cat(" * Using default fhigh =",fhigh,"\n")
      } 
     if(is.null(roll))
      {
        roll = 10^3
        if(verbose>0) cat(" * Using default roll =",roll,"\n")
      }   

     if(!is.null(targetP)) targetP <- as.double(targetP)
     if(is.null(targetP))
      {
# the four dominant precession peaks, based on spectral analysis of 
#   Laskar et al. (2004), 0-10 Ma
        targetP <- double(4)
        targetP[1] = 23.62069
        targetP[2] = 22.31868
        targetP[3] = 19.06768 
        targetP[4] = 18.91979 
        if(verbose>0) cat(" * Using default precession target periods (ka)=",targetP,"\n")
      }
   }

# next for short eccentricity modulations
  if(fit == 2)
   {
     if(is.null(flow))
      {
        flow = 0.007
        if(verbose>0) cat(" * Using default flow =",flow,"\n")
      }  
     if(is.null(fhigh))
      {
        fhigh = 0.0115
        if(verbose>0) cat(" * Using default fhigh =",fhigh,"\n")
      }
     if(is.null(roll))
      {
        roll = 10^5
        if(verbose>0) cat(" * Using default roll =",roll,"\n")
      }
     if(!is.null(targetP)) if(verbose) cat("\n**** WARNING: targetP is defined but will not be used in fitting!\n")
   }  
   
  if(!is.null(targetE)) targetE <- as.double(targetE)

  if(is.null(targetE))
   {
# the five domintant eccentricity peaks based on spectral analysis of LA10d solution 
#   (Laskar et al., 2011), 0-20 Ma
     targetE <- double(5)
     targetE[1] = 405.6795
     targetE[2] = 130.719
     targetE[3] = 123.839
     targetE[4] = 98.86307
     targetE[5] = 94.87666
     if(verbose) cat(" * Using default eccentricity target periods (ka)=",targetE,"\n")     
   }   

# targetTot is for plotting, and fitting if precession modulations assessed
  if(fit == 1 && fitModPwr) targetTot = c(targetE,targetP)
  if(fit == 1 && !fitModPwr) targetTot = c(targetP)
  if(fit == 2 && fitModPwr) targetTot = c(targetE)
  if(fit == 2 && !fitModPwr) targetTot = c(targetE[-1])

# remove duplicate values if present (this added to avoid user error)
  targetTot=unique(targetTot)

  if(flow < 0.5/max(targetTot) && verbose>0) cat("\n**** NOTE: flow is less than half of the smallest target frequency. Did you intend this?\n")
  if(fhigh > 2/min(targetTot) && verbose>0) cat("\n**** NOTE: fhigh is 2 times greater than the largest target frequency. Did you intend this?\n")  

# check minimum and maximum sedimentation rates. sedmin is now in m/ka, dx is in meters.
  NyqFreq=sedmin/(2*dx)
# fhigh is the upper half-power point for the filter   
  if(fhigh>NyqFreq)
   {
     if(limit) 
      {
        sedmin = 2*dx*fhigh
        if(verbose) cat("\n**** WARNING: minimum sedimentation rate is too low for full signal recovery.\n")
        if(verbose) cat("              sedmin reset to",100*sedmin,"cm/ka\n\n")
      }
     if(verbose>0 && !limit) 
      {
# note: when sedimentation rates exceed this value, the filter starts to behave more-and-more like a high-pass filter
       cat("\n**** WARNING: minimum sedimentation rate is too low for full signal recovery, but it will still be evaluated.\n")
       cat("              The high frequency half-power point of filter is exceeded when sedrate=",100*2*dx*fhigh,"cm/ka\n")
       cat("              USE WITH CAUTION!\n\n")
       if(fit == 1)
        {
          cat("              The shortest precession period is undetectable when sedrate <",100*2*dx/min(targetP),"cm/ka\n")
          cat("              The longest precession period is undetectable when sedrate <",100*2*dx/max(targetP),"cm/ka\n\n")
        }
        if(fit == 2)
         {   
           cat("              The shortest eccentricity period is undetectable when sedrate <",100*2*dx/min(targetE),"cm/ka\n")
           cat("              The longest eccentricity period is undetectable when sedrate <",100*2*dx/max(targetE),"cm/ka\n\n")
         }  
      }
   }

# check maximum sedimentation rate. sedmax is in m/ka. dx is in meters.
  RayFreq = sedmax/(npts*dx)
# freqLow identifies the frequency of the longest period in the target.
  freqLow=1/max(targetE)
  if(RayFreq>freqLow)
   {
    if(limit) 
     {
       sedmax = npts*dx*freqLow
       if(verbose>0) cat("\n**** WARNING: maximum sedimentation rate is too high for full signal recovery.\n")
       if(verbose>0) cat("              sedmax reset to",100*sedmax,"cm/ka\n\n")
     }
    if(!limit && verbose>0) 
     {
       cat("\n**** WARNING: maximum sedimentation rate is too high for full signal recovery, but it will still be evaluated.\n")
       cat("              The longest eccentricity period is undetectable when sedrate >",100*npts*dx*freqLow,"cm/ka\n")
       cat("              The shortest eccentricity period is undetectable when sedrate >",100*npts*dx/min(targetE),"cm/ka\n")
       cat("              USE WITH CAUTION!\n\n")
     }  
   }

  if(sedmin>sedmax)
   {
     cat("\n**** ERROR: sedmin > sedmax\n")
     stop("**** TERMINATING NOW!")
   }

# make a copy of template for testModel and makeModel
  scaled <- template

# when running in parallel, do not generate progress plots
  if(ncores>1 && genplot==2) genplot=1
# view progress
  if(genplot==2)
   {
      dev.new(height=6.5,width=9)
      par(mfrow=c(2,2))
   }    


#######################################################################################
# (2)  Definition of FUNCTIONS: genCycles, fitIt, testModel, makeModel
#######################################################################################
# genCycles called from fitIt. it is a function to generate cos (real) and sin (imaginary) 
#  terms for each target period, and convert them to spatial cycles, given a particular 
#  sedimentation rate template (in m/ka). note that timeModel is derived from function 
#  makeModel; it contains the sample times, given the modeled sedimentation history
genCycles <- function(timeModel, targetIn, n) 
  {
# generate cycle model    
    x <- matrix(0, n, 2*length(targetIn))
    for (i in 1:length(targetIn)) 
      {
        x[,2*i-1] <- cos( (2*pi)/(targetIn[i]) * timeModel[,1])
        x[,2*i] <- sin( (2*pi)/(targetIn[i]) * timeModel[,1])
      }
    return(x)
  }

# makeModel called from fitIt. note that srMax is always mapped to the maximum model value
makeModel <- function(srMin,srMax)  
  {
    slope = (srMax-srMin)/(max(template[,2])-min(template[,2]))
    b = srMax-(slope*max(template[,2]))
    scaled[,2] = (slope*template[,2])+b
# sedrate2time is expecting sedrates in cm/ka
    scaled[,2] = scaled[,2] * 100     
    spaceTimeMap = sedrate2time(scaled,check=F,verbose=F,genplot=F)
    return(spaceTimeMap)
   }
   
# testModel called from fitIt. we will optimize on the function testModel in fitIt, using 
#  Brent's method (note that "L-BFGS-B" doesn't converge)
#  dat, srMax, totTime and npts passed into function transparently
#  note that srMax is always mapped to the maximum model value
testModel <- function(srMin)  
  {
    slope = (srMax-srMin)/(max(template[,2])-min(template[,2]))
    b = srMax-(slope*max(template[,2]))
    scaled[,2] = (slope*template[,2])+b
# sedrate2time is expecting sedrates in cm/ka
    scaled[,2] = scaled[,2] * 100     
    spaceTimeMap = sedrate2time(scaled,check=F,verbose=F,genplot=F)
    T2=spaceTimeMap[npts,2]
# experimented with different powers here, to see how it influences tolerance for optimization       
    fit = (T2-totTime)^10
    return(fit)
  }
   
# function to perform fitting and calculate r-squared.
#  fit, targetTot, targetE, srAve, totTime, dat and dx, roll, difmin, space passed into function transparently
fitIt <- function(srMax) 
  {
# make srMax a global variable, so testModel can see it    
    srMax <<- srMax
# find optimal srMin given srMax and totTime. These sedrates are in m/ka; set maximum at srAve
    res1=optim(par=c(srAve),fn=testModel,method=c("Brent"),lower=difmin,upper=c(srAve),control=list(abstol=.Machine$double.eps))
    srMin=res1$par
    spaceTimeMap=makeModel(srMin,srMax)
    timeSeries=tune(dat,spaceTimeMap,check=F,genplot=F,verbose=F)   
    timeSeries2=linterp(timeSeries,check=F,genplot=F,verbose=F)
# missing one point after tuning sometimes?     
    npts2=length(timeSeries2[,2])
# GENERATE TIME-CALIBRATED MODEL for spectral power assessment (applies for fit=1 or fit=2) 
# there is an advantage to constructing the model at the original times, then interpolating.
#  but we will use this simplification
    xm.0 <- genCycles(timeSeries2, targetTot, npts2)
# measure fit 
    lm.0 <- lm(timeSeries2[,2] ~ xm.0)
# bandpass precession or short eccentricity band
    bp = taner(timeSeries2,padfac=2,flow=flow,fhigh=fhigh,roll=roll,demean=T,detrend=F,addmean=F,check=F,genplot=F,verbose=F)
# hilbert transform for instantaneous amplitude
    hil = hilbert(bp,padfac=2,demean=T,detrend=F,check=F,addmean=F,genplot=F,verbose=F)
# GENERATE TIME-CALIBRATED MODEL for envelope assessment
    if(fit==1) xm.1 <- genCycles(timeSeries2, targetE, npts2)
    if(fit==2) xm.1 <- genCycles(timeSeries2, targetE[1], npts2)
# measure fit
    lm.1 <- lm(hil[,2] ~ xm.1)
    if(cormethod==1) 
     {
       rsq_pwr = cor(timeSeries2[,2],lm.0$fitted,method=c("pearson"))^2
       rsq_mod = cor(hil[,2],lm.1$fitted,method=c("pearson"))^2
     }  
    if(cormethod==2) 
     {
       rsq_pwr = cor(timeSeries2[,2],lm.0$fitted,method=c("spearman"))^2
       rsq_mod = cor(hil[,2],lm.1$fitted,method=c("spearman"))^2
     }  
    rsq_tot = rsq_pwr * rsq_mod

##### PART BELOW IS TO VISUALIZE PROGRESS OF OPTIMIZATION, IF DESIRED #####      
    if(genplot==2)
     {       
       plot(timeSeries2,xlab="Time (ka)",ylab="Value",main="Calibrated Record", col="red",type="l")
       slope2 = (srMax-srMin)/(max(template[,2])-min(template[,2]))
       b2 = srMax-(slope2*max(template[,2]))
       scaled[,2] = (slope2*template[,2])+b2
# sedrate2time is expecting sedrates in cm/ka
       scaled[,2] = scaled[,2] * 100    
       plot(scaled,xlab="Depth (m)",ylab="Sedimentation rate (cm/ka)",main="Sedimentation Rates", col="red",type="l",ylim=c(0,plotMax))
       mtext(adj=0,round(srAve*100,digits=4),col="red")
# plot hil and ecc estimated by least squares fitting
       plot(hil,cex=.5,xlab="Time (ka)",ylab="Value",main="Hilbert AM vs. Reconstructed Eccentricity Fit", col="red",type="l")
       mtext(adj=1,round(rsq_mod,digits=2),col="red")
       par(new = TRUE)
       plot(hil[,1],lm.1$fitted,type="l",xaxt="n",yaxt="n",xlab="",ylab="")     
       fft = periodogram( timeSeries2, output=1, verbose=F,genplot=F)
# remove f(0) for log plot of power
       fft = subset(fft,(fft[,1] > 0))
       plot(fft[,1],fft[,3],xlim=c(0,0.1),type="l",xlab="Frequency (cycles/ka)",ylab="Power",main="Power Spectrum (black=linear; gray=log)")
       mtext(adj=0,round(rsq_pwr,digits=2),col="red")
# plot a second y-axis
       par(new=TRUE)
       plot(fft[,1],log(fft[,3]),xlim=c(0,0.1),type="l",yaxt="n",col="gray",xlab="",ylab="")
       abline(v=1/targetTot, col="red",lty=3)
       abline(v=c(flow,fhigh), col="blue",lty=2)
     }
##### PART ABOVE IS TO VISUALIZE PROGRESS OF OPTIMIZATION, IF DESIRED #####  

# at present, verbose=2 is the only means to report on convergence of testModel.
#  useful for testing purposes.
    if(iopt == 1) 
     {
       if(verbose>1)
        {
          if(res1$converge == 0) cat("SEDRATE (Min/Ave/Max):",srMin*100,100*space/totTime,srMax*100," r2=", rsq_mod," (converged)\n")
          if(res1$converge != 0) cat("SEDRATE (Min/Ave/Max):",srMin*100,100*space/totTime,srMax*100," r2=", rsq_mod," (DID NOT converge)\n")
        }  
       return(rsq_mod)
     }  
    if(iopt == 2) 
     {
       if(verbose>1)
        {
          if(res1$converge == 0) cat("SEDRATE (Min/Ave/Max):",srMin*100,"/",100*space/totTime,"/",srMax*100," r2=", rsq_pwr," (converged)\n")
          if(res1$converge != 0) cat("SEDRATE (Min/Ave/Max):",srMin*100,"/",100*space/totTime,"/",srMax*100," r2=", rsq_pwr," (DID NOT converge)\n")
        }
       return(rsq_pwr)
     }
    if(iopt == 3) 
     {
       if(verbose>1)
        {
          if(res1$converge == 0) cat("SEDRATE (Min/Ave/Max):",srMin*100,100*space/totTime,srMax*100," r2=", rsq_tot," (converged)\n")
          if(res1$converge != 0) cat("SEDRATE (Min/Ave/Max):",srMin*100,100*space/totTime,srMax*100," r2=", rsq_tot," (DID NOT converge)\n") 
        }
       return(rsq_tot)
     }     
# end function fitIt     
  }


#######################################################################################
# (3)  Set up sedimentation rate grid to evaluate
#######################################################################################
  sedrate <-double(numsed)
  if(linLog==0)
   {
     sedinc = (sedmax-sedmin)/(numsed-1)  
     sedrate = sedmin + 0:(numsed-1)*sedinc
    }
   
  if(linLog==1)
   {
     sedinc = (log10(sedmax)-log10(sedmin))/(numsed-1)
     sedrate = log10(sedmin) + 0:(numsed-1)*sedinc
     sedrate = 10^sedrate
   }


#######################################################################################
# (4)  Conduct grid search on one processor. This option is required if you want
#      to view the optimization in progress.
#######################################################################################
if(ncores==1)
{

# set up average sedimentation rate grid array, dimension appropriately
#  ansMessage<-character(numsed)
  ansSrMax<- double(numsed)
  ansConverge<- double(numsed)
  rPwr<-double(numsed)

  if(verbose==1) 
   {  
     cat("\n * PLEASE WAIT: Performing Optimization\n")
     cat("\n0%       25%       50%       75%       100%\n")
# create a progress bar
     progress = utils::txtProgressBar(min = 0, max = numsed, style = 1, width=43)
   }

# begin sedimentation rate loop
  i=0
  for(ii in 1:numsed) 
   {
     if(verbose==1) utils::setTxtProgressBar(progress, ii)
     if(verbose>1) cat("sedrate",ii,"=",sedrate[ii]*100,"\n")
     i=i+1
     srAve=sedrate[ii]
# convert average sedimentation rate to total time
     totTime = (max(dat[,1])-min(dat[,1]))/srAve
# execute functions
# note that fitIt is searching and optmizing on srMax, given srAve
#  optimal srMin is determined in testModel, given each tested srMax and totTime
# let's set maximum sedimentation rate to fac*srAve. default of 5 is based
#   on experimentation. If larger than this, risk getting into local minimum during fit.
# May 17, 2016: modified to allow for difmax; modified again Dec 21, 2017
     if(bound==1) setmax=difmax
     if(bound==2) setmax=fac*srAve
     res=optim(par=c(srAve),fn=fitIt,method=c("Brent"),lower=c(srAve),upper=c(setmax),control=list(fnscale=-1,abstol=.Machine$double.eps))
     if(res$converge != 0) cat("fitIt did not converge!\n")

     ansSrMax[i] <- res$par
# ansConverge should be zero if all is well     
     ansConverge[i] <- res$convergence
#     ansMessage[i] <- res$message
     rPwr[i] <- res$value
     
# end sedimentation rate loop
  }
  if(verbose==1) close(progress) 
}


#######################################################################################
# (5)  Conduct grid search in parallel
#######################################################################################
if(ncores>1)
{
# LOAD libraries for parallel processing and set up parallel backend for your ncores.
# NOTE: the packages foreach and doParallel (or doSNOW) must be loaded if running this 
#  as a stand-alone script. uncomment the relevant lines below
#  library(foreach)
#  library(doParallel)
#  library(doSNOW)

# set up cluster
  cl<-makeCluster(as.integer(ncores))

# IMPORTANT: select the pacakge you will use for parallel processing
# uncommment below if using doParallel
  registerDoParallel(cl)
# uncomment below if using doSNOW
#  registerDoSNOW(cl)

  if(verbose==1) 
   {  
     cat("\n * PLEASE WAIT: Performing Optimization (this could take a while)\n")
# IMPORTANT: uncomment the next line ONLY if using doSNOW
#     cat("\n0%       25%       50%       75%       100%\n")
   }  
# create a progress bar
# IMPORTANT: uncomment the following three lines ONLY if you are using doSNOW
#  pbar <- txtProgressBar(max = numsed, style = 1, width=43)
#  progress <- function(n) setTxtProgressBar(pbar, n)
#  opts <- list(progress = progress)

# begin parallel simulation loop
# IMPORTANT: ucomment the following line if you are using doParallel
  resParallel<-foreach(ii=1:numsed,.combine=rbind) %dopar% {
# uncomment the following line if you are using doSNOW
#  resParallel<-foreach(ii=1:numsed,.options.snow = opts,.combine=rbind) %dopar% {
  
# NOTE: following line must be uncommented if not running this in astrochron
#    require("astrochron")
# NOTE: alternatively, you can use the following in foreach call above: .packages=c("astrochron")

#######################################################################################
# This is a duplicate of functions from section 2, without iteration updates 
#  (console output and plots).
#######################################################################################
genCycles <- function(timeModel, targetIn, n) 
  {
# generate cycle model    
    x <- matrix(0, n, 2*length(targetIn))
    for (i in 1:length(targetIn)) 
      {
        x[,2*i-1] <- cos( (2*pi)/(targetIn[i]) * timeModel[,1])
        x[,2*i] <- sin( (2*pi)/(targetIn[i]) * timeModel[,1])
      }
    return(x)
  }

# makeModel called from fitIt. note that srMax is always mapped to the maximum model value.
makeModel <- function(srMin,srMax)  
  {
    slope = (srMax-srMin)/(max(template[,2])-min(template[,2]))
    b = srMax-(slope*max(template[,2]))
    scaled[,2] = (slope*template[,2])+b
# sedrate2time is expecting sedrates in cm/ka
    scaled[,2] = scaled[,2] * 100     
    spaceTimeMap = sedrate2time(scaled,check=F,verbose=F,genplot=F)
    return(spaceTimeMap)
   }
   
testModel <- function(srMin)  
  {
    slope = (srMax-srMin)/(max(template[,2])-min(template[,2]))
    b = srMax-(slope*max(template[,2]))
    scaled[,2] = (slope*template[,2])+b
# sedrate2time is expecting sedrates in cm/ka
    scaled[,2] = scaled[,2] * 100     
    spaceTimeMap = sedrate2time(scaled,check=F,verbose=F,genplot=F)
    T2=spaceTimeMap[npts,2]
# experimenting here with different powers, to see how it influences tolerance for optimization       
    fit = (T2-totTime)^10
    return(fit)
  }

fitIt <- function(srMax) 
  {
# make srMax a global variable, so testModel can see it    
    srMax <<- srMax
# find optimal srMin given srMax and totTime. These sedrates are in m/ka; set maximum at srAve
    res1=optim(par=c(srAve),fn=testModel,method=c("Brent"),lower=difmin,upper=c(srAve),control=list(abstol=.Machine$double.eps))
    srMin=res1$par
    spaceTimeMap=makeModel(srMin,srMax)
    timeSeries=tune(dat,spaceTimeMap,check=F,genplot=F,verbose=F)
# missing one point after tuning?    
    timeSeries2=linterp(timeSeries,check=F,genplot=F,verbose=F)
    npts2=length(timeSeries2[,2])
# GENERATE TIME-CALIBRATED MODEL for spectral power assessment (applies for fit=1 or fit=2) 
# there is an advantage to constructing the model at the original times, then interpolating.
#  but we will simplify this for now...
    xm.0 <- genCycles(timeSeries2, targetTot, npts2)
# measure fit 
    lm.0 <- lm(timeSeries2[,2] ~ xm.0)
# bandpass precession or short eccentricity band
    bp = taner(timeSeries2,padfac=2,flow=flow,fhigh=fhigh,roll=roll,demean=T,detrend=F,addmean=F,check=F,genplot=F,verbose=F)
# hilbert transform for instantaneous amplitude
    hil = hilbert(bp,padfac=2,demean=T,detrend=F,check=F,addmean=F,genplot=F,verbose=F)
# GENERATE TIME-CALIBRATED MODEL for envelope assessment
    if(fit==1) xm.1 <- genCycles(timeSeries2, targetE, npts2)
    if(fit==2) xm.1 <- genCycles(timeSeries2, targetE[1], npts2)
# measure fit
    lm.1 <- lm(hil[,2] ~ xm.1)
    if(cormethod==1) 
     {
       rsq_pwr = cor(timeSeries2[,2],lm.0$fitted,method=c("pearson"))^2
       rsq_mod = cor(hil[,2],lm.1$fitted,method=c("pearson"))^2
     }  
    if(cormethod==2) 
     {
       rsq_pwr = cor(timeSeries2[,2],lm.0$fitted,method=c("spearman"))^2
       rsq_mod = cor(hil[,2],lm.1$fitted,method=c("spearman"))^2
     }  
    rsq_tot = rsq_pwr * rsq_mod
    
    if(iopt == 1) return(rsq_mod) 
    if(iopt == 2) return(rsq_pwr)
    if(iopt == 3) return(rsq_tot)
  }
##### PART ABOVE is copied from section 2 #############################################  

    srAve=sedrate[ii]
# convert average sedimentation rate to total time
    totTime = (max(dat[,1])-min(dat[,1]))/srAve
# execute functions
# note that fitIt is searching and optmizing on srMax, given srAve
#  optimal srMin is determined in testModel, given each tested srMax and totTime
# let's set maximum sedimentation rate to fac*srAve. default of 5 is based
#   on experimentation. If larger than this, risk getting into local minimum during fit.
# May 17, 2016: modified to allow for difmax; modified again Dec 21, 2017
     if(bound==1) setmax=difmax
     if(bound==2) setmax=fac*srAve
    res=optim(par=c(srAve),fn=fitIt,method=c("Brent"),lower=c(srAve),upper=c(setmax),control=list(fnscale=-1,abstol=.Machine$double.eps))

    return(c(res$par,res$convergence,res$value))
# end foreach loop
   }

# IMPORTANT: uncomment the following line (to close progress bar) ONLY if you are using doSNOW 
#  close(pbar)
# shut down the cluster
  stopCluster(cl)  

# parse results    
  ansSrMax <- resParallel[,1]
  ansConverge <- resParallel[,2]
  rPwr <- resParallel[,3]

  if(any(ansConverge!=0)) cat("fitIt did not converge for some results!\n")

# end ncores>1
}
       

#######################################################################################
# (6)  Postprocess
#######################################################################################
# find sedimentation rate with maxima in r*p
  ptRp=which.max(rPwr)
# check for multiple maxima
  sumMax = sum( rPwr == max(rPwr) )
  if(sumMax > 1 && verbose>0) cat("\n**** WARNING: Multiple maxima detected.\n") 

# recalculate srMin, also needed for bandpasses and Hilbert of optimal sedimentation rate for plotting
  srAve=sedrate[ptRp]   
  srMax=ansSrMax[ptRp]
  totTime = (max(dat[,1])-min(dat[,1]))/srAve
  res2=optim(par=c(srAve),fn=testModel,method=c("Brent"),lower=difmin,upper=c(srAve),control=list(abstol=.Machine$double.eps))
  srMin=res2$par

  if(verbose>0) 
   {
     cat(" * Maximum r-squared =", rPwr[ptRp],"at an average sedimentation rate of", sedrate[ptRp]*100,"cm/ka \n")
     cat("   Total duration =", totTime,"ka \n")
     cat("    srMin (min. sedrate) =", srMin*100, " and  srMax (max. sedrate) =",srMax*100,"\n")
   }  
  
  if(output == 2 || genplot>0)
   { 
# Plotting
     spaceTimeMap=makeModel(srMin,srMax)
     timeSeries=tune(dat,spaceTimeMap,check=F,genplot=F,verbose=F)
     timeSeries2=linterp(timeSeries,check=F,genplot=F,verbose=F)
# missing one point after tuning sometimes?   
     npts2=length(timeSeries2[,2])

# bandpass precession or short eccentricity band
     bp = taner(timeSeries2,padfac=2,flow=flow,fhigh=fhigh,roll=roll,demean=T,detrend=F,check=F,addmean=F,genplot=F,verbose=F)
# hilbert transform for instantaneous amplitude
     hil = hilbert(bp,padfac=2,demean=T,detrend=F,addmean=F,check=F,genplot=F,verbose=F)
# center precHil for plotting
#    hil[2] = hil[2] - colMeans(hil[2])
 
# generate modulation cycles at optimal sedimentation rate for plotting
     if(fit == 1) xm <- genCycles(timeSeries2, targetE, npts2)
     if(fit == 2) xm <- genCycles(timeSeries2, targetE[1], npts2)
     lm.3 <- lm(hil[,2] ~ xm) 
   }
 
 
#######################################################################################
# (7)  Summary plots
#######################################################################################
  if(genplot>0)
   { 
     dev.new()
     par(mfrow=c(3,2))
# plot least squares fitting results. Note, cex is the factor by which to increase or decrease default symbol size
     plot(sedrate*100,rPwr,cex=.75,col="red",xlab="Average Sedimentation Rate (cm/ka)",ylab="r-squared",main="Fit")
# sedrate plot
     slope = (srMax-srMin)/(max(template[,2])-min(template[,2]))
     b = srMax-(slope*max(template[,2]))
     scaled[,2] = (slope*template[,2])+b
# sedrate2time is expecting sedrates in cm/ka
     scaled[,2] = scaled[,2] * 100    
     plot(scaled,xlab="Depth (m)",ylab="Sedimentation rate (cm/ka)",main="Sedimentation Rates", col="red",type="l")
# plot bp of precession, hil and ecc estimated by least squares fitting
     plot(bp,col="blue",cex=.5,xlab="Time (ka)",ylab="Value",main="Hilbert vs. Bandpass",ylim=c(min(bp),max(hil[,2],lm.3$fitted)))
     lines(bp,col="blue")
     lines(hil,col="red")
     lines(hil[,1],-1*hil[,2],col="red")
     lines(hil[,1],lm.3$fitted,lwd=1.5)
     abline(h=0,col="black",lty=3)
# plot of tuned series
     plot(timeSeries2,type="l",xlab="Time (ka)",ylab="Value",main="Astronomically tuned record")
# cross plot of hil and ecc
     plot(hil[,2],lm.3$fitted, main="Hilbert vs. Optimal Eccentricity Fit", xlab="Hilbert",ylab="Fitted Eccentricity")
# periodogram of stratigraphic series
     fft = periodogram( timeSeries2, output=1, verbose=F,genplot=F)
# remove f(0) for log plot of power
     fft = subset(fft,(fft[,1] > 0))
     plot(fft[,1],fft[,3],xlim=c(0,0.1),type="l",xlab="Frequency (cycles/ka)",ylab="Power",main="Power Spectrum (black=linear; gray=log)")
# plot a second y-axis
     par(new=TRUE)
     plot(fft[,1],log(fft[,3]),xlim=c(0,0.1),type="l",yaxt="n",col="gray",xlab="",ylab="")
     abline(v=1/targetTot, col="red",lty=3)
     abline(v=flow, col="blue",lty=2)
     abline(v=fhigh, col="blue",lty=2)   
# end genplot
  }

#######################################################################################
# (8)  Output results
#######################################################################################
# return sedimentation rate grid, r*p  
  if(output == 1) 
   {
     out <- data.frame (cbind(100*sedrate,rPwr) )
     colnames(out) <- c("sedrate","r2")
     return(out)
   }
       
# return optimal time series, hilbert and fitted periods
  if(output == 2) 
   {
     out <- data.frame (cbind(timeSeries2[,1],timeSeries2[,2],bp[,2],hil[,2],lm.3$fitted) )
     if(fit==1) colnames(out) <- c("time","value","filtered_precession","precession_envelope","reconstructed_ecc_model")
     if(fit==2) colnames(out) <- c("time","value","filtered_short_ecc","short_ecc_envelope","reconstructed_ecc_model")
     return(out)
   }  

### END function timeOptTemplate
}
