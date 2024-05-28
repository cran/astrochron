### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2023 Stephen R. Meyers
###
###########################################################################
### function timeOptSimPwrLaw - (SRM: August 21, 2021; September 13, 2023)
###########################################################################

# modified from FUNCTION-timeOptSim_v9.R

timeOptSimPwrLaw <- function (dat,numsim=2000,beta=NULL,sedrate=NULL,sedmin=0.5,sedmax=5,numsed=100,linLog=1,limit=T,fit=1,r2max=1,fitModPwr=T,flow=NULL,fhigh=NULL,roll=NULL,targetE=NULL,targetP=NULL,detrend=T,ncores=2,output=0,genplot=T,check=T,verbose=T)
{

if(verbose) cat("\n----- TimeOpt Monte Carlo Simulation -----\n")

cormethod=1
if(!is.null(sedrate) && verbose) 
 {
   cat("\n**** WARNING: you are only investigating one sedimentation rate.\n")
   cat("       sedrate option takes precedence over sedmin/sedmax/numsed\n\n")
   sedmin=sedrate
   sedmax=sedrate
   numsed=1
 }

# prepare data array
   dat = data.frame(dat)      
   npts <- length(dat[,1]) 
   dx <- dat[2,1]-dat[1,1]

if(check)
 {
# error checking 
   if(dx<0)
     { 
       if (verbose) cat("\n * Sorting data into increasing height/depth/time, removing empty entries\n")
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

### what is the estimated beta (power law) coefficient?
   if(is.null(beta))
    {
      spec1=periodogram(dat,padfac=1,output=1,verbose=F,genplot=F)
      spec2=data.frame(cbind(spec1[,1],spec1[,3]))
# remove Nyquist    
      spec2=spec2[-length(spec2[,1]),]
      beta=pwrLawFit(spec2,output=2,verbose=F,genplot=F)[,1]
      if(verbose) cat(" * Estimated beta =",beta,"\n")
    }  

   res=timeOpt(dat,sedmin=sedmin,sedmax=sedmax,numsed=numsed,linLog=linLog,limit=limit,fit=fit,r2max=r2max,fitModPwr=fitModPwr,flow=flow,fhigh=fhigh,roll=roll,targetE=targetE,targetP=targetP,detrend=detrend,output=1,title=NULL,genplot=F,verbose=F,check=F)

   if (r2max==1) datCorMax = max(res[,4])
   if (r2max==2) datCorMax = max(res[,3])
   if (r2max==3) datCorMax = max(res[,2])
 
if(verbose)
 {  
  if (r2max==1) cat(" * (Envelope r^2) x (Spectral Power r^2) =", datCorMax,"\n")
  if (r2max==2) cat(" * (Spectral Power r^2) =", datCorMax,"\n")
  if (r2max==3) cat(" * (Envelope r^2) =", datCorMax,"\n")
 }
   
#######################################################################################
# Monte Carlo Simulation

  if(verbose) cat("\n * PLEASE WAIT: Performing", numsim,"simulations using", ncores,"cores\n")
# IMPORTANT: uncomment the next line ONLY if using doSNOW
#  if(verbose)  cat("\n0%       25%       50%       75%       100%\n")
 
# section below is parallelized
# LOAD libraries for parallel processing
# and set up parallel backend for your ncores
  if(ncores<2) stop("WARNING: number of ncores must be greater than one!")

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

# IMPORTANT: uncomment the following three lines ONLY if you are using doSNOW
#  pbar <- txtProgressBar(max = numsim, style = 1, width=43)
#  progress <- function(n) setTxtProgressBar(pbar, n)
#  opts <- list(progress = progress)

   
# begin parallel simulation loop
# IMPORTANT: ucomment the following line if you are using doParallel
  resParallel<-foreach(icount(numsim),.combine=rbind) %dopar% {
# uncomment the following line if you are using doSNOW
#  resParallel<-foreach(icount(numsim),.options.snow = opts,.combine=rbind) %dopar% {

# NOTE: following line must be uncommented if running this as a stand-alone script.
#  require("astrochron")
# NOTE: alternatively, you can use the following in foreach call above: .packages=c("astrochron")

# generate power law noise
  sim = pwrLaw(npts, dx, mean=0, sdev=1, beta=beta, genplot=F, verbose=F)
# recenter and standardize
  sim[2]=sim[2]-colMeans(sim[2])
  sim[2]=sim[2]/sapply(sim[2],sd)
               
  simres=timeOpt(sim,sedmin=sedmin,sedmax=sedmax,numsed=numsed,linLog=linLog,limit=limit,fit=fit,r2max=r2max,fitModPwr=fitModPwr,flow=flow,fhigh=fhigh,roll=roll,targetE=targetE,targetP=targetP,detrend=detrend,output=1,title=NULL,genplot=F,verbose=F,check=F)

  if(r2max==1) simMax = max(simres[,4])
  if(r2max==2) simMax = max(simres[,3])
  if(r2max==3) simMax = max(simres[,2])

# return results
  simMax
  
# end foreach loop
}
# shut down the cluster
  stopCluster(cl)  

# IMPORTANT: uncomment the following line (to close progress bar) ONLY if you are using doSNOW 
#  close(pbar)

  simres=resParallel

# now sort results, determine how many have values > your result
    numgt = sum(simres>datCorMax)
    pval=numgt/numsim 
    if(pval < (10/numsim) && (10/numsim) <=1 ) pval= 10/numsim
    if(pval >= (10/numsim) && (10/numsim) <=1 ) pval=pval   
    if((10/numsim) > 1 ) pval=1
     
    if(verbose) 
     {
       if (r2max==1) cat("\n * (Envelope r^2) * (Spectral Power r^2) p-value =",pval, "\n")
       if (r2max==2) cat("\n * (Spectral Power r^2) p-value =",pval, "\n")
       if (r2max==3) cat("\n * (Envelope r^2) p-value =",pval, "\n")
     }
     
    if(genplot)
     { 
      dev.new(title = paste("TimeOpt Monte Carlo Results"), height = 5, width = 6)
      par(mfrow=c(1,1))
      if (r2max==1) plot(density(simres), col="black",xlim=c(0,1),type="l",xlab=expression(paste({"r"^2}["opt"])),main=expression(bold(paste({"r"^2}["opt"]," Monte Carlo Results"))),cex.lab=1.1,lwd=2)
      if (r2max==2) plot(density(simres), col="black",xlim=c(0,1),type="l",xlab=expression(paste({"r"^2}["spectra"])),main=expression(bold(paste({"r"^2}["spectral"]," Monte Carlo Results"))),cex.lab=1.1,lwd=2)
      if (r2max==1) plot(density(simres), col="black",xlim=c(0,1),type="l",xlab=expression(paste({"r"^2}["envelope"])),main=expression(bold(paste({"r"^2}["envelope"]," Monte Carlo Results"))),cex.lab=1.1,lwd=2)
      polygon(density(simres),col="red",border=NA)
#      grid()
      abline(v=datCorMax,col="blue",lwd=2,lty=3)
      mtext(round(datCorMax,digits=5),side=3,line=0,at=datCorMax,cex=1,font=4,col="blue")
     }
     
# output = (0) nothing, (1) envelope*spectral power r^2 p-value, (2) output simulation r^2 results
     if(output == 1) return(data.frame(pval))
     if(output == 2) return(data.frame(simres))
     
### END function timeOptSimPwrLaw
}
