### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2024 Stephen R. Meyers
###
###########################################################################
### linterpLH13: linearly interpolate data using approach of Laepple and Huybers (2013)
###                                    (SRM: Nov 3-5, 2015; Aug 25-26, 2023; 
###                                     Sept 26-27, 2023; Sept 3, 2024; Sept 16, 2024;
###                                     Nov 10, 2024)
### (1) generate 1/f stochastic signal (Laepple and Huybers use beta=1)
### (2) subsample 1/f noise on dataset sampling grid
### (3) interpolate 1/f noise to finest sample step of dataset (default) or specified fine grid
### (4) Divide resampled power spectrum by theoretical spectrum and smooth
### (5) Highest reliable frequency (fr) determined when ratio crosses value of 0.7 (or thresh)
### (6) Optimal interpolation resolution is set using fr
### (7) interpolate data to 10 times optimal resolution
### (8) filter using lowpass filter with cutoff of 1.2/(optimal resolution)
### (9) resample at optimal resolution
### Option to conduct nsim simulations
###########################################################################

linterpLH13 <- function (dat,dt=NULL,start=NULL,dtMin=NULL,thresh=0.7,beta=1,tbw=20,smooth=0.05,logF=F,antialias=T,roll=10^10,nsim=1,ncores=2,output=T,maxF=NULL,pl=1,genplot=T,check=T,verbose=1)
{
  if(verbose==0) 
    {
      verbose1=F
      verbose2=F
    }  
  if(verbose==1) 
    {
      verbose1=T
      verbose2=F
    }    
  if(verbose==2) 
    {
      verbose1=T
      verbose2=T
    }    
  if(verbose1)
    {
      cat("\n----- APPLYING PIECEWISE-LINEAR INTERPOLATION TO STRATIGRAPHIC SERIES -----\n")
      cat("              USING THE APPROACH OF LAEPPLE AND HUYBERS (2013) \n")
    }  
  dat <- data.frame(dat)

  if(check)
   {
    dt2 <- dat[2,1]-dat[1,1]     
    if(dt2<0)
     { 
      if (verbose1) cat("\n * Sorting data into increasing height/depth/time, removing empty entries\n")
      dat <- dat[order(dat[1], na.last = NA, decreasing = F), ]
     }
   }
     
  npts <- length(dat[,1]) 
  if(verbose1) cat("\n * Number of samples=", npts,"\n")
 
  dt3=dat[2:npts,1]-dat[1:(npts-1),1] 
  if(check && length(unique(dt3)) == 1)
    {
      cat("\n**** ERROR: Data series is already on an even sample grid\n")
      stop("**** TERMINATING NOW!")
    }
 
### first define points for interpolation
  if(!is.null(start)) start2 <- start
  if(is.null(start)) start2 <- dat[1,1]
  
  if(verbose1) cat(" * Minimum, Median, Maximum sampling interval =",min(dt3),median(dt3),max(dt3),"\n")

# find optimal dt
   if(is.null(dt))
     {
       if(verbose1) cat("\n * Determining optimal sampling interval for series using",nsim,"simulation(s)\n")
# determine length and minimum sampling interval
       duration <- abs(dat[npts,1]-dat[1,1])
       if(is.null(dtMin)) dtMin=min(dt3)
    
       if(nsim==1)
         {
# (1) generate 1/f stochastic signal (Laepple and Huybers use beta=1)
#     add extra data point to ensure that you capture the entire interval of dat
           npts2=as.integer(duration/dtMin) + 2 
           noise=pwrLaw(dt=dtMin,npts=npts2,beta=beta,genplot=F,verbose=verbose2)
           noise[1]=noise[1]+min(dat[,1])
# (2) subsample 1/f noise on dataset sampling grid
           noiseR=resample(noise,dat[,1],genplot=F,check=F,verbose=verbose2)
# (3) interpolate 1/f noise to finest sample step of dataset
           noiseR2=linterp(noiseR,dt=dtMin,genplot=F,verbose=verbose2)
# (4) Divide resampled power spectrum by theoretical spectrum
# first full noise spectrum (before resampling)
           if(verbose2) cat("\n----- PERFORMING Multitaper Spectral Analysis -----\n")
           mtmNoise=mtm(noise,padfac=1,tbw=tbw,output=1,genplot=F,verbose=F)
           pwrNoise=data.frame(cbind(mtmNoise[,1],mtmNoise[,2]))
# noise spectrum for resampled and interpolated 1/f noise
           mtmNoiseR=mtm(noiseR2,padfac=1,tbw=tbw,output=1,genplot=F,verbose=F)
           pwrNoiseR=data.frame(cbind(mtmNoiseR[,1],mtmNoiseR[,2]))
# resample theoretical signal on frequency grid from noiseR, for ratio estimate
           theoryPwr=resample(pwrNoise,pwrNoiseR[,1],genplot=F,check=F,verbose=verbose2)
# calculate ratio
           ratio=pwrNoiseR
           ratio[2]=pwrNoiseR[2]/theoryPwr[2]
           if(verbose1 && max(ratio[2]) > 1) cat("\n * NOTE: Maximum in R is",max(ratio[2]),"\n")
# perform Gaussian kernel smoothing on ratio (using frequency or log10-frequency)
           ratio2=ratio
           if(logF) ratio2[1]=log10(ratio[1])
           if(smooth>0) smoothRatio=noKernel(ratio2,smooth=smooth,output=2,genplot=F,verbose=verbose2) 
           if(smooth==0) smoothRatio=ratio2
           if(logF) smoothRatio[1]=10^smoothRatio[1]
           if(verbose1 && max(smoothRatio[2]) > 1) cat(" * NOTE: Maximum in R2 is",max(smoothRatio[2]),"\n")
# (5) Highest reliable frequency (fr) determined when ratio crosses value of 0.7 or thresh
           ii=min(which(smoothRatio[2]<thresh))
           fr=smoothRatio[ii,1]
# (6) Optimal interpolation resolution is set using fr
           optDt=0.5/fr  
# end nsim==1
         }   

       if(nsim>1)
         {
# NOTE: the packages foreach and doParallel must be loaded if running this 
#  as a stand-alone script. uncomment the relevant lines below
#           library(foreach)
#           library(doParallel)
           
# set up cluster
           cl<-makeCluster(as.integer(ncores))
           registerDoParallel(cl)
   
# begin parallel simulation loop
           resParallel<-foreach(icount(nsim),.combine=rbind) %dopar% {

# NOTE: following line must be uncommented if running this as a stand-alone script.
#           require("astrochron")
# NOTE: alternatively, you can use the following in foreach call above: .packages=c("astrochron")

# turn off verbose output for multiple simulations
           verboseSwitch=F
           if (verbose2) { verboseSwitch = T ; verbose2 = F }
         
# (1) generate 1/f stochastic signal (Laepple and Huybers use beta=1)
#     add extra data point to ensure that you capture the entire interval of dat
           npts2=as.integer(duration/dtMin) + 2 
           noise=pwrLaw(dt=dtMin,npts=npts2,beta=beta,genplot=F,verbose=verbose2)
           noise[1]=noise[1]+min(dat[,1])
# (2) subsample 1/f noise on dataset sampling grid
           noiseR=resample(noise,dat[,1],genplot=F,check=F,verbose=verbose2)
# (3) interpolate 1/f noise to finest sample step of dataset
           noiseR2=linterp(noiseR,dt=dtMin,genplot=F,verbose=verbose2)
# (4) Divide resampled power spectrum by theoretical spectrum
# first full noise spectrum (before resampling)
           mtmNoise=mtm(noise,padfac=1,tbw=tbw,output=1,genplot=F,verbose=F)
           pwrNoise=data.frame(cbind(mtmNoise[,1],mtmNoise[,2]))
# noise spectrum for resampled and interpolated 1/f noise
           mtmNoiseR=mtm(noiseR2,padfac=1,tbw=tbw,output=1,genplot=F,verbose=F)
           pwrNoiseR=data.frame(cbind(mtmNoiseR[,1],mtmNoiseR[,2]))
# resample theoretical signal on frequency grid from noiseR, for ratio estimate
           theoryPwr=resample(pwrNoise,pwrNoiseR[,1],genplot=F,check=F,verbose=verbose2)
# calculate ratio
           ratio=pwrNoiseR
           ratio[2]=pwrNoiseR[2]/theoryPwr[2]
# perform Gaussian kernel smoothing on ratio (using frequency or log10-frequency)
           ratio2=ratio
           if(logF) ratio2[1]=log10(ratio[1])
           if(smooth>0) smoothRatio=noKernel(ratio2,smooth=smooth,output=2,genplot=F,verbose=verbose2) 
           if(smooth==0) smoothRatio=ratio2 
           if(logF) smoothRatio[1]=10^smoothRatio[1]
# (5) Highest reliable frequency (fr) determined when ratio crosses value of 0.7 or thresh
           ii=min(which(smoothRatio[2]<thresh))
           frN=smoothRatio[ii,1]
           if(verboseSwitch) verbose2 = T
# return results
           frN
# end foreach loop
           }
# shut down the cluster
           stopCluster(cl)  
# (6) Optimal interpolation resolution is set using fr
           frN=resParallel
           optDtN=0.5/frN
           optDt=median(optDtN)
           fr=0.5/optDt
           if(verbose1) 
             {
               cat("\n * Median simulated optimal cutoff frequency =",fr,"\n")
               cat(" * Median simulated optimal interpolation sampling interval =",optDt,"\n")
             }  
# end nsim>1
         }

       if(verbose1 && nsim==1) 
          {
            cat("\n * Optimal cutoff frequency =",fr,"\n")
            cat(" * Optimal interpolation sampling interval =",optDt,"\n")
          }  

       if(genplot)
         {
           dev.new()
           par(mfrow=c(3,2))
           if(nsim == 1)
             { 
               if(is.null(maxF)) maxF=max(ratio[,1])
               plot(ratio,col="blue",cex=0.5,ylab="Spectral Ratio",xlab="Frequency",main="Sampled Spectrum/Theoretical Spectrum",xlim=c(0,maxF))
               lines(smoothRatio,lwd=2,col="black")
               abline(h=thresh,lty=3)
               abline(v=fr,lty=4)
               abline(v=1/(2*median(dt3)),lty=3,col="blue")
               if(pl==1) plLog="y"
               if(pl==2) plLog=""
               plot(pwrNoiseR,type="l",log=plLog,ylab="Spectral Power",xlab="Frequency",main="Sampled (black), Theory (red), Interpolated (green)",xlim=c(0,maxF),ylim=c(min(pwrNoiseR[,2],theoryPwr[,2]),max(pwrNoiseR[,2],theoryPwr[,2])))
               lines(theoryPwr,col="#FF00003C",lwd=2)
               abline(v=fr,lty=4)
             }  
           if(nsim > 1)
             { 
               hist(optDtN,freq=F,xlab="Sample spacing",main="Distribution of simulated optimal dt"); lines(density(optDtN, bw="nrd"),col="red");grid()
             }  
         }
# end optimal dt section
     } 

# error checking 
  if(!is.null(dt))
    {
      if(genplot) 
        {
          dev.new()
          par(mfrow=c(2,2))
        }   
      if(dt==0) 
        {
          if (verbose1) cat("\n**** ERROR: dt is ZERO!\n")
          stop("**** TERMINATING NOW!")
        }  
      if(dt<0) 
        {
          if (verbose1) cat("\n**** WARNING: dt must be positive. Will take absolute value.\n")
          dt=abs(dt)
        }  
      optDt=dt  
    }  

# note that due to the anti-alias filtering approach, interpolated data values that fall
#  on the original sampling grid typically do not exactly match the original values.
#  this is a consequence of the fact that some of the data variance has been removed 
#  by the filtering.
#
#  but importantly, the first and last value of the interpolated series have a tendency to
#  be strongly biased when these points are close to the original sampling grid. this can
#  be demonstrated with the following code (run it 10+ times to evaluate, but first
#  comment the linterpLH13 section that removes points):
#  ex = ar1(npts=501,dt=10,genplot=F);ex2=delPts(ex,sample.int(501,80),genplot=F);res=linterpLH13(ex2,dt=10,genplot=F);pl(2);plot(res,xlim=c(0,200));points(ex2,col="red");plot(res,xlim=c(4800,5000));points(ex2,col="red")
#  to address this issue, the first and last interpolated points are removed if they 
#  are within 0.5*optDt (the optimal interpolation interval) of the original sampling grid.

  if(!antialias) low=dat
  if(antialias)
    {
# (7) interpolate data to 10 times optimal resolution
      dat2=linterp(dat,dt=0.1*optDt,genplot=F,check=F,verbose=verbose2)
# (8) filter using lowpass filter with cutoff of 1.2/(optimal resolution)
      low=taner(dat2,fhigh=1.2/optDt,roll=roll,genplot=F,check=F,verbose=verbose2)
    }

# (9) resample at optimal resolution
  res=linterp(low,dt=optDt,start=start2,genplot=F,check=F,verbose=verbose2)
  
  if(antialias)
    {
# remove first and last interpolated data point if they are within 0.5 optDt 
# (interpolation interval) of the original sampling grid
      if(abs(res[1,1]-dat[1,1]) < 0.5*optDt) 
        {
          res = res[-1,] 
          if (verbose1) cat("\n * First interpolated value removed due to potential bias.\n")
        }  
      if(abs(res[length(res[,1]),1]-dat[npts,1]) < 0.5*optDt) 
        {
          res = res[-length(res[,1]),]
          if (verbose1) cat("\n * Last interpolated value removed due to potential bias.\n") 
        }  
    }
    
  if(genplot && is.null(dt) && nsim==1)
    {     
# do the same thing to the resampled noise
       if(!antialias) low2=noiseR
       if(antialias)
         {
           noiseFinal=linterp(noiseR,dt=0.1*optDt,genplot=F,check=F,verbose=verbose2)
           low2=taner(noiseFinal,fhigh=1.2/optDt,genplot=F,check=F,verbose=verbose2)
         }
       resNoise=linterp(low2,dt=optDt,genplot=F,check=F,verbose=verbose2)   
       if(abs(resNoise[1,1]-noiseR[1,1]) < 0.5*optDt)  resNoise = resNoise[-1,] 
       if(abs(resNoise[length(resNoise[,1]),1]-noiseR[npts,1]) < 0.5*optDt) resNoise = resNoise[-length(resNoise[,1]),]
       mtmRes=mtm(resNoise,padfac=1,tbw=tbw,output=1,genplot=F,verbose=F)
       pwrRes=data.frame(cbind(mtmRes[,1],mtmRes[,2]))      
       lines(pwrRes,col="#0064003C",lwd=3)
    }

### assign interpolated data to data.frame d
  d <- as.data.frame(res)
  npts3=length(d[,1])
  if(verbose) cat("\n * New number of samples=",npts3,"\n")
  colnames(d)[1] <- colnames(dat[1])
  colnames(d)[2] <- colnames(dat[2])

### plot data
  if(genplot)
    {
      plot(d,cex=0.5, xlab="Location",ylab=colnames(dat[2]),main="Raw (black) and Interpolated (red) Data",col="red"); lines(dat)
### plot the denisty and the histogram together
      hist(d[,2],freq=F,xlab="Interpolated Value",main="Distribution of Interpolated Values"); lines(density(d[,2], bw="nrd"),col="red");grid()
### boxplot
      boxplot(d[,2],ylab="Interpolated Value",main="Boxplot of Interpolated Values")
### Normal probabilty plot (Normal Q-Q Plot)
      qqnorm(d[,2]); qqline(d[,2], col="red"); grid()
    } 
    
  if(output) return(d)

### END function linterpLH13
}
