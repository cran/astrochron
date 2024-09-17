### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2018 Stephen R. Meyers
###
###########################################################################
### timeOptTemplatePlot: generate summary figure for timeOptTemplate 
###                      (SRM: January 3-4, 2018; Jan 8, 2018; October 23, 2018;
###                            November 22, 2018; November 24, 2018; September 10, 2024)
###
###########################################################################

timeOptTemplatePlot <- function (dat=NULL,template=NULL,detrend=T,detrendTemplate=F,flipTemplate=F,srMin=NULL,srMax=NULL,res1=NULL,simres=NULL,fit=1,flow=NULL,fhigh=NULL,roll=NULL,targetE=NULL,targetP=NULL,xlab="Depth (m)",ylab="Proxy Value",fitR=NULL,output=0,verbose=T)
{

if(verbose) cat("\n----- Generating summary plot for TimeOptTemplate analysis----\n")

# error checking (ADD MORE ERROR CHECKS HERE!)
if(is.null(dat)) stop("**** ERROR: missing dat. Terminating now!")
if(is.null(res1)) stop("**** ERROR: missing res1. Terminating now!")
if(is.null(simres)) stop("**** ERROR: missing simres. Terminating now!")
if(is.null(targetE)) stop("**** ERROR: missing targetE. Terminating now!")
if(is.null(targetP) && fit==1) cat("**** WARNING: missing targetP.")

# Prepare data array
  dat = data.frame(dat)      
  npts = length(dat[,1]) 
  dx = dat[2,1]-dat[1,1]
  
# save copy for plotting
  dat2=dat
    
# convert from cm/ka to m/ka for processing
  srMin=srMin/100
  srMax=srMax/100
  
# convert average sedimentation rate to total time
#  totTime = (max(dat[,1])-min(dat[,1]))/srAve  

plf=T
if(is.null(flow) || is.null(fhigh) || is.null(roll)) 
 {
   if(verbose) cat("**** WARNING: flow, fhigh and/or roll not specified, will omit from plot")
   plf=F
 }  

if(is.null(fitR) && verbose) cat("**** WARNING: fitR not specified, will omit from plot")
 
# detrend data series
  if(detrend) 
    {
      lm.1 <- lm(dat[,2] ~ dat[,1])
      dat[2] <- dat[2] - (lm.1$coeff[2]*dat[1] + lm.1$coeff[1])
      if(verbose>0) cat(" * Linear trend subtracted from data series. m=",lm.1$coeff[2],"b=",lm.1$coeff[1],"\n")
    }

# standardize data series
  dat[2]=dat[2]-colMeans(dat[2])
  dat[2]=dat[2]/sapply(dat[2],sd)

  if(flipTemplate)
    {
      template[2]=template[2]*-1
      if(verbose>0) cat("\n * Sedimentation template flipped.\n")
    } 

# detrend template
  if(detrendTemplate) 
    {
      lm.2 <- lm(template[,2] ~ template[,1])
      template[2] <- template[2] - (lm.2$coeff[2]*template[1] + lm.2$coeff[1])
      if(verbose>0) cat(" * Linear trend subtracted from sedimentation template. m=",lm.2$coeff[2],"b=",lm.2$coeff[1],"\n")
    }
 
# make a copy of template for makeModel
  scaled <- template
 
 
# determine sedimentation rate history, given average sedrate and srMax 
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
   
makeModel2 <- function(srMin,srMax)  
  {
    slope = (srMax-srMin)/(max(template[,2])-min(template[,2]))
    b = srMax-(slope*max(template[,2]))
    scaled[,2] = (slope*template[,2])+b
# sedrate2time is expecting sedrates in cm/ka
    scaled[,2] = scaled[,2] * 100     
    return(scaled)
   }   
   
# determine sedimentation rate history, make spaceTimeMap, tune series
    sedrates=makeModel2(srMin,srMax)
    spaceTimeMap=makeModel(srMin,srMax)
    timeSeries=tune(dat,spaceTimeMap,check=F,genplot=F,verbose=F)
    timeSeries2=linterp(timeSeries,check=F,genplot=F,verbose=F)
# missing one point after tuning sometimes?    
    npts2=length(timeSeries2[,2]) 

# also tune dat2 
    dat2tuned=tune(dat2,spaceTimeMap,check=F,genplot=F,verbose=F)
 
# determine res2 by making call to timeOpt
#  do not detrend again
res2=timeOpt(timeSeries2,detrend=FALSE,fit=fit,flow=flow,fhigh=fhigh,roll=roll,targetE=targetE,targetP=targetP,numsed=1,sedmin=100,sedmax=100,output=2,genplot=F)
 
 
# plot 
pl(r=4,c=2,mar=c(3,3,1,3))

# PLOT 1: (A) Stratigraphic Series
plot(dat2[,1], dat2[,2], col = "black", xlab = "", ylab = "", main = "",type="l",lwd=1.2)
mtext(ylab,side=2,line=1.8,cex=0.8)            
mtext(xlab,side=1,line=2,cex=0.8)   

            
# PLOT 2: (E) timeOpt fit
plot(res1[, 1], res1[, 2], type="l", lwd=2, col = "black", xlab = "", ylab = "", main = "")
mtext("Envelope & Spectral Fit",side=2,line=2,cex=0.8)            
mtext("Average Sedimentation Rate (cm/ka)",side=1,line=2,cex=0.8)   


# PLOT 3: (B) Sedimentation Rates
plot(sedrates[,1], sedrates[,2], col = "black", xlab = "", ylab = "", main = "",type="l",lwd=1.2)
mtext("Sed. Rate (cm/ka)",side=2,line=1.8,cex=0.8)            
mtext(xlab,side=1,line=2,cex=0.8)   


# PLOT 4: (F) Monte Carlo
denPlot=density(simres[,1])
yRange=max(denPlot$y)-min(denPlot$y)
maxR=max(max(denPlot$x),fitR)
xlim2=c(0,maxR+maxR*0.5)
if(max(xlim2) > 1) xlim2=c(0,1)
plot(denPlot,xlim=xlim2,main="",xlab="")
polygon(denPlot,col="black")
mtext("# Simulations",side=2,line=2,cex=0.8)            
mtext("Envelope & Spectral Fit",side=1,line=2,cex=0.8)   
if(!is.null(fitR))
 {
   abline(v=fitR,lwd=2,lty=4,col="red")
   text(fitR,max(denPlot$y)-0.4*yRange,expression(paste("r"^"2",""[opt],"=")),cex=1.3,font=2,pos=4) 
   text(fitR,max(denPlot$y)-0.52*yRange,round(fitR,digits=3),cex=1.2,font=2,pos=4) 
 }


# PLOT 5: (C) Tuned Series
plot(dat2tuned[,1], dat2tuned[,2], col = "black", xlab = "", ylab = "", main = "",type="l",lwd=1.2)
mtext(ylab,side=2,line=1.8,cex=0.8)            
mtext("Elapsed Time (ka)",side=1,line=2,cex=0.8)   


# PLOT 6: (G) Tuned spectrum
fft = periodogram(data.frame(cbind(res2[, 1], res2[, 2])), output = 1, verbose = F, genplot = F)
fft = subset(fft, (fft[, 1] > 0))
xlim2=c(0, min(0.1,max(fft[,1])))
plot(fft[, 1], fft[, 3], xlim = xlim2, type = "l",lwd=0, xlab = "", ylab = "", main = "")
if(plf) 
 {
     bpWin= taner(data.frame(cbind(res2[, 1], res2[, 2])),padfac=2,flow=flow,fhigh=fhigh,roll=roll,demean=T,detrend=F,addmean=F,genplot=F,verbose=F,output=2)
     polygon(bpWin[,1],bpWin[,2]*max(fft[,3]),col="#FFFF005A",border=NA)
 } 
par(new = TRUE)
plot(fft[, 1], fft[, 3], xlim = xlim2, type = "l",lwd=1.5, xlab = "", ylab = "", main = "",axes=F)
par(new = TRUE)
plot(fft[, 1], log(fft[, 3]), xlim = xlim2,type = "l",lwd=1.5, yaxt = "n", col = "gray", xlab = "", ylab = "")
if(fit==1) targetTot=append(targetE,targetP)
if(fit==2 && is.null(targetP)) targetTot=append(targetE,targetP)
if(fit==2 && !is.null(targetP)) targetTot=append(targetE)
abline(v = 1/targetTot, col = "red", lty = 3,lwd=1.5)
mtext("Spectral Power",side=2,line=2,cex=0.8)            
mtext("Frequency (cycles/ka)",side=1,line=2,cex=0.8)   


# plot 7: plot bp of precession, hil and ecc estimated by least squares fitting
plot(res2[,1],res2[,3],type="l",col="blue",cex=.5,xlab="",ylab="",main="",ylim=c(min(res2[,3]),max(res2[,4],res2[,5])))
lines(res2[,1],res2[,4],col="red")
lines(res2[,1],-1*res2[,4],col="red")
lines(res2[,1],res2[,5],lwd=1.5)
abline(h=0,col="black",lty=3)
mtext("Standardized Value",side=2,line=1.8,cex=0.8)            
mtext("Elapsed Time (ka)",side=1,line=2,cex=0.8)   


# plot 8: (H)
plot(res2[,4], res2[,5], main = "", xlab = "", ylab = "")
mtext("Fitted Eccentricity",side=2,line=2,cex=0.8)            
if(fit==1) mtext("Data Precession Envelope",side=1,line=2,cex=0.8)   
if(fit==2) mtext("Data Short Eccentricity Envelope",side=1,line=2,cex=0.8)   
abline(0,1,lty=2,lwd=2,col="red")

if(output==1) return(res2)
if(output==2) return(spaceTimeMap)
if(output==3) return(sedrates)

}