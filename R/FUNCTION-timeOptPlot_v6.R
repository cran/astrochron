### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2018 Stephen R. Meyers
###
###########################################################################
### timeOptPlot: generate summary figure for timeOpt (SRM: April 23, 2017;
###                      July 18-19, 2017; August 13, 2017; May 10, 2018;
###                      October 23, 2018)
###
###########################################################################

timeOptPlot <- function (dat=NULL,res1=NULL,res2=NULL,simres=NULL,fit=1,fitModPwr=T,flow=NULL,fhigh=NULL,roll=NULL,targetE=NULL,targetP=NULL,xlab="Depth (m)",ylab="Proxy Value",fitR=NULL,verbose=T)
{

if(verbose) cat("\n----- Generating summary plot for TimeOpt analysis----\n")

# error checking
if(is.null(dat)) stop("**** ERROR: missing dat. Terminating now!")
if(is.null(res1)) stop("**** ERROR: missing res1. Terminating now!")
if(is.null(res2)) stop("**** ERROR: missing res2. Terminating now!")
if(is.null(simres)) stop("**** ERROR: missing simres. Terminating now!")
if(is.null(targetE)) stop("**** ERROR: missing targetE. Terminating now!")
if(is.null(targetP) && fit==1) cat("**** WARNING: missing targetP.")

plf=T
if(is.null(flow) || is.null(fhigh) || is.null(roll)) 
 {
   if(verbose) cat("**** WARNING: flow, fhigh and/or roll not specified, will omit from plot")
   plf=F
 }  

if(is.null(fitR) && verbose) cat("**** WARNING: fitR not specified, will omit from plot")
 
pl(r=4,c=2,mar=c(3,3,1,3))
# plot 1: (A)
plot(dat[,1], dat[,2], col = "black", xlab = "", ylab = "", main = "",type="l",lwd=1.2)
mtext(ylab,side=2,line=1.8,cex=0.8)            
mtext(xlab,side=1,line=2,cex=0.8)   
            
# plot 2: (E)
plot(res1[, 1], res1[, 2], cex = 1.2, col = "red", xlab = "", ylab = "", main = "",pch=16)
par(new = T)
plot(res1[, 1], res1[, 3], col = "#00000096", xlab = "", ylab = "", type = "l", axes = F, lwd = 2)
axis(4, ylim = c(0, max(res1[, 3])), lwd = 1, col = "black")
mtext("Envelope Fit",side=2,line=2,cex=0.8,col="red")            
mtext("Sedimentation Rate (cm/ka)",side=1,line=2,cex=0.8)   
mtext("Spectral Power Fit",side=4,line=2,cex=0.8)  

# plot 3: (B)
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
if(fit == 1 && fitModPwr) targetTot = c(targetE,targetP)
if(fit == 1 && !fitModPwr) targetTot = c(targetP)
if(fit == 2 && fitModPwr) targetTot = c(targetE)
if(fit == 2 && !fitModPwr) targetTot = c(targetE[-1])
abline(v = 1/targetTot, col = "red", lty = 3,lwd=1.5)
mtext("Spectral Power",side=2,line=2,cex=0.8)            
mtext("Frequency (cycles/ka)",side=1,line=2,cex=0.8)   

# plot 4: (F)
plot(res1[, 1], res1[, 4], type="l", lwd=2, col = "black", xlab = "", ylab = "", main = "")
mtext("Envelope & Spectral Fit",side=2,line=2,cex=0.8)            
mtext("Sedimentation Rate (cm/ka)",side=1,line=2,cex=0.8)   

# plot 5: (C)   
ylim=c(min(res2[,3],res2[,4],-1*res2[,4]),max(res2[,3],res2[,4],-1*res2[,4])) 
plot(res2[,1], res2[,3], col = "black", cex = 0.5, xlab = ")", ylab = "", main = "",type="l",ylim=ylim)
lines(res2[,1], res2[,4], col = "red",lwd=2)
lines(res2[,1], -1 * res2[,4], col = "red",lwd=2)
abline(h = 0, col = "black", lty = 3)
mtext("Standardized Value",side=2,line=2,cex=0.8)            
mtext("Time (ka)",side=1,line=2,cex=0.8)   

# plot 6: (G)
denPlot=density(simres[,1])
yRange=max(denPlot$y)-min(denPlot$y)
maxR=max(max(denPlot$x),fitR)
xlim2=c(0,maxR+maxR*0.5)
if(max(xlim2) > 1) xlim2=c(0,1)
plot(denPlot,xlim=xlim2,main="")
polygon(denPlot,col="black")
mtext("# Simulations",side=2,line=2,cex=0.8)            
mtext("Envelope & Spectral Fit",side=1,line=2,cex=0.8)   
if(!is.null(fitR))
 {
   abline(v=fitR,lwd=2,lty=4,col="red")
   text(fitR,max(denPlot$y)-0.4*yRange,expression(paste("r"^"2",""[opt],"=")),cex=1.3,font=2,pos=4) 
   text(fitR,max(denPlot$y)-0.52*yRange,round(fitR,digits=3),cex=1.2,font=2,pos=4) 
 }
 
# plot 7: (D)
ylim=c(min(res2[,4],res2[,5]),max(res2[,4],res2[,5])) 
plot(res2[,1],res2[,4], cex = 0.5, xlab = "Time (ka)", ylab = "Value",main = "", col = "red", type = "l",lwd=2,ylim=ylim)
lines(res2[,1], res2[,5],lwd=2)
mtext("Standardized Value",side=2,line=2,cex=0.8)            
mtext("Time (ka)",side=1,line=2,cex=0.8)   

# plot 8: (H)
plot(res2[,4], res2[,5], main = "", xlab = "", ylab = "")
mtext("Fitted Eccentricity",side=2,line=2,cex=0.8)            
if(fit==1) mtext("Data Precession Envelope",side=1,line=2,cex=0.8)   
if(fit==2) mtext("Data Short Eccentricity Envelope",side=1,line=2,cex=0.8)   
abline(0,1,lty=2,lwd=2,col="red")

}