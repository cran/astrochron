### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2023 Stephen R. Meyers
###
###########################################################################
### function insoSeries : wrapper function to calculate insolation using 
###                     palinsol package, following astrochron syntax 
###                     (SRM: October 9, 2020; March 30, 2022; November 8, 2022;
###                           July 22, 2023)
###
###########################################################################


insoSeries <- function (tmin=0,tmax=1000,dt=1,opt=1,long=90,lat=65,threshold=400,l1=0,l2=70,S0=1365,genplot=TRUE,verbose=TRUE) 
 {

  if(verbose) 
   {
    cat("\n----- CALCULATE INSOLATION USING LASKAR ET AL. (2004) SOLUTION -----\n")
    if(opt==1)  cat("       Calculating Insolation for Given Day and Latitude\n")
    if(opt==2)  cat("       Calculating Caloric Insolation for Given Latitude\n")
    if(opt==3)  cat("       Calculating Integrated Insolation for Given Latitude\n")
    if(opt==4)  cat("       Calculating Integrated Insolation for Given Latitude, with Threshold\n")
   }

# longitude values: pi/2 = June solstice; 3*pi/2 = December solstice
  
# error checking
  if(lat > 90 || lat < -90)
   {
     cat("\n**** ERROR: the assigned value for lat is not valid. -90 <= lat <= 90 \n")
     stop("**** TERMINATING NOW!")
   }     

# error checking
  if(long > 360 || long < 0)
   {
     cat("\n**** ERROR: the assigned value for long is not valid. 0 <= long <= 360 \n")
     stop("**** TERMINATING NOW!")
   }    

  if(tmin>tmax) 
   {
    if(verbose) cat("\n**** WARNING: tmin is greater than tmax. The values will be switched.\n")
    tmax2=tmax
    tmin2=tmin
    tmin=tmax2
    tmax=tmin2
   } 
  
# restrict to 21 million years into future (>= -21000). 
  if(tmin < (-21000))
   {
    if(verbose) cat("\n**** WARNING: This solution can only be calculated to 21000 kiloyears in the future. Resetting tmin to -20000.\n")
    tmin=-21000
   }   

# restrict to 51 million years into past (<= 51000).
  if(tmax > 51000)
   {
    if(verbose) cat("\n**** WARNING: This solution can only be calculated to 51000 kiloyears in the past. Resetting tmax to 50000.\n")
    tmax=51000
   }    

  if(dt<1)
   {
    if(verbose) cat("\n**** WARNING: dt is less than 1 ka. It will be increased to 1 ka.\n")
    dt = 1
   } 


# convert longitude and latitude from degrees to radians 
  longRad=long*pi/180
  latRad=lat*pi/180
  l1Rad=l1*pi/180
  l2Rad=l2*pi/180
 
  times <- seq(from=-1000*tmax,to=-1000*tmin,by=dt*1000)

# Calculating Mean Insolation for Given Day and Latitude, using function 'Insol'
  if(opt==1) res=sapply(times, function(tt) Insol(orbit=la04(tt),long=longRad,lat=latRad,S0=S0))
# Calculating Caloric Insolation for Given Latitude, using function 'calins'
  if(opt==2) res=sapply(times, function(tt) calins(orbit=la04(tt),lat=latRad,S0=S0))
# Calculating Integrated Insolation for Given Latitude, using function 'Insol_l1l2'
  if(opt==3) res=sapply(times, function(tt) Insol_l1l2(orbit=la04(tt),lat=latRad,l1=l1Rad,l2=l2Rad,S0=S0))
# Calculating Integrated Insolation for Given Latitude, with Threshold, using function 'thrins'
  if(opt==4) res=sapply(times, function(tt) thrins(orbit=la04(tt),lat=latRad,threshold=threshold))
  
  out=data.frame(cbind(times/-1000,res))
  colnames(out) = c("Time_ka","Insol")
# sort into increasing order for output
  out <- out[order(out[,1], na.last = NA, decreasing = FALSE),]

if(genplot)
 {
   par(mfrow = c(2, 1))
   plot(out, type="l",col="red", main="", xlab="", ylab="",bty="n")
   mtext(paste("Insolation at", lat, "degrees"), side=3,line=1.5,font=2,cex=1.5)
   mtext(expression(Insolation (W/m^2)),side=2,line=2,cex=1.3)
   mtext("Time (kyr)",side=1,line=2,cex=1.3)
   
   insoSpec=periodogram(out,detrend=TRUE,output=1,genplot=FALSE,verbose=FALSE)
   plot(insoSpec[,1],insoSpec[,2],xlim=c(0,0.1),type="l",xlab="",ylab="")
   mtext("Insolation Spectrum",side=3,line=2,font=2,cex=1.5)
   mtext("Frequency",side=1,line=2.5,cex=1.3)
   mtext("Amplitude",side=2,line=2.5,cex=1.3)   
 }
 
  return(out)

### END function insoSeries
}

