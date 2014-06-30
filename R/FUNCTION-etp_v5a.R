### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2014 Stephen R. Meyers
###
###########################################################################
### etp: generate combined eccentricity-tilt-precession curve 
###          (SRM: March 23, 2012; May 3, 2012; October 11, 2012; April 23, 2013; April 29, 2013
###                May 20, 2013; July 10, 2013; August 15, 2013; April 7, 2014; June 26-27, 2014)
###
###########################################################################

etp <- function (laskar=NULL,tmin=0,tmax=1000,dt=1,eWt=1,oWt=1,pWt=1,sol=5,esinw=F,standardize=T,genplot=T,verbose=T)
{

  if(verbose) cat("\n----- GENERATING ECCENTRICITY-TILT-PRECESSION SERIES -----\n")
  if(is.null(laskar))      
   {
     cat("\n**** ERROR: You must input the Laskar astronomical solutions\n")
     stop("    TERMINATING NOW!")
   }  

  if(sol == 1 ) 
    {
      la11 <- data.frame(cbind(laskar[,1],laskar[,4],laskar[,3],laskar[,2]))
    }  
  if(sol == 2 ) 
    {
      la11 <- data.frame(cbind(laskar[,1],laskar[,5],laskar[,3],laskar[,2]))
    }      
  if(sol == 3 ) 
    {
      la11 <- data.frame(cbind(laskar[,1],laskar[,6],laskar[,3],laskar[,2]))
    }    
  if(sol == 4 ) 
    {
      la11 <- data.frame(cbind(laskar[,1],laskar[,7],laskar[,3],laskar[,2]))
    }  
  if(sol == 5 ) 
    {
      la11 <- data.frame(cbind(laskar[,1],laskar[,8],laskar[,3],laskar[,2]))
    }
    
# replace negative time with postive time
  la11[1] <- la11[1]*-1
# replace preccession angle with sin(precession angle)
  la11[4] <- sin(la11[4])
# replace precession angle with esin(precession angle)
  if(esinw) la11[4] <- la11[2]*la11[4]
# isolate portion of record that is needed
  la11 <- subset(la11, (la11[1] >= tmin) & (la11[1] <= tmax) ) 
# calculate mean values for eccentricity, obliquity, precession
  la11mean <- colMeans(la11)
# calculate std deviations for eccentricity, obliquity, precession
  la11stdev <- sapply(la11,sd)
# standardize eccentricity, tilt, precession
  if(standardize)
    {
     ecc <- ( la11[2]-la11mean[2] ) / la11stdev[2]
     obl <- ( la11[3]-la11mean[3] ) / la11stdev[3]
     prec <- ( la11[4]-la11mean[4] ) / la11stdev[4]
    }
  if(!standardize)
    {
     ecc <- la11[2]
     obl <- la11[3]
     prec <- la11[4]
    }

# now weight each term as desired
  ecc <- ecc * eWt
  obl <- obl * oWt
  prec <- prec * pWt
# combine eccentricity, tilt and precession
  combined <- ecc + obl + prec
  out <- data.frame( cbind(la11[1],combined) )

### then interpolate as needed
#  xout <- seq(la11[1,1],la11[length(la11[,1]),1],by=dt)
  xout <- seq(tmin,tmax,by=dt)
### redefine npts
  npts <- length(xout)
  interp <- approx(out[,1],out[,2],xout,method="linear",n=npts)
  interp <- data.frame(interp)

if (genplot)
  {
### plots
   par(mfrow=c(2,2))
   plot(la11[,1],la11[,2],cex=0.5,xlab="Time (ka BP)",ylab="Eccenticity",main="Eccentricity",type="l")
   plot(la11[,1],la11[,3],cex=0.5,xlab="Time (ka BP)",ylab="Obliquity (radians)",main="Obliquity",type="l")
   if(esinw) plot(la11[,1],la11[,4],cex=0.5,xlab="Time (ka BP)",ylab="Eccentricity*sin(angle)",main="Precession",type="l")
   if(!esinw) plot(la11[,1],la11[,4],cex=0.5,xlab="Time (ka BP)",ylab="sin(angle)",main="Precession",type="l")
   plot(interp,cex=0.25,xlab="Time (ka BP)",ylab="Value",main="ETP");lines(interp,col="red")
  }
  
  return(interp)

### END function etp
}
