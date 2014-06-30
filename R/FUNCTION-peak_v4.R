### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2014 Stephen R. Meyers
###
###########################################################################
### function peak : find maxima of peaks in series, report those that exceed
###                  a threshold value - (SRM: March 1-29, 2012; 
###                                  April 25, 2012; May 22, 2013; May 23, 2013; 
###                                  June 5, 2013; June 14, 2013)
###
###########################################################################

peak <- function (dat,level=NULL,genplot=T,verbose=T) 
{

if(verbose) cat("\n----- FINDING MAXIMA OF PEAKS, FILTERING AT THRESHOLD VALUE -----\n")

dat <- data.frame(dat)
npts <- as.integer(dim(dat)[1])
ncols <- as.integer(dim(dat)[2])
if(verbose) cat(" * Number of data points=", npts,"\n")
if(verbose) cat(" * Number of columns=", ncols,"\n")

if(ncols == 1) { y <- dat[,1] }
if(ncols == 2) { x <- dat[,1]; y <- dat[,2] }

# FORTRAN wrapper
peakID <- function (npts,y) 
 { 
    F_dat = .Fortran('peak_r',PACKAGE='astrochron',
    
    npts=as.integer(npts),y=as.double(y),
    
    loc=integer(as.integer(npts)),iplat=integer(as.integer(npts)),
    numpeak=integer(1),numplat=integer(1)
    )    
# return the results
    return(F_dat)
 }

if(verbose) cat(" * Identifying maxima of peaks\n")
# identify maxima of peaks
res <- peakID(npts,y)

numpeak <- res$numpeak
numplat <- res$numplat
loc <- res$loc[1:numpeak]
iplat <- res$iplat[1:numplat]

if(verbose) cat(" * Number of peaks identified=",numpeak,"\n")

if(numplat>0) 
 {
   if(verbose) cat(" * Number of plateau points detected=",numplat,"\n")
   if(ncols == 1)
    {
      plats = data.frame( cbind(iplat,y[iplat]))
      colnames(plats)[1] <- 'Index'
      colnames(plats)[2] <- 'Plateau_Value'
    }
   
   if(ncols == 2)
    {
      plats = data.frame( cbind(iplat,x[iplat],y[iplat]))
      colnames(plats)[1] <- 'Index'
      colnames(plats)[2] <- 'Plateau_Location'
      colnames(plats)[3] <- 'Plateau_Value'
    }
# this warning should be output regardless of whether verbose selected.
      cat("\n**** WARNING: The following plateaus were not evaluated!:\n")
      print(plats)  
 }
  
if(ncols == 1) { ymax <- y[loc] }
if(ncols == 2) { xloc <- x[loc] ; ymax <- y[loc] }

if( is.null(level) ) 
 {
    if(verbose) cat("\n * No filtering of peaks applied.\n")
    filt_loc <- loc
    if(ncols == 2) filt_xloc <- xloc
    filt_ymax <- ymax
  }    
      
if( !is.null(level) )
 {
   if(verbose) cat(" * Filtering peaks at threshold of",level,"\n")
   filt_loc <- loc[ymax >= level]
   if(ncols == 2) filt_xloc <- xloc[ymax >= level]
   filt_ymax <- ymax[ymax >= level]
   if(verbose) cat(" * Number of peaks >=",level,":", length(filt_ymax),"\n")
 }
  
if(numpeak>0 && ncols == 1) 
  {
        out = data.frame( cbind(filt_loc,filt_ymax) )
        colnames(out)[1] = 'ID'
        colnames(out)[2] = 'Peak_Value'
  }
if(numpeak>0 && ncols == 2)
  {
        out = data.frame( cbind(filt_loc,filt_xloc,filt_ymax) )
        colnames(out)[1] = 'ID'
        colnames(out)[2] = 'Location'
        colnames(out)[3] = 'Peak_Value'
  }
 
if(numpeak>0 && genplot)
 {
   setline = 0.25
   par(mfrow=c(1,1))
   if(ncols == 1)
     {
      plot(1:npts,y,cex=0.5,main="Data with Peak Maxima Identified",xlab="Point Number",ylab="Value",bty="n")
      lines(1:npts,y,col="orange")
      for(j in 1:length(filt_loc))
        {
         abline(v=filt_loc[j],col="blue",lty=22)
         points(filt_loc[j],filt_ymax[j],pch=1,col='blue')
         mtext(filt_loc[j], side=3,line=setline,at=filt_loc[j],cex=0.5,font=4,col="blue")
         setline=setline*-1
        }
      if(numplat > 0)
       {
       for(j in 1:numplat)
        {
         mtext(iplat[j], side=3,line=setline,at=iplat[j],cex=0.5,font=4,col="red")
         points(iplat[j],y[iplat[j]],pch=19,col='red')
         setline=setline*-1
        }
       } 
    }   
   if(ncols == 2)
    {
      plot(dat,cex=0.5,main="Data with Peak Maxima Identified",xlab="Location",ylab="Value",bty="n")
      lines(dat,col="orange")
      for(j in 1:length(filt_loc))
        {
         abline(v=filt_xloc[j],col="blue",lty=22)
         points(filt_xloc[j],filt_ymax[j],pch=1,col='blue')
         mtext(filt_xloc[j], side=3,line=setline,at=filt_xloc[j],cex=0.5,font=4,col="blue")
         setline=setline*-1
        }
      if(numplat > 0)
       {
       for(j in 1:numplat)
        {
         mtext(iplat[j], side=3,line=setline,at=x[iplat[j]],cex=0.5,font=4,col="red")
         points(x[iplat[j]],y[iplat[j]],pch=19,col='red')
         setline=setline*-1
        }
       } 
    }   
# end genplot section
 }
 
if(numpeak>0) return( out )

### END function peak
}

