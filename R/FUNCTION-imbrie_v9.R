### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2024 Stephen R. Meyers
###
###########################################################################
### imbrie: Imbrie and Imbrie (1980) ice model. Model follows convention
###         used in Analyseries.
###          (SRM: April 14-16, 2015; June 19, 2018; December 5-6, 2018; 
###                January 7, 2019; December 3, 2020; January 14, 2021;
###                September 21, 2022; November 7, 2022; December 6, 2024)
###
###########################################################################

imbrie <- function (insolation=NULL,Tm=17,b=0.6,times=NULL,initial=0,burnin=100,standardize=T,output=T,genplot=1,check=T,verbose=T)
{

  if(verbose) cat("\n----- GENERATING IMBRIE and IMBRIE (1980) ICE MODEL -----\n")
      
  if(is.null(insolation))
   {
     if(verbose) cat(" * Using default insolation series: 65 deg. North, 21 June\n" )
     insolation=getLaskar("insolation",verbose=F)
     insolation=iso(insolation,xmin=0,xmax=1100,genplot=F,verbose=F)
    }  

  if(!is.null(insolation)) insolation=data.frame(insolation)
  npts=length(insolation[,1])

  if(check && is.null(times)) 
   {
     if(length(Tm) != length(b) || length(Tm) != 1)
      {
        cat("\n**** ERROR: vectors Tm and b do not have matching entries\n")
        stop("    TERMINATING NOW!")
      }
   }
   
   if(check && !is.null(times))
    {
      if(length(Tm) != length(times) || length(b) != length(times))
       {
        cat("\n**** ERROR: vectors Tm, b and times do not have matching entries\n")
        stop("    TERMINATING NOW!")
      }
      if(insolation[npts,1] > max(times))
        {
        cat("\n**** ERROR: maximum time for model (in vector times) must be >= maximum time in insolation\n")
        stop("    TERMINATING NOW!")
      }
   }     

# flip so we are going from past to present
  insolation <- insolation[order(insolation[,1], na.last = NA, decreasing = T),]
  insolation[1] <- insolation[1]*-1  

# center the insolation to mean of zero
  insolation[2] <- insolation[2]-colMeans(insolation[2])

# follow convention of Analyseries  
  x<-insolation[,2]*-1
  
  if(is.null(times))
   {
     b2=b
     Tm2=Tm
     ipts=1
   }
     
  if(!is.null(times))
   {
     ice=data.frame(cbind(times,Tm,b))    
# also flip sequence of ice model parameters to match insolation direction
     ice <- ice[order(ice[,1], na.last = NA, decreasing = T),]
     ice[1] <- ice[1]*-1
     getTime <- double(npts)
# assign values to each model step
     for(i in 1:npts) getTime[i]=max(which( ice[1] <= insolation[i,1] ))
     b2=ice[getTime,3]
     Tm2=ice[getTime,2]
     ipts=npts
   }     

runMod <- function (npts,x,ipts,Tm2,b2,yinit)
 {
    F_dat = .Fortran( 'imbrie_r',
    
    npts=as.integer(npts),x=as.double(x),ipts=as.integer(ipts),
    Tm2=as.double(Tm2),b2=as.double(b2),yinit=as.double(yinit),
    
    y=double(npts)
    )

# return the results
    return(F_dat)
 }

  resMod=runMod(npts=npts,x=x,ipts=ipts,Tm2=Tm2,b2=b2,yinit=initial)

# remove burn in points, and make a data frame
  out<-data.frame(cbind(insolation[burnin:npts,1],resMod$y[burnin:npts]))

# flip to original direction
  out <- out[order(out[,1], na.last = NA, decreasing = T),]
  out[1] <- out[1]*-1  

  if(standardize) 
   {
     out[2] <- out[2]-min(out[2])
     out[2] <- out[2]/max(out[2])
   }  

if(genplot==1 || genplot==2 || genplot==3)
  {
# remove burn in points, and make a data frame
  insolation<-data.frame(cbind(insolation[burnin:npts,1],insolation[burnin:npts,2]))

# flip to original direction
  insolation <- insolation[order(insolation[,1], na.last = NA, decreasing = T),]
  insolation[1] <- insolation[1]*-1  
  
### plots
# 2 plots, no animation
   if(genplot==1) 
    {
     par(mfrow=c(2,1))
     plot(insolation,cex=0.5,xlab="Time (ka BP)",ylab="Insolation",main="Insolation model",type="l",col="red")
#    plot(out,cex=0.5,xlab="Time (ka BP)",ylab="Ice volume",main="Ice volume model",type="l",ylim=c(max(out[,2]),min(out[,2])))
     plot(out,cex=0.5,xlab="Time (ka BP)",ylab="Ice volume",main="Ice volume model",type="l",ylim=c(min(out[,2]),max(out[,2])))
     polygon(c(min(out[,1]),out[,1],max(out[,1])),c(min(out[,2]),out[,2],min(out[,2])),col="cornflowerblue",border=NA)
   }

# three plots with animation
   if(genplot==2) 
    {
      dev.new(height=5.5,width=10)
      layout(matrix(c(1,2,3,3), 2, 2, byrow = FALSE))
      pts2=length(out[,2])
      xlim=c(min(insolation[1:pts2,2]),max(insolation[1:pts2,2]))
      ylim=c(min(out[1:pts2,2]),max(out[1:pts2,2]))
      for(i in pts2:1)
       {
         plot(insolation[i:pts2,1],insolation[i:pts2,2],cex=0.5,xlab="Time (ka BP)",ylab="Insolation",main="Insolation model",font.lab=2,type="l",ylim=xlim,col="red",lwd=2)
         plot(out[i:pts2,1],out[i:pts2,2],cex=0.5,xlab="Time (ka BP)",ylab="Ice volume",main="Ice volume model",font.lab=2,type="l",ylim=ylim,col="cornflowerblue",lwd=2)
         plot(insolation[i:pts2,2],out[i:pts2,2],type="l",xlab="Insolation",ylab="Ice volume",font.lab=2,col="darkgreen",lty=2,xlim=xlim,ylim=ylim,lwd=2)
         points(insolation[i,2],out[i,2],col="black",pch=16)
         Sys.sleep(0.15)
       }
    }    

# three plots, no animation
   if(genplot==3) 
    {
      dev.new(height=5.5,width=10)
      layout(matrix(c(1,2,3,3), 2, 2, byrow = FALSE))
      plot(insolation,cex=0.5,xlab="Time (ka BP)",ylab="Insolation",main="Insolation model",type="l",col="red",lwd=2)
      plot(out,cex=0.5,xlab="Time (ka BP)",ylab="Ice volume",main="Ice volume model",type="l",col="cornflowerblue",lwd=2)
      polygon(c(min(out[,1]),out[,1],max(out[,1])),c(min(out[,2]),out[,2],min(out[,2])),col="cornflowerblue",border=NA)
      plot(insolation[,2],out[,2],type="l",xlab="Insolation",ylab="Ice volume",col="darkgreen",lty=2,lwd=2)
      points(insolation[,2],out[,2],col="darkgreen",pch=16)
    }    
    
# END all genplot
  }
  
  if(output) return(out)

### END function imbrie
}
